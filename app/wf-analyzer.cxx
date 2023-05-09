#include <iostream>
#include <sstream>
#include <cstdlib>
#include <stdlib.h> /* atoi, atof */
#include <map>

#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>

#define NCH 6
int NCHmax = NCH;

std::map<int,std::map<int,double> > ADCrangeMap;
std::map<int,std::map<int,double> > ADCoffsetMap;

int verbose = 0;

double ADC2voltage(int ADC, int range, int offset){
    double voltage = ADC*range/127.-offset;
    return voltage;
}

double voltage2ADC(double voltage, int range, int offset){
    int ADC = std::round((voltage+offset)*127/range);
    return ADC;
}

int getWidth(const std::vector<float> * pwf, double ped, int left, int right, double thr, double max){
    if (verbose>10) std::cout<<"Searching for peak width over "<<thr<<", left "<<left<<" right "<<right<<" max "<<max<<std::endl;
    if (thr>=max) return 0;
    int l = left;
    int r = right;
    if (-(pwf->at(l)-ped)<thr){
        if (verbose>10) std::cout<<" pwf[left] "<<-(pwf->at(l)-ped)<<" is lower than threshold "<<thr<<std::endl;
        while (l+1<=right&&-(pwf->at(l)-ped)<thr){
            if (verbose>10) std::cout<<"   pwf["<<l<<"] "<<-(pwf->at(l)-ped)<<" < "<<thr<<", l++"<<std::endl;
            l++;
        }
    }
    else{
        if (verbose>10) std::cout<<" pwf[left] "<<-(pwf->at(l)-ped)<<" is higher than threshold "<<thr<<std::endl;
        while (l-1>=0&&pwf->at(l-1)>thr){
            if (verbose>10) std::cout<<"   pwf["<<l-1<<"] "<<pwf->at(l-1)<<" > "<<thr<<", l++"<<std::endl;
            l--;
        }
    }
    if (-(pwf->at(r)-ped)<thr){
        if (verbose>10) std::cout<<" pwf[right] "<<-(pwf->at(r)-ped)<<" is lower than threshold "<<thr<<std::endl;
        while (r-1>=l&&-(pwf->at(r)-ped)<thr){
            if (verbose>10) std::cout<<"   pwf["<<r<<"] "<<-(pwf->at(r)-ped)<<" < "<<thr<<", r--"<<std::endl;
            r--;
        }
    }
    else{
        if (verbose>10) std::cout<<" pwf[right] "<<-(pwf->at(r)-ped)<<" is higher than threshold "<<thr<<std::endl;
        while (r+1<pwf->size()&&pwf->at(r+1)>thr){
            if (verbose>10) std::cout<<"   pwf["<<r+1<<"] "<<pwf->at(r+1)<<" > "<<thr<<", r++"<<std::endl;
            r++;
        }
    }
    if (verbose>10){
        std::cout<<" width is "<<r<<" - "<<l<<" + 1= "<<r-l+1<<std::endl;
        std::cout<<" "<<(l-1>=0?l-1:0)<<" | "<<l<<" "<<r<<" | "<<(r+1>=pwf->size()?pwf->size()-1:r+1)<<std::endl;
        std::cout<<" "<<-(pwf->at((l-1>=0?l-1:0))-ped)<<" | "<<-(pwf->at(l)-ped)<<" "<<-(pwf->at(r)-ped)<<" | "<<-(pwf->at((r+1>=pwf->size()?pwf->size()-1:r+1))-ped)<<std::endl;
    }
    return (r-l+1);
};

bool loadCalibADC(const char* filename){
    TChain * ichain = new TChain("t");
    ichain->Add(filename);
    if (ichain->GetEntries()==0) return true;
    int range[NCH];
    int offset[NCH];
    int run;
    ichain->SetBranchAddress("run",&run);
    for (int ch = 0; ch<NCH; ch++){
        ichain->SetBranchAddress(Form("range%c",'A'+ch),&range[ch]);
        ichain->SetBranchAddress(Form("offset%c",'A'+ch),&offset[ch]);
    }
    for (int e = 0; e<ichain->GetEntries(); e++){
        ichain->GetEntry(e);
        for (int ch = 0; ch<NCH; ch++){
            ADCrangeMap[run][ch] = range[ch];
            ADCoffsetMap[run][ch] = offset[ch];
        }
    }
    return false;
};

void print_usage(char * prog_name){
    fprintf(stderr,"Usage %s [options] runname\n",prog_name);
    fprintf(stderr,"[options]\n");
};

int main(int argc, char** argv)
{

    int nSamplesBeforeSignal = 300;
    double                ps_dt = 0.8; // ns
    int                   ps_size = 12500;
    bool                  use3rmsAsThrshold = false;
    double                threshold[NCH] = {0};
    double                resistance = 50;// Ohm
    double   ADCrange[NCH] = {1000,1000,1000,1000,1000,1000};
    double   ADCoffset[NCH] = {800,800,800,800,800,800};
    double   ADCstep = 0; // will be set according to the largest step calculated by ADCrange and ADCoffset
    int      ADCmax = 26;
    int      ADCmin = -127;
    TString runname = "test";
    TString comment = "";
    bool debug = false;
    int tZoomLeft = 0;
    int tZoomRight = 500;
    bool isOldData = false;

    // Load options
    int    opt_result;
    std::stringstream stream;
    int run_begin = 0;
    int run_end = -1; // keep run_end < run_begin if you don't want the following codes to automatically load runs
    while((opt_result=getopt(argc,argv,"OA:C:T:b:e:v:r:c:z:dh"))!=-1){
        stream.clear();
        if(optarg) stream.str(optarg);
        switch(opt_result){
            /* INPUTS */
            case 'O':
                isOldData = true; // Februaray data with only 4 channels for picoscope: no pre-T0 yet
                break;
            case 'C':
                std::cerr<<"Loading ADCrange and ADCoffset from the calibration file \""<<optarg<<"\""<<std::endl;
                if (loadCalibADC(optarg)){
                    std::cerr<<"ERROR: didn't find any entry in this calibration file!"<<std::endl;
                }
                break;
            case 'T':
                if (stream.str()=="3rms"){
                    use3rmsAsThrshold = true;
                    std::cerr<<"Thresholds set to 3xRMS in each event"<<std::endl;
                }
                else{
                    stream>>threshold[0]>>threshold[1]>>threshold[2]>>threshold[3]>>threshold[4]>>threshold[5];
                    std::cerr<<"Thresholds set to "<<threshold[0]<<" "<<threshold[1]<<" "<<threshold[2]<<" "<<threshold[3]<<threshold[4]<<" "<<threshold[5]<<" mV"<<std::endl;
                }
                break;
            case 'z':
                stream>>tZoomLeft>>tZoomRight;
                std::cerr<<"Time zoom range to "<<tZoomLeft<<" ~ "<<tZoomRight<<std::endl;
                break;
            case 'A':
                stream>>ADCmin>>ADCmax;
                std::cerr<<"ADC range set to "<<ADCmin<<" ~ "<<ADCmax<<std::endl;
                break;
            case 'b':
                run_begin = atoi(optarg);
                std::cerr<<"begin of the run number to "<<run_begin<<std::endl;
                break;
            case 'e':
                run_end = atoi(optarg);
                std::cerr<<"end of the run number to "<<run_end<<std::endl;
                break;
            case 'd':
                debug = true;
                std::cerr<<"This is debugging mode. Will generate multiple width vectors counting over 10~100 mV for each channel"<<std::endl;
                break;
            case 'v':
                verbose = atoi(optarg);
                std::cerr<<"verbose level set to "<<verbose<<std::endl;
                break;
            case 'r':
                runname = optarg;
                std::cerr<<"Run name \""<<runname<<"\""<<std::endl;
                break;
            case 'c':
                comment = optarg;
                std::cerr<<"comment \""<<comment<<"\""<<std::endl;
                break;
            case '?':
                fprintf(stderr,"Wrong option! optopt=%c, optarg=%s\n", optopt, optarg);
            case 'h':
            default:
                print_usage(argv[0]);
                return 1;
        }
    }
    if (isOldData) NCHmax = 4;

    TChain* chain = new TChain("tr");
    for (int i = optind; i<argc; i++){
        chain->Add(argv[i]);
    }
    for (int run = run_begin; run<=run_end; run++){
        chain->Add(Form("decdata/rawdata%05d.root",run));
    }
    if(0 == chain->GetNtrees()) return 1;

    // Set Branch Addresses
    Int_t run = 0;
    Int_t event = 0;
    std::vector< std::vector<Float_t>* > pwfs(NCHmax, nullptr);
    // muon beam monitor
    std::vector<unsigned char> * mbm_x_ch = 0;
    std::vector<double>        * mbm_x_time = 0;
    std::vector<unsigned char> * mbm_y_ch = 0;
    std::vector<double>        * mbm_y_time = 0;

    // load the first run number
    chain->SetBranchAddress("run",&run);
    for (Long64_t e = 0; e<chain->GetEntries(); e++){
        chain->GetEntry(e);
        bool good = true;
        for (int ch = 0; ch<NCHmax; ch++){
            if (ADCrangeMap[run][ch]==0) good = false;
        }
        if (!good) continue;
        double adcstep = 0;
        for (int ch = 0; ch<NCHmax; ch++){
            ADCrange[ch] = ADCrangeMap[run][ch];
            ADCoffset[ch] = ADCoffsetMap[run][ch];
            adcstep = ADCrange[ch]/127.;
            if (adcstep>ADCstep) ADCstep = adcstep;
            std::cerr<<"  "<<ch<<" "<<ADCrange[ch]<<" "<<ADCoffset[ch]<<std::endl;
        }
        break;
    }

    chain->SetBranchAddress("run", &run);
    chain->SetBranchAddress("event", &event);
    for(int ch=0; ch<NCHmax; ch++){
        chain->SetBranchAddress(Form("pwf%c", 'A'+ch), &pwfs[ch]);
    }
    bool hasMBM = true;
    if (!chain->GetListOfBranches()->FindObject("mbm.x.ch")) hasMBM = false;
    if (!chain->GetListOfBranches()->FindObject("mbm.x.time")) hasMBM = false;
    if (!chain->GetListOfBranches()->FindObject("mbm.y.ch")) hasMBM = false;
    if (!chain->GetListOfBranches()->FindObject("mbm.y.time")) hasMBM = false;
    if (hasMBM){
        chain->SetBranchAddress("mbm.x.ch",&mbm_x_ch);
        chain->SetBranchAddress("mbm.x.time",&mbm_x_time);
        chain->SetBranchAddress("mbm.y.ch",&mbm_y_ch);
        chain->SetBranchAddress("mbm.y.time",&mbm_y_time);
    }

    gStyle->SetOptStat(0);
    TCanvas* canvas = new TCanvas("canv","canv");
    if (isOldData){
        canvas->SetCanvasSize(1080,1200);
    }
    else{
        canvas->SetCanvasSize(1620,1200);
    }
    TPad* pads[NCH];
    canvas->cd();
    for(int i=0; i<NCHmax; i++){
        if (isOldData){
            pads[i] = new TPad(Form("p%d", NCHmax+i), "",
                    i/2 * 0.5,    (1-i%2)*0.45,
                    (i/2+1)*0.5,  (2-i%2)*0.45);
        }
        else{ // FIXME: shall change to 6 pads mode
            pads[i] = new TPad(Form("p%d", NCHmax+i), "",
                    i/2 * 1./3,    (1-i%2)*0.45,
                    (i/2+1)*1./3,  (2-i%2)*0.45);
        }
        pads[i]->Draw();
        pads[i]->SetGridx(1);
        pads[i]->SetGridy(1);
    }
    TLatex * text = new TLatex(0.02,0.92,"");
    text->SetTextFont(42);
    text->SetTextSize(0.023);
    text->Draw();

    TFile * ofile = new TFile(Form("anaout/%s.hit.root",runname.Data()),"recreate");
    TTree * otree = new TTree("t","t");
    TH2D * hwf[NCH];
    TH1D * hrms[NCH];
    TH1D * hped[NCH];
    TH1D * harea[NCH];
    TH1D * hheight[NCH];
    TH2D * hheightVStime[NCH];
    TH2D * hheightVStimeZoom[NCH];
    Int_t wfCount[NCH] = {0};
    Int_t peakCount[NCH] = {0};
    double rms[NCH] = {0};
    double pedestal[NCH] = {0};
    int    nPeaks[NCH] = {0};
    std::vector<double> * start[NCH] = {0}; // ns
    std::vector<int> * width[NCH] = {0}; // #samples
    std::vector<int> * widthOverCuts[NCH][10] = {0}; // #samples
    std::vector<int> * fwhm[NCH] = {0}; // #samples
    std::vector<double> * height[NCH] = {0}; // mV over pedestal
    std::vector<double> * area[NCH] = {0}; // pC over pedestal
    otree->Branch("run",&run);
    otree->Branch("event",&event);
    for (int ch = 0; ch<NCHmax; ch++){
        otree->Branch(Form("rms%c",'A'+ch),&rms[ch]);
        otree->Branch(Form("ped%c",'A'+ch),&pedestal[ch]);
        otree->Branch(Form("nPeaks%c",'A'+ch),&nPeaks[ch]);
        otree->Branch(Form("start%c",'A'+ch),&start[ch]);
        otree->Branch(Form("width%c",'A'+ch),&width[ch]);
        if (debug){
            for (int i = 0; i<10; i++){
                otree->Branch(Form("width%c%d",'A'+ch,(i+1)*10),&widthOverCuts[ch][i]);
                widthOverCuts[ch][i] = new std::vector<int>;//(1000); // TODO: this 1000 is to keep storage space for each vector. Might not be enough in very noisy runns
            }
        }
        otree->Branch(Form("fwhm%c",'A'+ch),&fwhm[ch]);
        otree->Branch(Form("height%c",'A'+ch),&height[ch]);
        otree->Branch(Form("area%c",'A'+ch),&area[ch]);
        start[ch] = new std::vector<double>;
        width[ch] = new std::vector<int>;
        fwhm[ch] = new std::vector<int>;
        height[ch] = new std::vector<double>;
        area[ch] = new std::vector<double>;
        hwf[ch] = new TH2D(Form("hwf%d",ch),Form("PicoScope ch%c;Time [ns];Voltage [mV]",'A'+ch),ps_size,-0.5*ps_dt,(ps_size-0.5)*ps_dt,ADCmax-ADCmin+1,(ADCmin-0.5)*ADCstep,(ADCmax+0.5)*ADCstep);
        hrms[ch] = new TH1D(Form("hrms%d",ch),Form("PicoScope ch%c;RMS in the first %d samples [mV];Count",'A'+ch,nSamplesBeforeSignal),128,0,20);
        hped[ch] = new TH1D(Form("hped%d",ch),Form("PicoScope ch%c;Pedestal in the first %d samples [mV];Count",'A'+ch,nSamplesBeforeSignal),256,-20,20);
        harea[ch] = new TH1D(Form("harea%d",ch),Form("PicoScope ch%c;Peak area over pedestal [pC];Count",'A'+ch),800,0,800);
        hheight[ch] = new TH1D(Form("hheight%d",ch),Form("PicoScope ch%c;Peak height over pedestal [mV];Count",'A'+ch),0-ADCmin,0,(0-ADCmin)*ADCstep);
        hheightVStime[ch] = new TH2D(Form("hheightVStime%d",ch),Form("PicoScope ch%c;Start time [ns];Peak height over pedestal [mV]",'A'+ch),ps_size/10,0,ps_size*ps_dt,0-ADCmin,0,(0-ADCmin)*ADCstep); // 10 pC, 10 TDC (8 ns)
        hheightVStimeZoom[ch] = new TH2D(Form("hheightVStimeZoom%d",ch),Form("PicoScope ch%c;Start time [ns];Peak height over pedestal [mV]",'A'+ch),(tZoomRight-tZoomLeft)/ps_dt+1,tZoomLeft-0.5*ps_dt,tZoomRight+0.5*ps_dt,0-ADCmin,0,(0-ADCmin)*ADCstep); // 10 mV, 1 TDC (0.8 ns)
    }
    Int_t wfPeakCount[NCH][100] = {0};
    TH2D * hwfPeak[NCH][100]; // Group peaks into 100 regions by height: 10 mV per region, from 0 to 1000 mV
    for (int ch = 0; ch<NCHmax; ch++){
        for (int i = 0; i<100; i++){
            hwfPeak[ch][i] = new TH2D(Form("hwfPeak%c%d",'A'+ch,i),";Time [ns];Voltage [mV]",151,-0.5*ps_dt,(150+0.5)*ps_dt,0-ADCmin,0,(0-ADCmin)*ADCstep);
        }
    }
    if (hasMBM){
        otree->Branch("mbm.x.ch",&mbm_x_ch);
        otree->Branch("mbm.x.time",&mbm_x_time);
        otree->Branch("mbm.y.ch",&mbm_y_ch);
        otree->Branch("mbm.y.time",&mbm_y_time);
    }

    Long64_t nEntries = chain->GetEntries();
    int run_min = 0;
    int run_max = 0;
    std::vector<double> wfpeak;
    TH1I * hADC = new TH1I("hADC","",255,-127.5,127.5);
    int prev_run = -1;
    for(Long64_t entry=0; entry<nEntries; entry++){
        chain->GetEntry(entry);
        if (entry==0) {
            run_max = run;
            run_min = run;
        }
        if (run!=prev_run){
            for (int ch = 0; ch<NCHmax; ch++){
                if (fabs(ADCrange[ch] - ADCrangeMap[run][ch])>1e-10||fabs(ADCoffset[ch] - ADCoffsetMap[run][ch])>1e-10){
                    std::cerr<<"WARNING! ADC parameters changed in run:"<<run<<" ch "<<ch
                        <<" initial: "<<ADCrange[ch]<<" "<<ADCoffset[ch]
                        <<" new: "<<ADCrangeMap[run][ch]<<" "<<ADCoffsetMap[run][ch]
                        <<", plots axis with mV unit will have improper binning for this run!"<<std::endl;
                    ADCrange[ch] = ADCrangeMap[run][ch];
                    ADCoffset[ch] = ADCoffsetMap[run][ch];
                }
            }
        }
        prev_run = run;
        if (run_min>run) run_min = run;
        if (run_max<run) run_max = run;

        for(int ch=0; ch<NCHmax; ch++){
            rms[ch] = 0;
            pedestal[ch] = 0;
            nPeaks[ch] = 0;
            start[ch]->clear();
            width[ch]->clear();
            if (debug){
                for (int i = 0; i<10; i++){
                    widthOverCuts[ch][i]->clear();
                }
            }
            fwhm[ch]->clear();
            height[ch]->clear();
            area[ch]->clear();
            if(pwfs[ch]){
                wfCount[ch]++;
                int nSamplesCounted = 0;
                hADC->Reset();
                for(size_t i=0; i<nSamplesBeforeSignal&&i<pwfs[ch]->size(); i++){
                    hADC->Fill(voltage2ADC(pwfs[ch]->at(i),ADCrange[ch],ADCoffset[ch]));
                }
                // first get pedestal
                double mpv = ADC2voltage(hADC->GetBinCenter(hADC->GetMaximumBin()),ADCrange[ch],ADCoffset[ch]);
                for(size_t i=0; i<nSamplesBeforeSignal&&i<pwfs[ch]->size(); i++){
                    if (fabs(pwfs[ch]->at(i) - mpv)>20) continue; // 20 mV is a hard coded value to reject peaks in the pedestal
                    rms[ch]+=pow(pwfs[ch]->at(i),2);
                    pedestal[ch]+=pwfs[ch]->at(i);
                    nSamplesCounted++;
                }
                pedestal[ch]/=nSamplesCounted;
                // then get rms
                nSamplesCounted=0;
                for(size_t i=0; i<nSamplesBeforeSignal&&i<pwfs[ch]->size(); i++){
                    rms[ch]+=pow(pwfs[ch]->at(i)-pedestal[ch],2);
                    nSamplesCounted++;
                }
                rms[ch] /= nSamplesCounted;
                rms[ch] = sqrt(rms[ch]);

                // tune the threshold in case of 3rms option
                if (use3rmsAsThrshold){
                    threshold[ch] = rms[ch]*3;
                }

                hrms[ch]->Fill(rms[ch]);
                hped[ch]->Fill(pedestal[ch]);

                bool inPeak = false;
                double cur_start = 0;
                double cur_width = 0;
                double cur_area = 0;
                double cur_height = 0;
                wfpeak.clear();
                for(size_t i=0; i<pwfs[ch]->size(); i++){
                    hwf[ch]->Fill(i*ps_dt, pwfs[ch]->at(i));
                    double voltage = -(pwfs[ch]->at(i) - pedestal[ch]);
                    if (verbose>5) std::cout<<i<<": pwfs = "<<pwfs[ch]->at(i)<<" mV, voltage = "<<voltage<<" mV"<<std::endl;
                    if (voltage>threshold[ch]){
                        if (verbose>5) std::cout<<"above threshold!";
                        if (!inPeak){
                            cur_start = i*ps_dt;
                            if (verbose>5) std::cout<<" new peak "<<nPeaks[ch]<<std::endl;
                        }
                        else{
                            if (verbose>5) std::cout<<" already in peak "<<nPeaks[ch]<<std::endl;
                        }
                        inPeak = true;
                        cur_width++;
                        if (cur_height<voltage) cur_height = voltage;
                        cur_area += voltage*1e-3*ps_dt*1e-9/resistance/1e-12; // V*second/Ohm/1e-12 -> pC
                        wfpeak.push_back(voltage);
                    }
                    else{
                        if (inPeak){
                            if (verbose>5) std::cout<<" peak "<<nPeaks[ch]<<" stop! start = "<<cur_start<<" height = "<<pedestal[ch]-cur_height<<std::endl;
                            nPeaks[ch]++;
                            // rewind to find fwhm
                            start[ch]->push_back(cur_start);
                            width[ch]->push_back(cur_width);
                            if (debug){
                                for (int iThr = 0; iThr<10; iThr++){
                                    widthOverCuts[ch][iThr]->push_back(getWidth(pwfs[ch],pedestal[ch],i-cur_width,i,(iThr+1)*10,cur_height));
                                }
                            }
                            fwhm[ch]->push_back(getWidth(pwfs[ch],pedestal[ch],i-cur_width,i,cur_height/2,cur_height));
                            height[ch]->push_back(cur_height);
                            area[ch]->push_back(cur_area);
                            harea[ch]->Fill(cur_area);
                            hheight[ch]->Fill(cur_height);
                            hheightVStime[ch]->Fill(cur_start,cur_height);
                            hheightVStimeZoom[ch]->Fill(cur_start,cur_height);
                            peakCount[ch]++;
                            int groupID = cur_height/10;
                            if (groupID>=100) groupID = 99;
                            for (int i = 0; i<wfpeak.size(); i++){
                                hwfPeak[ch][groupID]->Fill(i*ps_dt,wfpeak[i]);
                            }
                            wfPeakCount[ch][groupID]++;
                            cur_start = 0;
                            cur_width = 0;
                            cur_height = 0;
                            cur_area = 0;
                            wfpeak.clear();
                        }
                        inPeak = false;
                    }
                }
            }
        }
        otree->Fill();
    }

    text->SetText(0.02,0.92,Form("%s, %lld Entries", comment.Data(), nEntries));
    canvas->cd();
    text->Draw();

    for(int ch=0; ch<NCHmax; ch++){
        pads[ch]->cd();
        pads[ch]->SetLogz(1);
        hwf[ch]->SetTitle(Form("PicoScope ch%c: %d waveforms",'A'+ch,wfCount[ch]));
        hwf[ch]->GetYaxis()->SetMaxDigits(2);
        hwf[ch]->Draw("COLZ");
        hwf[ch]->Write();
    }
    canvas->SaveAs(Form("pictures/00wf.%s.png",runname.Data()));
    //
    for(int ch=0; ch<NCHmax; ch++){
        pads[ch]->cd();
        pads[ch]->SetLogz(1);
        hwf[ch]->GetXaxis()->SetRangeUser(tZoomLeft,tZoomRight);
        hwf[ch]->Draw("COLZ");
        hwf[ch]->Write();
    }
    canvas->SaveAs(Form("pictures/01wfzoom.%s.png",runname.Data()));
    //
    for(int ch=0; ch<NCHmax; ch++){
        pads[ch]->cd();
        pads[ch]->SetLogy(1);
        hped[ch]->GetYaxis()->SetRangeUser(0.5,nEntries*1.1);
        hped[ch]->SetTitle(Form("PicoScope ch%c, pedestal in the first %d samples: %.1f #pm %.1f mV",'A'+ch,nSamplesBeforeSignal,hped[ch]->GetBinCenter(hped[ch]->GetMaximumBin()),hped[ch]->GetRMS()));
        hped[ch]->Draw();
        hped[ch]->Write();
    }
    canvas->SaveAs(Form("pictures/02ped.%s.png",runname.Data()));
    //
    for(int ch=0; ch<NCHmax; ch++){
        pads[ch]->cd();
        pads[ch]->SetLogy(1);
        hrms[ch]->GetYaxis()->SetRangeUser(0.5,nEntries*1.1);
        hrms[ch]->SetTitle(Form("PicoScope ch%c, RMS in the first %d samples: %.1f #pm %.1f mV",'A'+ch,nSamplesBeforeSignal,hrms[ch]->GetBinCenter(hrms[ch]->GetMaximumBin()),hrms[ch]->GetRMS()));
        hrms[ch]->Draw();
        hrms[ch]->Write();
    }
    canvas->SaveAs(Form("pictures/03rms.%s.png",runname.Data()));
    //
    for (int group = 5; group<=20; group*=2){
        double height_min = group*10; // mV
        double height_max = height_min+10;
        for(int ch=0; ch<NCHmax; ch++){
            pads[ch]->cd();
            pads[ch]->SetLogy(0);
            pads[ch]->SetLogz(1);
            hwfPeak[ch][group]->SetTitle(Form("Picoscope ch%c, %d peaks with height from %.0f mV to %.0f mV out of %d waveforms",'A'+ch,wfPeakCount[ch][group],height_min,height_max,wfCount[ch]));
            hwfPeak[ch][group]->GetYaxis()->SetRangeUser(0,250);
            hwfPeak[ch][group]->GetXaxis()->SetRangeUser(-0.4,60.4);
            hwfPeak[ch][group]->Draw("COLZ");
            hwfPeak[ch][group]->Write();
        }
        canvas->SaveAs(Form("pictures/04wfpeak%03.0fmV.%s.png",height_min,runname.Data()));
    }
    //
    for(int ch=0; ch<NCHmax; ch++){
        pads[ch]->cd();
        pads[ch]->SetLogy(1);
        harea[ch]->SetTitle(Form("PicoScope ch%c: %d waveforms, %d peaks (over %.0f mV)",'A'+ch,wfCount[ch],peakCount[ch],threshold[ch]));
        harea[ch]->Draw();
        harea[ch]->Write();
    }
    canvas->SaveAs(Form("pictures/06area.%s.png",runname.Data()));
    //
    for(int ch=0; ch<NCHmax; ch++){
        pads[ch]->cd();
        harea[ch]->GetXaxis()->SetRangeUser(0,100);
        harea[ch]->GetYaxis()->SetRangeUser(0.5,harea[ch]->GetMaximum()*2);
        harea[ch]->Draw();
        harea[ch]->Write();
    }
    canvas->SaveAs(Form("pictures/07areaZoom.%s.png",runname.Data()));
    //
    for(int ch=0; ch<NCHmax; ch++){
        pads[ch]->cd();
        pads[ch]->SetLogy(1);
        hheight[ch]->SetTitle(Form("PicoScope ch%c: %d waveforms, %d peaks (over %.0f mV)",'A'+ch,wfCount[ch],peakCount[ch],threshold[ch]));
        hheight[ch]->Draw();
        hheight[ch]->Write();
    }
    canvas->SaveAs(Form("pictures/08height.%s.png",runname.Data()));
    //
    for(int ch=0; ch<NCHmax; ch++){
        pads[ch]->cd();
        pads[ch]->SetLogy(0);
        hheight[ch]->GetXaxis()->SetRangeUser(0,100);
        hheight[ch]->GetYaxis()->SetRangeUser(0.5,hheight[ch]->GetMaximum()*2);
        hheight[ch]->Draw();
        hheight[ch]->Write();
    }
    canvas->SaveAs(Form("pictures/09heightZoom.%s.png",runname.Data()));
    //
    for(int ch=0; ch<NCHmax; ch++){
        pads[ch]->cd();
        pads[ch]->SetLogy(0);
        pads[ch]->SetLogz(1);
        hheightVStime[ch]->SetTitle(Form("PicoScope ch%c: %d waveforms, %d peaks (over %.0f mV)",'A'+ch,wfCount[ch],peakCount[ch],threshold[ch]));
        hheightVStime[ch]->Draw("COLZ");
        hheightVStime[ch]->Write();
    }
    canvas->SaveAs(Form("pictures/12heightVStime.%s.png",runname.Data()));
    //
    for(int ch=0; ch<NCHmax; ch++){
        pads[ch]->cd();
        pads[ch]->SetLogy(0);
        pads[ch]->SetLogz(1);
        hheightVStimeZoom[ch]->SetTitle(Form("PicoScope ch%c: %d waveforms, %d peaks (over %.0f mV)",'A'+ch,wfCount[ch],peakCount[ch],threshold[ch]));
        hheightVStimeZoom[ch]->Draw("COLZ");
        hheightVStimeZoom[ch]->Write();
    }
    canvas->SaveAs(Form("pictures/13heightVStimeZoom.%s.png",runname.Data()));
    //

    for (int ch = 0; ch<NCHmax; ch++){
        for (int i = 0; i<100; i++){
            hwfPeak[ch][i]->Write();
        }
    }

    otree->Write();
    ofile->Close();
};
