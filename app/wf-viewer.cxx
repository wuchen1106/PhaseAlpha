#include <iostream>
#include <cstdlib>
#include <stdlib.h> /* atoi, atof */

#include <TVectorF.h>
#include <TParameter.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TApplication.h>

void print_usage(char * prog_name){
    fprintf(stderr,"Usage %s [options] runname\n",prog_name);
    fprintf(stderr,"[options]\n");
};

int main(int argc, char** argv)
{

    bool flip = false;
    int iEntryStart = 0;
    int nEntriesToScan = 0;
    TString runname = "test";
    TString fOutputFile = "";
    double zoomTimeLeft = 0;
    double zoomTimeRight = 0;
    bool batchMode = false;
    TChain * ichain_peak = 0;
    double t0[4] = {0,1.68,1.45,8.13};
    double hmin[4] = {35,35,25,25};

    // Load options
    int    opt_result;
    std::stringstream stream;
    while((opt_result=getopt(argc,argv,"F:Z:r:s:n:o:p:bh"))!=-1){
        stream.clear();
        if(optarg) stream.str(optarg);
        switch(opt_result){
            /* INPUTS */
            case 'F':
                flip = atoi(optarg);
                std::cerr<<"Flip the waveform? "<<(flip?"yes":"no")<<std::endl;
                break;
            case 'Z':
                stream>>zoomTimeLeft>>zoomTimeRight;
                std::cerr<<"Zoom within time "<<zoomTimeLeft<<" ~ "<<zoomTimeRight<<" ns"<<std::endl;
                break;
            case 'r':
                runname = optarg;
                std::cerr<<"run name set to \""<<runname<<"\""<<std::endl;
                break;
            case 's':
                iEntryStart = atoi(optarg);
                std::cerr<<"Starting from entry "<<iEntryStart<<std::endl;
                break;
            case 'o':
                fOutputFile = optarg;
                std::cerr<<"Will save the histogram to the output file \""<<fOutputFile<<"\""<<std::endl;
                break;
            case 'p':
                if (!ichain_peak) ichain_peak = new TChain("t");
                ichain_peak->Add(optarg);
                std::cerr<<"Loading peaks from file "<<optarg<<std::endl;
            case 'n':
                nEntriesToScan = atoi(optarg);
                std::cerr<<"Will scan "<<nEntriesToScan<<" entries"<<std::endl;
                break;
            case 'b':
                batchMode = true;
                break;
            case '?':
                fprintf(stderr,"Wrong option! optopt=%c, optarg=%s\n", optopt, optarg);
            case 'h':
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    TChain* chain = new TChain("tr");
    for (int i = optind; i<argc; i++){
        chain->Add(argv[i]);
    }
    if(0 == chain->GetNtrees()) return 1;

    if (!batchMode) new TApplication("app", &argc, argv);

    // Set Branch Addresses
    Int_t run, bunch, event;
    std::vector< std::vector<Float_t>* > pwfs(4, nullptr);

    chain->SetBranchAddress("run", &run);
    //chain->SetBranchAddress("bunch", &bunch);
    chain->SetBranchAddress("event", &event);
    for(int ch=0; ch<4; ch++){
        chain->SetBranchAddress(Form("pwf%c", 'A'+ch), &pwfs[ch]);
    }

    std::vector<double> * start[4] = {0};
    std::vector<double> * height[4] = {0};
    int run_peak = 0;
    int event_peak = 0;
    double ped[4] = {0};
    if (ichain_peak){
        ichain_peak->SetBranchAddress("run",&run_peak);
        ichain_peak->SetBranchAddress("event",&event_peak);
        for (int ch = 0; ch<4; ch++){
            ichain_peak->SetBranchAddress(Form("start%c",'A'+ch),&start[ch]);
            ichain_peak->SetBranchAddress(Form("height%c",'A'+ch),&height[ch]);
            ichain_peak->SetBranchAddress(Form("ped%c",'A'+ch),&ped[ch]);
        }
    }

    Float_t               ps_dt = 0.8; // ns
    TH1F*                 hwf[4] = {};
    TH1F*                 hwfx10[4] = {};

    TCanvas* canvas = new TCanvas("canv","",1024,768);
    TPad* pads[4];
    canvas->cd();
    for(int i=0; i<4; i++){
        pads[i] = new TPad(Form("p%d", 4+i), "",
                i%2 * 0.5,    (1-i/2)*0.45,
                (i%2+1)*0.5,  (2-i/2)*0.45);
        pads[i]->SetGridx(1);
        pads[i]->SetGridy(1);
        pads[i]->Draw();
    }
    TLatex * text = new TLatex(0.1,0.92,"");
    text->SetTextFont(42);
    text->SetTextSize(0.023);
    text->Draw();

    TFile * ofile = 0;
    if (fOutputFile != ""){
        ofile = new TFile(fOutputFile,"recreate");
    }

    TTree* lastTree = nullptr;
    Long64_t iEntryStop = chain->GetEntries()-1;
    if (nEntriesToScan&&iEntryStart+nEntriesToScan-1<iEntryStop) iEntryStop = iEntryStart+nEntriesToScan-1;
    Long64_t entry_peak = 0;
    for(Long64_t entry=iEntryStart; entry<=iEntryStop; entry++){
        chain->GetEntry(entry);
        if (ichain_peak){
            while((run_peak!=run||event_peak!=event)&&entry_peak<ichain_peak->GetEntries()){
                ichain_peak->GetEntry(entry_peak);
                entry_peak++;
            }
            if (run_peak!=run||event_peak!=event){
                std::cerr<<"Cannot find run "<<run<<" event "<<event<<" from the input peak file!"<<std::endl;
                return -1;
            }
        }
        /*
           std::cout << "Run:" << run
           << " Bunch:" << bunch
           << " Event:" << event << std::endl;
         */
        //canvas->SetTitle(Form("Run%04d Bunch%06d Event%04d", run, bunch, event));
        //canvas->SetTitle(Form("Run%04d Bunch??? Event%04d", run, event));
        text->SetText(0.1,0.92,Form("Run%04d Bunch??? Event%04d Entry %08lld", run, event, entry));
        canvas->cd();text->Draw();

        //--------------------------------------------
        // Load user info
        //--------------------------------------------
        if(chain->GetTree() != lastTree){
            TTree* tree  = lastTree = chain->GetTree();
            TList* uinfo = tree->GetUserInfo();

            // PicoScope time step
            for(int ch=0; ch<4; ch++){
                if(pwfs[ch]){
                    if(not hwf[ch]){
                        hwf[ch] = new TH1F(Form("hwf_pico_ch%c", 'A'+ch), Form("Pico ch%c", 'A'+ch),
                                pwfs[ch]->size(), -t0[ch], -t0[ch]+ps_dt * pwfs[ch]->size());
                        hwfx10[ch] = new TH1F(Form("hwfx10_pico_ch%c", 'A'+ch), Form("Pico ch%c", 'A'+ch),
                                pwfs[ch]->size(), -t0[ch], -t0[ch]+ps_dt * pwfs[ch]->size());
                        hwfx10[ch]->SetLineColor(kGray);
                        hwf[ch]->SetStats(0);
                        hwfx10[ch]->SetStats(0);
                        if (flip){
                            if (ch>=2){
                                hwf[ch]->GetYaxis()->SetRangeUser(-100, 500);
                                hwfx10[ch]->GetYaxis()->SetRangeUser(-100, 500);
                            }
                            else{
                                hwf[ch]->GetYaxis()->SetRangeUser(-100, 1000);
                                hwfx10[ch]->GetYaxis()->SetRangeUser(-100, 1000);
                            }
                        }
                        else{
                            if (ch>=2){
                                hwf[ch]->GetYaxis()->SetRangeUser(-500, 100);
                                hwfx10[ch]->GetYaxis()->SetRangeUser(-500, 100);
                            }
                            else{
                                hwf[ch]->GetYaxis()->SetRangeUser(-1000, 100);
                                hwfx10[ch]->GetYaxis()->SetRangeUser(-1000, 100);
                            }
                        }
                    }
                    //hwf[ch]->GetXaxis()->Set(pwfs[ch]->size(), 0, ps_dt * pwfs[ch]->size());
                    //hwfx10[ch]->GetXaxis()->Set(pwfs[ch]->size(), 0, ps_dt * pwfs[ch]->size());
                    if (zoomTimeLeft<zoomTimeRight){
                        hwf[ch]->GetXaxis()->SetRangeUser(zoomTimeLeft, zoomTimeRight);
                        hwfx10[ch]->GetXaxis()->SetRangeUser(zoomTimeLeft, zoomTimeRight);
                    }
                    else{
                        hwf[ch]->GetXaxis()->UnZoom();
                        hwfx10[ch]->GetXaxis()->UnZoom();
                    }
                    pads[ch]->cd();
                    hwfx10[ch]->Draw("hist");
                    hwf[ch]->Draw("histsame");
                }
            }
        }


        for(int ch=0; ch<4; ch++){
            hwf[ch]->Reset();
            hwf[ch]->SetTitle(Form("Pico ch%c, event %lld, %d samples;Time [ns];voltage [mV]", 'A'+ch,entry,(pwfs[ch]?pwfs[ch]->size():0)));
            hwfx10[ch]->Reset();
            hwfx10[ch]->SetTitle(Form("Pico ch%c, event %lld, %d samples;Time [ns];voltage [mV]", 'A'+ch,entry,(pwfs[ch]?pwfs[ch]->size():0)));

            if(pwfs[ch]){
                for(size_t i=0; i<pwfs[ch]->size(); i++){
                    if (flip){
                        hwf[ch]->SetBinContent(i+1, -pwfs[ch]->at(i));
                        hwfx10[ch]->SetBinContent(i+1, -10*pwfs[ch]->at(i));
                    }
                    else{
                        hwf[ch]->SetBinContent(i+1, pwfs[ch]->at(i));
                        hwfx10[ch]->SetBinContent(i+1, 10*pwfs[ch]->at(i));
                    }
                }
                pads[ch]->Modified();
                pads[ch]->Update();
            }
        }
        if (ofile){
            for(int ch=0; ch<4; ch++){
                hwf[ch]->Write();
            }
        }
        if (ichain_peak){
            for (int ch = 0; ch<4; ch++){
                pads[ch]->cd();
                for (int i = 0; i<start[ch]->size(); i++){
                    double t = start[ch]->at(i)-t0[ch];
                    if (zoomTimeLeft<zoomTimeRight&&(t<zoomTimeLeft||t>zoomTimeRight)) continue;
                    double y0 = ped[ch];
                    double y1 = ped[ch]-height[ch]->at(i);
                    if (flip) {
                        y0*=-1; y1*=-1;
                    }
                    TLine * line = new TLine(t,y0,t,y1);
                    TLatex * text = new TLatex(t,y1,Form("%.0f %.0f",t,height[ch]->at(i)));
                    if (height[ch]->at(i)>hmin[ch]){
                        line->SetLineColor(kRed);
                        text->SetTextColor(kRed);
                    }
                    else{
                        line->SetLineColor(kBlue);
                        text->SetTextColor(kBlue);
                    }
                    line->Draw("SAME");
                    text->SetTextFont(42);
                    text->SetTextSize(0.027);
                    text->Draw("SAME");
                    double ythr = ped[ch]-hmin[ch];
                    if (flip) ythr*=-1;
                    double xleft = -t0[ch];
                    double xright = -t0[ch]+ps_dt*pwfs[ch]->size();
                    if (zoomTimeLeft<zoomTimeRight){
                        xleft = zoomTimeLeft;
                        xright = zoomTimeRight;
                    }
                    line = new TLine(xleft,ythr*10,xright,ythr*10);
                    line->SetLineColor(kGray);
                    line->SetLineStyle(2);
                    line->Draw("SAME");
                }
            }
        }
        canvas->Update();
        canvas->SaveAs(Form("results/%s_%d_%d.png",runname.Data(),run,event));
        if (!batchMode) canvas->WaitPrimitive();
    }

    if (ofile){
        ofile->Close();
    }
}
