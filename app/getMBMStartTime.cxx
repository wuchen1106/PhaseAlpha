#include <iostream>
#include <sstream>
#include <cstdlib>
#include <stdlib.h> /* atoi, atof */
#include <map>

#include <TTree.h>
#include <TF1.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>

int verbose = 0;
bool fDrawAllSpills = false;
bool fWithSmooth = false;

double AccCycle = 9.6e9;

TTree * otree = 0; 
int    spillID = 0;
double spillOffset = 0;
double spillWidth = 0;
int    spillNHitsBO = 0; // BO means within the beam on period
int    spillNHits = 0;
int    spillNTrigBO = 0;
int    spillNTrig = 0;
double spillTotalGaps = 0;
double spillMaxGap = 0;
double spillMinGap = 0;
bool   spillValid = 0;
int    maxGoodSpill = 0;
std::vector<int>          * hitEvent = 0;
std::vector<short int>    * hitType = 0;
std::vector<short int>    * hitChannel = 0;
std::vector<double> * hitTime = 0;
std::vector<double> * hitRawTime = 0;
std::vector<double> * eventTrigTime = 0;
std::vector<double> * eventLength = 0;
std::vector<int>    * eventNHits = 0;

int run = 0;

struct Hit{
    double t;
    short int ch;
    short int type;
    bool operator < (const Hit &b)const {
        return (t<b.t);
    };
    bool operator > (const Hit &b)const {
        return (t>b.t);
    };
    bool operator == (const Hit &b)const {
        return (t==b.t);
    };
    bool operator < (const Hit &b){
        return (t<b.t);
    };
    bool operator > (const Hit &b){
        return (t>b.t);
    };
    bool operator == (const Hit &b){
        return (t==b.t);
    };
    Hit(double at, short int ach, short int atype){
        t = at;
        ch = ach;
        type = atype;
    };
};

TH1D * anaHist(TH1D * hist, double & left, double & right, double ratioL = 0.01, double ratioR = 0.01){
    // now analyze the histogram to get the starting time and ending time
    TH1D * theHist = hist;
    if (fWithSmooth){
        theHist = new TH1D(*hist);
        theHist->Smooth();
    }
    double maximum = theHist->GetMaximum();
    bool inSpill = false;
    for (int ibin = 1; ibin<=theHist->GetNbinsX(); ibin++){
        double content = theHist->GetBinContent(ibin);
        if (!inSpill){
            if (content>maximum*ratioL*3){
                for (int jbin = ibin; jbin>=1; jbin--){
                    double content2 = theHist->GetBinContent(jbin);
                    if (content2<maximum*ratioL){
                        left = theHist->GetBinLowEdge(jbin);
                        inSpill = true;
                        break;
                    }
                    else if (jbin==1){
                        std::cerr<<"WARNING! Cannot find the left side of the spill "<<spillID<<" in run "<<run<<std::endl;
                        left = theHist->GetBinLowEdge(jbin);
                        inSpill = true;
                    }
                }
            }
        }
        else{
            if (content<maximum*ratioR){
                bool isEnd = true;
                for (int jbin = ibin; jbin<=theHist->GetNbinsX(); jbin++){
                    double content2 = theHist->GetBinContent(jbin);
                    if (content2>maximum*ratioR*3){
                        isEnd = false;
                        break;
                    }
                }
                if (isEnd){
                    right = theHist->GetBinLowEdge(ibin);
                    break;
                }
            }
        }
    }
    if (fWithSmooth) return theHist;
    return 0;
}

void drawHist(TH1D * hist, TH1D * histSmoothed, TH1D * histOff, TH1D * histT, int run, int spillID, double left, double right){
    TCanvas * canv = new TCanvas("canv","");
    canv->SetCanvasSize(1600,1200);
    canv->Divide(1,3);

    double maximum;

    canv->cd(1);
    int il = histT->GetXaxis()->FindBin(left);
    int ir = histT->GetXaxis()->FindBin(right);
    double NTrigTotal = histT->GetEntries();
    double NTrigInSpill = histT->Integral(il,ir);
    histT->SetTitle(Form("Spill #%d of run %d, %.0f triggers (%.2f%% off spill), %.2f Hz within spill",spillID,run,NTrigTotal,(NTrigTotal-NTrigInSpill)/NTrigTotal*100,NTrigInSpill/0.6));
    histT->SetLineColor(kGreen);
    histT->SetFillStyle(1001);
    histT->SetFillColor(kGreen);
    if (spillID!=-1) histT->GetXaxis()->SetRangeUser(-1e8,8e8);
    histT->Draw();
    maximum = histT->GetMaximum();
    TLine * line3 = new TLine(left,0,left,maximum);
    line3->SetLineColor(kRed);
    line3->Draw();
    TLine * line4 = new TLine(right,0,right,maximum);
    line4->SetLineColor(kRed);
    line4->Draw();
    canv->cd(2); //gPad->SetLogy(1);
    if (histSmoothed) maximum = histSmoothed->GetMaximum();
    else maximum = hist->GetMaximum();
    double NHitsTotal = hist->GetEntries();
    il = hist->GetXaxis()->FindBin(left);
    ir = hist->GetXaxis()->FindBin(right);
    double NHitsInSpill = hist->Integral(il,ir);
    hist->SetTitle(Form("Spill #%d of run %d, %.0f MBM hits (%.2f%% off spill), left %.2f s, right %.2f s, width %.2f s, %.2f kHz within spill",spillID,run,hist->GetEntries(),(NHitsTotal-NHitsInSpill)/NHitsTotal*100,left*1e-9,right*1e-9,(right-left)*1e-9,NHitsInSpill/0.6e3));
    hist->SetLineColor(kGray);
    hist->SetFillStyle(1001);
    hist->SetFillColor(kGray);
    if (spillID!=-1) hist->GetXaxis()->SetRangeUser(-1e8,8e8);
    hist->Draw();
    TLine * line1 = new TLine(left,0,left,maximum);
    line1->SetLineColor(kRed);
    line1->Draw();
    TLine * line2 = new TLine(right,0,right,maximum);
    line2->SetLineColor(kRed);
    line2->Draw();
    if (histSmoothed){
        histSmoothed->SetFillStyle(0);
        histSmoothed->SetLineColor(kBlue);
        histSmoothed->Draw("SAME");
    }
    canv->cd(3); gPad->SetLogy(1);
    histOff->SetTitle(Form("Sideband of Spill #%d of run %d, %.0f MBM hits, %.2f Hz",spillID,run,histOff->GetEntries(),histOff->GetEntries()/(AccCycle-1e9)*1e9));
    histOff->SetLineColor(kGray);
    histOff->SetFillStyle(1001);
    histOff->SetFillColor(kGray);
    histOff->Draw();
    maximum = histOff->GetMaximum();
    TLine * line5 = new TLine(left,0,left,maximum);
    line5->SetLineColor(kRed);
    line5->Draw();
    TLine * line6 = new TLine(right,0,right,maximum);
    line6->SetLineColor(kRed);
    line6->Draw();
    if (spillID==-2)
        canv->SaveAs(Form("pictures/GoodSpills.run%d.png",run));
    else if (spillID==-1)
        canv->SaveAs(Form("pictures/FirstSpill.run%d.png",run));
    else
        canv->SaveAs(Form("pictures/run%d.spill%d.png",run,spillID));
    delete line1;
    delete line2;
    delete line3;
    delete line4;
    delete line5;
    delete line6;
    delete canv;
    return;
}

void newSpill(int run, int& spillID, TH1D * hist, TH1D * histOff, TH1D * histT){
    double left = 0;
    double right = 0;
    TH1D * histSmoothed = anaHist(hist,left,right);
    if (fDrawAllSpills) drawHist(hist,histSmoothed,histOff,histT,run,spillID,left,right);
    if (histSmoothed) delete histSmoothed;
    spillOffset = left;
    spillWidth = right-left;
    int il = hist->GetXaxis()->FindBin(left);
    int ir = hist->GetXaxis()->FindBin(right);
    spillNHitsBO = hist->Integral(il,ir);
    spillNHits = hist->GetEntries();
    spillNTrigBO = histT->Integral(il,ir);
    spillNTrig = histT->GetEntries();
    // check if this spill is valid
    if (left<-1e9||spillWidth<0.3e9||spillWidth>0.72e9) spillValid = false;
    if (spillMaxGap>50e6) spillValid = false;
    if (spillValid) maxGoodSpill++;
    otree->Fill();
    // reset
    spillID++;
    spillMaxGap = 0;
    spillMinGap = 1e10;
    hitRawTime->clear();
    hitTime->clear();
    hitType->clear();
    hitEvent->clear();
    hitChannel->clear();
    eventTrigTime->clear();
    eventLength->clear();
    eventNHits->clear();
    hist->Reset();
    histOff->Reset();
    histT->Reset();
}

void print_usage(char * prog_name){
    fprintf(stderr,"Usage %s [options] input\n",prog_name);
    fprintf(stderr,"[options]\n");
};

int main(int argc, char** argv)
{
    // Load options
    int    opt_result;
    std::stringstream stream;
    int theRun = 0;
    while((opt_result=getopt(argc,argv,"r:d:s:v:h"))!=-1){
        stream.clear();
        if(optarg) stream.str(optarg);
        switch(opt_result){
            /* INPUTS */
            case 'r':
                theRun = atoi(optarg);
                std::cerr<<"run number "<<theRun<<std::endl;
                break;
            case 'd':
                fDrawAllSpills = atoi(optarg);
                std::cerr<<"draw hists for all spills? "<<(fDrawAllSpills?"yes":"no")<<std::endl;
                break;
            case 's':
                fWithSmooth = atoi(optarg);
                std::cerr<<"smoothing hists? "<<(fWithSmooth?"yes":"no")<<std::endl;
                break;
            case 'v':
                verbose = atoi(optarg);
                std::cerr<<"verbose level set to "<<verbose<<std::endl;
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
    chain->Add(Form("decdata/rawdata%05d.root",theRun));
    if(0 == chain->GetNtrees()) return 1;

    int event = 0;
    std::vector<unsigned char> * mbm_x_ch = 0;
    std::vector<unsigned char> * mbm_y_ch = 0;
    std::vector<double> * mbm_x_time = 0;
    std::vector<double> * mbm_y_time = 0;

    chain->SetBranchAddress("run",&run);
    chain->SetBranchAddress("event",&event);
    bool hasMBM = true;
    if (!chain->GetListOfBranches()->FindObject("mbm.x.ch")) hasMBM = false;
    if (!chain->GetListOfBranches()->FindObject("mbm.x.time")) hasMBM = false;
    if (!chain->GetListOfBranches()->FindObject("mbm.y.ch")) hasMBM = false;
    if (!chain->GetListOfBranches()->FindObject("mbm.y.time")) hasMBM = false;
    if (!hasMBM){
        std::cerr<<"There are no mbm related branches (or incomplete) in run "<<theRun<<std::endl;
        return -1;
    }
    chain->SetBranchAddress("mbm.x.ch",&mbm_x_ch);
    chain->SetBranchAddress("mbm.y.ch",&mbm_y_ch);
    chain->SetBranchAddress("mbm.x.time",&mbm_x_time);
    chain->SetBranchAddress("mbm.y.time",&mbm_y_time);

    //===============================================================================================================
    // First of all, analyze the first spill to get left0 if it's not provided
    //   make a histogram to check the position of the first event in the run relative to the current spill
    // also, check the maximum spill ID
    bool firstEntry = true;// don't add offset for the first entry
    double t0 = 0;
    std::vector<Hit> * mbm_all_hits = new std::vector<Hit>;
    TH1D * histHits = new TH1D("histHits",Form("run %d;Connected Time [ns];Count per 10 msec",theRun),200,-1e9,1e9);
    TH1D * histHitsOff = new TH1D("histHitsOff",Form("run %d;Connected Time [ns];Count per 100 msec",theRun),100,-5e9,5e9);
    TH1D * histTrig = new TH1D("histTrig",Form("run %d;Connected Time [ns];Count per 10 msec",theRun),2000,-1e9,1e9);
    spillID = 0;
    for (Long64_t entry = 0; entry<chain->GetEntries(); entry++){
        chain->GetEntry(entry);
        if (run!=theRun){
            std::cerr<<"ERROR! theRun "<<theRun<<", run @ entry "<<entry<<" is "<<run<<std::endl;
            return -1;
        }

        // load all the valid time stamps from x channels and y channels
        mbm_all_hits->clear();
        for (unsigned int i = 0; i<mbm_x_time->size(); i++){
            if (mbm_x_time->at(i)>0) continue;
            mbm_all_hits->push_back(Hit(mbm_x_time->at(i),mbm_x_ch->at(i),0));
        }
        for (unsigned int i = 0; i<mbm_y_time->size(); i++){
            if (mbm_y_time->at(i)>0) continue;
            mbm_all_hits->push_back(Hit(mbm_y_time->at(i),mbm_y_ch->at(i),1));
        }
        if (mbm_all_hits->size()<=0) continue;
        // sort the stamps
        std::sort(mbm_all_hits->begin(),mbm_all_hits->end());
        // use the negative value of the first time stamp as the offset for the following stamps
        if (!firstEntry){
            t0 += -mbm_all_hits->at(0).t;
            if (mbm_all_hits->size()>1){
                t0 += mbm_all_hits->at(1).t-mbm_all_hits->at(0).t;
            }
            else if (mbm_all_hits->size()>0){
                t0 += -mbm_all_hits->at(0).t/mbm_all_hits->size();
            }
        }
        firstEntry = false;
        // check if it's already outside of the spill
        for (unsigned int i = 0; i<mbm_all_hits->size(); i++){
            double localTime = mbm_all_hits->at(i).t + t0-spillID*AccCycle;
            if (localTime>AccCycle/2){
                spillID++;
            }
            // save the stamps of the first spill to the histogram
            if (spillID==0){
                histHits->Fill(localTime);
            }
        }
    }
    int maxSpill = spillID;
    // TODO: add in the calibrated values
    double left0, right0;
    TH1D * histNew = anaHist(histHits,left0,right0);
    if (fDrawAllSpills) drawHist(histHits,histNew,histHitsOff,histTrig,run,-1,left0,right0);
    if (histNew) delete histNew;

    //===============================================================================================================
    // Secondly, visit the whole run again, save information in spills
    // each spill needs to be re visited
    TFile * ofile = new TFile(Form("results/run%05d.spill.root",theRun),"recreate");
    otree = new TTree("t","t");
    otree->Branch("run",&run);
    otree->Branch("spill",&spillID);
    otree->Branch("offset",&spillOffset);
    otree->Branch("width",&spillWidth);
    otree->Branch("nHits",&spillNHitsBO);
    otree->Branch("nHitsAll",&spillNHits);
    otree->Branch("nTrig",&spillNTrigBO);
    otree->Branch("nTrigAll",&spillNTrig);
    otree->Branch("totalGaps",&spillTotalGaps);
    otree->Branch("maxGap",&spillMaxGap);
    otree->Branch("minGap",&spillMinGap);
    otree->Branch("valid",&spillValid);
    otree->Branch("hitEvent",&hitEvent);
    otree->Branch("hitType",&hitType);
    otree->Branch("hitChannel",&hitChannel);
    otree->Branch("hitTime",&hitTime);
    otree->Branch("hitRawTime",&hitRawTime); // FIXME: commented for now to save disk space, but this might be useful for debugging
    otree->Branch("eventTrigTime",&eventTrigTime);
    otree->Branch("eventLength",&eventLength);
    otree->Branch("eventNHits",&eventNHits);
    hitTime = new std::vector<double>;
    hitEvent = new std::vector<int>;
    hitType = new std::vector<short int>;
    hitChannel = new std::vector<short int>;
    hitRawTime = new std::vector<double>;
    eventNHits = new std::vector<int>;
    eventTrigTime = new std::vector<double>;
    eventLength = new std::vector<double>;
    firstEntry = true;// don't add offset for the first entry
    spillID = 0;
    spillTotalGaps = 0;
    spillMaxGap = 0;
    spillMinGap = 1e10;
    spillValid = true;
    t0 = -left0;
    delete histHits;
    histHits = new TH1D("histHits",Form("run %d;Connected Time [ns];Count per 1 msec;Spill ID",theRun),2000,-1e9,1e9);
    TH2D * histHits2D = new TH2D("histHits2D",Form("run %d;Connected Time [ns];Count per 1 msec;Spill ID",theRun),2000,-1e9,1e9,maxSpill,0,maxSpill);
    TH2D * histHitsOff2D = new TH2D("histHitsOff2D",Form("run %d;Connected Time [ns];Count per 100 msec;Spill ID",theRun),100,-5e9,5e9,maxSpill,0,maxSpill);
    TH2D * histTrig2D = new TH2D("histTrig2D",Form("run %d;Connected Time [ns];Count per 1 msec;Spill ID",theRun),2000,-1e9,1e9,maxSpill,0,maxSpill);
    TH2D * histTrigOff2D = new TH2D("histTrigOff2D",Form("run %d;Connected Time [ns];Count per 100 msec;Spill ID",theRun),100,-5e9,5e9,maxSpill,0,maxSpill);
    TH2D * hSpillVSGap = new TH2D("hSpillVSGap",";log_{10}(Event Gap [ns]);Spill ID",100,1,10,maxSpill,0,maxSpill);
    TH1D * hGaps = new TH1D("hGaps",";log_{10}(Event Gap [ns])",100,1,10);
    TH1D * hLength = new TH1D("hLength",";log_{10}(Event Length [ns])",100,1,10);
    maxGoodSpill = -1;
    for (Long64_t entry = 0; entry<chain->GetEntries(); entry++){
        chain->GetEntry(entry);
        if (run!=theRun){
            std::cerr<<"ERROR! theRun "<<theRun<<", run @ entry "<<entry<<" is "<<run<<std::endl;
            return -1;
        }
        
        // load all the valid time stamps from x channels and y channels
        mbm_all_hits->clear();
        for (unsigned int i = 0; i<mbm_x_time->size(); i++){
            if (mbm_x_time->at(i)>0) continue;
            mbm_all_hits->push_back(Hit(mbm_x_time->at(i),mbm_x_ch->at(i),0));
        }
        for (unsigned int i = 0; i<mbm_y_time->size(); i++){
            if (mbm_y_time->at(i)>0) continue;
            mbm_all_hits->push_back(Hit(mbm_y_time->at(i),mbm_y_ch->at(i),1));
        }
        if (mbm_all_hits->size()<=0) continue;
        // sort the stamps
        std::sort(mbm_all_hits->begin(),mbm_all_hits->end());
        // use the negative value of the first time stamp as the offset for the following stamps
        if (!firstEntry){
            t0 += -mbm_all_hits->at(0).t;
            double additionalGap = 0;
            if (mbm_all_hits->size()>1){
                additionalGap = mbm_all_hits->at(1).t-mbm_all_hits->at(0).t;
            }
            else if (mbm_all_hits->size()>0){
                additionalGap = -mbm_all_hits->at(0).t/mbm_all_hits->size();
            }
            hSpillVSGap->Fill(log(additionalGap)/log(10),spillID);
            hGaps->Fill(log(additionalGap)/log(10));
            t0+=additionalGap;
            if (additionalGap>spillMaxGap) spillMaxGap = additionalGap;
            if (additionalGap<spillMinGap) spillMinGap = additionalGap;
            if (verbose>1){
                std::cerr<<"entry "<<entry<<" run "<<run<<" event "<<event<<" spill "<<spillID<<" additionalGap = "<<additionalGap<<std::endl;
            }
            spillTotalGaps += additionalGap;
        }
        // save the stamps to the histogram
        for (unsigned int i = 0; i<mbm_all_hits->size(); i++){
            double localTime = mbm_all_hits->at(i).t + t0-spillID*AccCycle;
            // new spill?
            if (localTime>AccCycle/2){
                if (verbose>1) std::cerr<<"Entry "<<entry<<" t "<<mbm_all_hits->at(i).t<<" local t "<<localTime<<" end of spill "<<spillID<<std::endl;
                newSpill(run,spillID,histHits,histHitsOff,histTrig);
            }
            hitTime->push_back(localTime);
            hitRawTime->push_back(mbm_all_hits->at(i).t);
            hitType->push_back(mbm_all_hits->at(i).type);
            hitEvent->push_back(event);
            hitChannel->push_back(mbm_all_hits->at(i).ch);
            if (i==mbm_all_hits->size()-1){
                double trigT = t0-spillID*AccCycle-600; // FIXME: the t0 for trigger shall be calibrated
                eventTrigTime->push_back(trigT);
                eventNHits->push_back(mbm_all_hits->size());
                eventLength->push_back(-mbm_all_hits->at(0).t);
                hLength->Fill(log(-mbm_all_hits->at(0).t)/log(10));
                histTrig->Fill(trigT);
                histTrig2D->Fill(trigT,spillID);
                histTrigOff2D->Fill(trigT,spillID);
                if (entry == chain->GetEntries()-1){
                    if (verbose>1) std::cerr<<"Last entry "<<entry<<" t "<<mbm_all_hits->at(i).t<<" local t "<<localTime<<" end of spill "<<spillID<<std::endl;
                    newSpill(run,spillID,histHits,histHitsOff,histTrig);
                }
            }
            if (firstEntry&&localTime<-AccCycle/2) continue; // Belong to the previous spill. don't use them
            histHits->Fill(localTime);
            histHits2D->Fill(localTime,spillID);
            histHitsOff->Fill(localTime);
            histHitsOff2D->Fill(localTime,spillID);
        }
        firstEntry = false;
    }

    TH1D * histHitsStacked = histHits2D->ProjectionX("histHitsStacked",1,maxGoodSpill+1);
    TH1D * histHitsOffStacked = histHitsOff2D->ProjectionX("histHitsOffStacked",1,maxGoodSpill+1);
    TH1D * histTrigStacked = histTrig2D->ProjectionX("histTrigStacked",1,maxGoodSpill+1);
    double left = 0;
    double right = 0;
    TH1D * histSmoothed = anaHist(histHitsStacked,left,right);
    drawHist(histHitsStacked,histSmoothed,histHitsOffStacked,histTrigStacked,run,-2,left,right); // -1 means all
    if (histSmoothed) delete histSmoothed;
    std::cout<<theRun<<" "<<maxSpill<<" "<<maxGoodSpill<<" "<<left0<<" "<<left<<" "<<right<<std::endl;

    TCanvas * canv = new TCanvas("canv2","");
    canv->SetCanvasSize(1600,1200);
    canv->Divide(2,2);
    histHits2D->GetXaxis()->SetRangeUser(-1e8,8e8);
    histTrig2D->GetXaxis()->SetRangeUser(-1e8,8e8);
    canv->cd(1); histHitsOff2D->Draw("COLZ");
    canv->cd(2); histHits2D->Draw("COLZ");
    canv->cd(3); histTrigOff2D->Draw("COLZ");
    canv->cd(4); histTrig2D->Draw("COLZ");
    canv->SaveAs(Form("pictures/2Dspill.run%d.png",run));

    canv = new TCanvas("canv3","");
    canv->SetCanvasSize(1600,1600);
    canv->Divide(2,2);
    canv->cd(1);
    hSpillVSGap->Draw("COLZ");
    otree->SetMarkerStyle(24); otree->Draw("spill+0.5:log(maxGap)/log(10)","","PSAME");
    otree->SetMarkerStyle(20); otree->Draw("spill+0.5:log(maxGap)/log(10)","valid","PSAME");
    canv->cd(2);
    TH2D * hSpillVSTotalGap = new TH2D("hSpillVSTotalGap",";Accumulated Gap [ns];Spill ID",1000,1,1e8,maxSpill,0,maxSpill);
    hSpillVSTotalGap->Draw();
    otree->SetMarkerStyle(24); otree->Draw("spill+0.5:totalGaps","","PLSAME");
    otree->SetMarkerStyle(20); otree->Draw("spill+0.5:totalGaps","valid","PSAME");
    canv->cd(3);
    hGaps->Draw();
    canv->cd(4);
    hLength->Draw();
    canv->SaveAs(Form("pictures/gaps.run%d.png",run));

    canv = new TCanvas("canv4","");
    canv->SetCanvasSize(1200,600);
    canv->Divide(3,1);
    TH1D * histNHits[3];
    int    nHitsCounter[3] = {0};
    int    stepSize[3] = {0};
    bool   firstBin[3] = {0};
    TF1 * f[3];
    for (int i = 0; i<3; i++){
        stepSize[i] = pow(10,i);
        int nHitsMax = 200*stepSize[i];
        firstBin[i] = true;
        f[i] = new TF1("f0","TMath::PoissonI(x,[0])*[1]",0,nHitsMax);
        f[i]->SetNpx(1000);
        histNHits[i] = new TH1D(Form("histNHits%d",i),Form("Number of hits per %d msec;Number of hits",stepSize[i]),200,0,nHitsMax);
    }
    int leftBin = histHits2D->GetXaxis()->FindBin(left+1e8);
    int rightBin = histHits2D->GetXaxis()->FindBin(right-1e8);
    for (int iBiny = 1; iBiny<=maxGoodSpill+1; iBiny++){
        for (int iBinx = leftBin; iBinx<=rightBin; iBinx++){
            for (int i = 0; i<3; i++){
                if ((iBinx-leftBin)%stepSize[i]==0&&!firstBin[i]){
                    histNHits[i]->Fill(nHitsCounter[i]);
                    nHitsCounter[i] = 0;
                }
                nHitsCounter[i] += histHits2D->GetBinContent(iBinx,iBiny);
                firstBin[i] = false;
            }
        }
    }
    for (int i = 0; i<3; i++){
        f[i]->SetParameters(histNHits[i]->GetMean(),1);
        f[i]->SetParameters(histNHits[i]->GetMean(),histNHits[i]->GetMaximum()/f[i]->Eval(histNHits[i]->GetMean()));
    }
    canv->cd(1);
    histNHits[0]->Draw();
    f[0]->Draw("SAME");
    canv->cd(2);
    histNHits[1]->Draw();
    f[1]->Draw("SAME");
    canv->cd(3);
    histNHits[2]->Draw();
    f[2]->Draw("SAME");
    canv->SaveAs(Form("pictures/nHits.run%d.png",run));

    otree->Write();
    ofile->Close();

    return 0;
}
