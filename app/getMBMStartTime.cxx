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

int verbose = 0;

double AccCycle = 9.61215e9;

TTree * otree = 0; 
TH1D * histHits = 0;
TH1D * histHitsOff = 0;
TH1D * histTrig = 0;
int    spillID = 0;
double spillOffset = 0;
double spillWidth = 0;
int    spillNHitsBO = 0; // BO means within the beam on period
int    spillNHits = 0;
int    spillNTrigBO = 0;
int    spillNTrig = 0;
std::vector<int>    * hitType = 0;
std::vector<int>    * hitChannel = 0;
std::vector<double> * hitTime = 0;
std::vector<double> * trigTime = 0;

int run = 0;

TH1D * anaHist(TH1D * hist, double & left, double & right){
    // now analyze the histogram to get the starting time and ending time
    TH1D * histNew = new TH1D(*hist);
    histNew->Smooth();
    double maximum = histNew->GetMaximum();
    bool inSpill = false;
    for (int ibin = 1; ibin<=histNew->GetNbinsX(); ibin++){
        double content = histNew->GetBinContent(ibin);
        if (!inSpill){
            if (content>maximum/3){
                for (int jbin = ibin; jbin>=1; jbin--){
                    double content2 = histNew->GetBinContent(jbin);
                    if (content2<maximum/100){
                        left = histNew->GetBinLowEdge(jbin);
                        inSpill = true;
                        break;
                    }
                }
            }
        }
        else{
            if (content<maximum/100){
                bool isEnd = true;
                for (int jbin = ibin; jbin<=histNew->GetNbinsX(); jbin++){
                    double content2 = histNew->GetBinContent(jbin);
                    if (content2>maximum/3){
                        isEnd = false;
                        break;
                    }
                }
                if (isEnd){
                    right = histNew->GetBinLowEdge(ibin);
                    break;
                }
            }
        }
    }
    return histNew;
}

void drawHist(TH1D * hist, TH1D * histNew, TH1D * histOff, TH1D * histTrig, int run, int spillID, double left, double right){
    TCanvas * canv = new TCanvas("canv","",1024,1024);
    canv->Divide(1,3);
    canv->cd(1);
    double maximum = histNew->GetMaximum();
    double NHitsTotal = hist->GetEntries();
    int il = hist->GetXaxis()->FindBin(left);
    int ir = hist->GetXaxis()->FindBin(right);
    double NHitsInSpill = hist->Integral(il,ir);
    hist->SetTitle(Form("Spill #%d of run %d, %.0f MBM hits (%.2f%% off spill), left %.2f s, right %.2f s, width %.2f s, %.2f kHz within spill",spillID,run,hist->GetEntries(),(NHitsTotal-NHitsInSpill)/NHitsTotal*100,left*1e-9,right*1e-9,(right-left)*1e-9,NHitsInSpill/0.6e3));
    hist->SetLineColor(kGray);
    hist->SetFillStyle(1001);
    hist->SetFillColor(kGray);
    hist->Draw();
    TLine * line1 = new TLine(left,0,left,maximum);
    line1->SetLineColor(kRed);
    line1->Draw();
    TLine * line2 = new TLine(right,0,right,maximum);
    line2->SetLineColor(kRed);
    line2->Draw();
    histNew->SetFillStyle(0);
    histNew->SetLineColor(kBlue);
    histNew->Draw("SAME");
    canv->cd(2);
    il = histTrig->GetXaxis()->FindBin(left);
    ir = histTrig->GetXaxis()->FindBin(right);
    double NTrigTotal = histTrig->GetEntries();
    double NTrigInSpill = histTrig->Integral(il,ir);
    histTrig->SetTitle(Form("Spill #%d of run %d, %.0f triggers (%.2f%% off spill), %.2f Hz within spill",spillID,run,NTrigTotal,(NTrigTotal-NTrigInSpill)/NTrigTotal*100,NTrigInSpill/0.6));
    histTrig->SetLineColor(kGreen);
    histTrig->SetFillStyle(1001);
    histTrig->SetFillColor(kGreen);
    histTrig->Draw();
    maximum = histTrig->GetMaximum();
    TLine * line3 = new TLine(left,0,left,maximum);
    line3->SetLineColor(kRed);
    line3->Draw();
    TLine * line4 = new TLine(right,0,right,maximum);
    line4->SetLineColor(kRed);
    line4->Draw();
    canv->cd(3);
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
    canv->SaveAs(Form("pictures/MBM/run%d.spill%d.png",run,spillID));
    delete line1;
    delete line2;
    delete line3;
    delete line4;
    delete line5;
    delete line6;
    delete canv;
    return;
}

void newSpill(){
    double left = 0;
    double right = 0;
    TH1D * histNew = anaHist(histHits,left,right);
    drawHist(histHits,histNew,histHitsOff,histTrig,run,spillID,left,right);
    delete histNew;
    spillID++;
    spillOffset = left;
    spillWidth = right-left;
    int il = histHits->GetXaxis()->FindBin(left);
    int ir = histHits->GetXaxis()->FindBin(right);
    spillNHitsBO = histHits->Integral(il,ir);
    spillNHits = histHits->GetEntries();
    spillNTrigBO = histTrig->Integral(il,ir);
    spillNTrig = histTrig->GetEntries();
    otree->Fill();
    // reset
    hitTime->clear();
    hitType->clear();
    hitChannel->clear();
    trigTime->clear();
    histHits->Reset();
    histHitsOff->Reset();
    histTrig->Reset();
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
    while((opt_result=getopt(argc,argv,"r:v:h"))!=-1){
        stream.clear();
        if(optarg) stream.str(optarg);
        switch(opt_result){
            /* INPUTS */
            case 'r':
                theRun = atoi(optarg);
                std::cerr<<"run number "<<theRun<<std::endl;
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
    chain->SetBranchAddress("mbm.x.ch",&mbm_x_ch);
    chain->SetBranchAddress("mbm.y.ch",&mbm_y_ch);
    chain->SetBranchAddress("mbm.x.time",&mbm_x_time);
    chain->SetBranchAddress("mbm.y.time",&mbm_y_time);

    TFile * ofile = new TFile(Form("results/run%05d.spill.root",theRun),"recreate");
    otree = new TTree("t","t");
    otree->Branch("offset",&spillOffset);
    otree->Branch("width",&spillWidth);
    otree->Branch("nHits",&spillNHitsBO);
    otree->Branch("nHitsAll",&spillNHits);
    otree->Branch("nTrig",&spillNTrigBO);
    otree->Branch("nTrigAll",&spillNTrig);
    otree->Branch("hitType",&hitType);
    otree->Branch("hitChannel",&hitChannel);
    otree->Branch("hitTime",&hitTime);
    otree->Branch("trigTime",&trigTime);
    hitTime = new std::vector<double>;
    hitType = new std::vector<int>;
    hitChannel = new std::vector<int>;
    trigTime = new std::vector<double>;
    std::vector<int> * mbm_all_type = new std::vector<int>;
    std::vector<int> * mbm_all_channel = new std::vector<int>;

    bool firstEntry = true;// don't add offset for the first entry
    double t0 = 0;
    std::vector<double> * mbm_all_time = new std::vector<double>;
    histHits = new TH1D("histHits",Form("run %d;Connected Time [ns];Count per msec",theRun),1000,-2e8,8e8);
    histHitsOff = new TH1D("histHitsOff",Form("run %d;Connected Time [ns];Count per 100 msec",theRun),100,-5e9,5e9);
    histTrig = new TH1D("histTrig",Form("run %d;Connected Time [ns];Count per 10 msec",theRun),100,-2e8,8e8);
    for (Long64_t entry = 0; entry<chain->GetEntries(); entry++){
        chain->GetEntry(entry);
        if (run!=theRun){
            std::cerr<<"ERROR! theRun "<<theRun<<", run @ entry "<<entry<<" is "<<run<<std::endl;
            return -1;
        }
        
        // load all the valid time stamps from x channels and y channels
        mbm_all_time->clear();
        mbm_all_type->clear();
        mbm_all_channel->clear();
        for (unsigned int i = 0; i<mbm_x_time->size(); i++){
            if (mbm_x_time->at(i)>0) continue;
            mbm_all_time->push_back(mbm_x_time->at(i));
            mbm_all_type->push_back(0);
            mbm_all_channel->push_back(mbm_x_ch->at(i));
        }
        for (unsigned int i = 0; i<mbm_y_time->size(); i++){
            if (mbm_y_time->at(i)>0) continue;
            mbm_all_time->push_back(mbm_y_time->at(i));
            mbm_all_type->push_back(1);
            mbm_all_channel->push_back(mbm_y_ch->at(i));
        }
        if (mbm_all_time->size()<=0) continue;
        // sort the stamps
        std::sort(mbm_all_time->begin(),mbm_all_time->end());
        // use the negative value of the first time stamp as the offset for the following stamps
        if (!firstEntry){
            t0 += -mbm_all_time->at(0);
            if (mbm_all_time->size()>0){
                t0 += -mbm_all_time->at(0)/mbm_all_time->size(); // use average
                //t0 += mbm_all_time->at(1)-mbm_all_time->at(0); // repeat the first recorded gap
            }
        }
        // save the stamps to the histogram
        for (unsigned int i = 0; i<mbm_all_time->size(); i++){
            double localTime = mbm_all_time->at(i) + t0-spillID*AccCycle;
            if (firstEntry&&localTime<-AccCycle/2) continue; // Belong to the previous spill. don't use them
            // new spill?
            if (localTime>AccCycle/2){
                std::cout<<"Entry "<<entry<<" t "<<mbm_all_time->at(i)<<" local t "<<localTime<<" end of spill "<<spillID<<std::endl;
                newSpill();
            }
            hitTime->push_back(localTime);
            hitType->push_back(mbm_all_type->at(i));
            hitChannel->push_back(mbm_all_channel->at(i));
            if (i==mbm_all_time->size()-1){
                double trigT = t0-spillID*AccCycle-600; // FIXME: the t0 for trigger shall be calibrated
                trigTime->push_back(trigT);
                histTrig->Fill(trigT);
                if (entry == chain->GetEntries()-1){
                    std::cout<<"Last entry "<<entry<<" t "<<mbm_all_time->at(i)<<" local t "<<localTime<<" end of spill "<<spillID<<std::endl;
                    newSpill();
                }
            }
            histHits->Fill(localTime);
            if (localTime<-2e8||localTime>8e8) histHitsOff->Fill(localTime);
        }
        firstEntry = false;
    }
    otree->Write();
    ofile->Close();

    return 0;
}
