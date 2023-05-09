#include <iostream>
#include <sstream>
#include <cstdlib>
#include <stdlib.h> /* atoi, atof */
#include <iomanip> // for setprecision

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

void print_usage(char * prog_name){
    fprintf(stderr,"Usage %s [options] runname\n",prog_name);
    fprintf(stderr,"[options]\n");
};

int main(int argc, char** argv)
{

    int iEntryStart = 0;
    int nSamplesBeforeSignal = 450;
    Float_t               ps_dt = 0.8; // ns
    int                   ps_size = 12500;
    double                threshold[4] = {1};
    double                thresholdSignal[4] = {1};
    double                resistance = 50;// Ohm
    TString runname = "";

    // Load options
    int    opt_result;
    std::stringstream stream;
    while((opt_result=getopt(argc,argv,"T:S:s:r:h"))!=-1){
        stream.clear();
        if(optarg) stream.str(optarg);
        switch(opt_result){
            /* INPUTS */
            case 'T':
                stream>>threshold[0]>>threshold[1]>>threshold[2]>>threshold[3];
                std::cerr<<"Thresholds set to "<<threshold[0]<<" "<<threshold[1]<<" "<<threshold[2]<<" "<<threshold[3]<<" mV"<<std::endl;
                break;
            case 'S':
                stream>>thresholdSignal[0]>>thresholdSignal[1]>>thresholdSignal[2]>>thresholdSignal[3];
                std::cerr<<"Thresholds for signal set to "<<thresholdSignal[0]<<" "<<thresholdSignal[1]<<" "<<thresholdSignal[2]<<" "<<thresholdSignal[3]<<" pC"<<std::endl;
                break;
            case 's':
                iEntryStart = atoi(optarg);
                std::cerr<<"Starting from entry "<<iEntryStart<<std::endl;
                break;
            case 'r':
                runname = optarg;
                std::cerr<<"Run name \""<<runname<<"\""<<std::endl;
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
    if (chain->GetEntries()==0) return 1;

    // Set Branch Addresses
    Int_t run = 0;
    Int_t event = 0;
    std::vector< std::vector<Float_t>* > pwfs(4, nullptr);

    chain->SetBranchAddress("run", &run);
    chain->SetBranchAddress("event", &event);
    for(int ch=0; ch<4; ch++){
        chain->SetBranchAddress(Form("pwf%c", 'A'+ch), &pwfs[ch]);
    }

    double step[4] = {1e9,1e9,1e9,1e9};
    double minimum[4] = {1e9,1e9,1e9,1e9};
    for (Long64_t entry = iEntryStart; entry<100&&entry<chain->GetEntries(); entry++){
        chain->GetEntry();
        for (int ch = 0; ch<4; ch++){
            if (pwfs[ch]){
                double min = 1e9;
                double next2min = 1e9;
                for (int i = 0; i<100&&i<pwfs[ch]->size(); i++){
                    double value = pwfs[ch]->at(i);
                    if (value<min) {next2min = min; min = value;}
                    if (value>min&&value<next2min) next2min = value;
                }
                minimum[ch] = min;
                step[ch] = next2min-min;
            }
        }
    }
    for (int ch = 0; ch<4; ch++){
        if (step[ch]==0) continue;
        while(minimum[ch]<0) minimum[ch]+=step[ch];
        while(minimum[ch]-step[ch]>0) minimum[ch]-=step[ch];
    }

    std::cout<<runname;
    for (int ch = 0; ch<4; ch++){
        std::cout<<" "<<std::setprecision(20)<<step[ch];
    }
    for (int ch = 0; ch<4; ch++){
        std::cout<<" "<<std::setprecision(20)<<minimum[ch];
    }
    std::cout<<std::endl;

    return 0;
}
