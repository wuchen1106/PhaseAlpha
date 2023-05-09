#include <iostream>
#include <sstream>
#include <cstdlib>
#include <stdlib.h> /* atoi, atof */

#include "TChain.h"
#include "TString.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"

void print_usage(char * prog_name){
    fprintf(stderr,"Usage %s [options] runname\n",prog_name);
    fprintf(stderr,"[options]\n");
};

int main(int argc, char ** argv){
    TString runname;
    TString outputname;

    // Load options
    int    opt_result;
    std::stringstream stream;
    while((opt_result=getopt(argc,argv,"r:o:h"))!=-1){
        stream.clear();
        if(optarg) stream.str(optarg);
        switch(opt_result){
            /* INPUTS */
            case 'r':
                runname = optarg;
                std::cerr<<"Run name \""<<runname<<"\""<<std::endl;
                break;
            case 'o':
                outputname = optarg;
                std::cerr<<"output name \""<<runname<<"\""<<std::endl;
                break;
            case '?':
                fprintf(stderr,"Wrong option! optopt=%c, optarg=%s\n", optopt, optarg);
            case 'h':
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    TChain* t = new TChain("t");
    for (int i = optind; i<argc; i++){
        t->Add(argv[i]);
    }
    if(0 == t->GetNtrees()) return 1;

    int foundTrigger = 0;
    int trig_n0 =0;
    int trig_np =0;
    int trig_n =0;
    double trig_h[4] = {0};
    double trig_t = 0;
    std::vector<int> * evt_n0 =0;
    std::vector<int> * evt_np =0;
    std::vector<int> * evt_n =0;
    std::vector<double> * evt_h[4] = {0};
    std::vector<double> * evt_t = 0;
    t->SetBranchAddress("foundTrigger",&foundTrigger);
    t->SetBranchAddress("trig_n0",&trig_n0);
    t->SetBranchAddress("trig_np",&trig_np);
    t->SetBranchAddress("trig_n",&trig_n);
    t->SetBranchAddress("trig_t",&trig_t);
    t->SetBranchAddress("evt_n0",&evt_n0);
    t->SetBranchAddress("evt_np",&evt_np);
    t->SetBranchAddress("evt_n",&evt_n);
    t->SetBranchAddress("evt_t",&evt_t);
    for (int ch = 0; ch<4; ch++){
        t->SetBranchAddress(Form("trig_h%c",'A'+ch),&trig_h[ch]);
        t->SetBranchAddress(Form("evt_h%c",'A'+ch),&evt_h[ch]);
    }

    TH2D * hist = new TH2D("hist","Event type",4,0,4,4,0,4);
    TH2D * hists[4][2] = {0};
    for (int iTrig = 0; iTrig<4; iTrig++){
        hists[iTrig][0] = new TH2D(Form("hist%d_0",iTrig),"Event type",4,0,4,4,0,4);
        hists[iTrig][1] = new TH2D(Form("hist%d_1",iTrig),"Event type",4,0,4,4,0,4);
    }
    TH2D * hist_trig = new TH2D("hist_trig","Trigger type",4,0,4,4,0,4);

    const char * labelx[4] = {"none","pre-T0","T0","pre-T0&T0"};
    const char * labely[4] = {"none","T1","T2","T1&T2"};

    for (int x = 0; x<4; x++){
        for (int y = 0; y<4; y++){
            hist->Fill(labelx[x],labely[y],1);
            hist_trig->Fill(labelx[x],labely[y],1);
            for (int iTrig = 0; iTrig<4; iTrig++){
                hists[iTrig][0]->Fill(labelx[x],labely[y],1);
                hists[iTrig][1]->Fill(labelx[x],labely[y],1);
            }
        }
    }
    hist->Reset();
    hist_trig->Reset();
    for (int iTrig = 0; iTrig<4; iTrig++){
        hists[iTrig][0]->Reset();
        hists[iTrig][1]->Reset();
    }

    for (Long64_t entry = 0; entry<t->GetEntries(); entry++){
        t->GetEntry(entry);
        if (!foundTrigger) continue;
        int xtrig = (trig_np>0)+(trig_n0>0);
        int ytrig = (trig_h[2]>0)+(trig_h[3]>0);
        if (xtrig==1){
            if (trig_np>0) xtrig = 1;
            else xtrig = 2;
        }
        else if (xtrig==2){
            xtrig = 3;
        }
        if (ytrig==1){
            if (trig_h[2]>0) ytrig = 1;
            else ytrig = 2;
        }
        else if (ytrig==2){
            ytrig = 3;
        }
        hist_trig->Fill(labelx[xtrig],labely[ytrig],1);
        for (int i = 0; i<evt_n->size(); i++){
            if (fabs(evt_t->at(i)-trig_t)<1e-6) continue;
            int xevt = (evt_np->at(i)>0)+(evt_n0->at(i)>0);
            int yevt = (evt_h[2]->at(i)>0)+(evt_h[3]->at(i)>0);
            if (xevt==1){
                if (evt_np->at(i)>0) xevt = 1;
                else xevt = 2;
            }
            else if (xevt==2){
                xevt = 3;
            }
            if (yevt==1){
                if (evt_h[2]->at(i)>0) yevt = 1;
                else yevt = 2;
            }
            else if (yevt==2){
                yevt = 3;
            }
            hist->Fill(labelx[xevt],labely[yevt],1);
            if (xtrig>=2){ // without pre-T0 requirement: T0 & preT0+T0
                hists[ytrig][0]->Fill(labelx[xevt],labely[yevt],1);
            }
            if (xtrig==3){ // with pre-T0 requirement: preT0+T0
                hists[ytrig][1]->Fill(labelx[xevt],labely[yevt],1);
            }
        }
    }

    gStyle->SetOptStat(0);
    TCanvas * canv = new TCanvas();
    canv->SetCanvasSize(1200,1200);
    hist->SetTitle(Form("Event type %s: %.0f (%.1e /trigger) events in entries with trig.",runname.Data(),hist->GetEntries(),hist->GetEntries()/hist_trig->GetEntries()));
    hist->Draw("COLTEXT");
    for (int x = 0; x<4; x++){
        for (int y = 0; y<4; y++){
            TLatex * text = new TLatex(x+0.2,y+0.3,Form("%.1e /en.",hist->GetBinContent(x+1,y+1)/hist_trig->GetEntries()));
            text->SetTextFont(42);
            text->SetTextSize(0.02);
            text->Draw("SAME");
        }
    }
    canv->SaveAs(Form("pictures/type/evt.%s.trigall.png",outputname.Data()));
    hist_trig->SetTitle(Form("Trigger type %s: %lld entries, %.0f (%.1f %%)  with trigger",runname.Data(),t->GetEntries(),hist_trig->GetEntries(),100*hist_trig->GetEntries()/t->GetEntries()));
    hist_trig->Draw("COLTEXT");
    for (int x = 0; x<4; x++){
        for (int y = 0; y<4; y++){
            TLatex * text = new TLatex(x+0.3,y+0.3,Form("%.2f %%",100*hist_trig->GetBinContent(x+1,y+1)/t->GetEntries()));
            text->SetTextFont(42);
            text->SetTextSize(0.02);
            text->Draw("SAME");
        }
    }
    canv->SaveAs(Form("pictures/type/trig.%s.png",outputname.Data()));
    for (int iTrig = 0; iTrig<4; iTrig++){
        for (int withPT0 = 0; withPT0<=1; withPT0++){
            if (hists[iTrig][withPT0]->GetEntries()==0) continue;
            TString trigTitle = "4-fold trig.";
            TString trigName = "trig4f";
            if (iTrig==0){
                trigTitle = "2-fold trig.";
                trigName = "trig2f";
            }
            else if (iTrig==1){
                trigTitle = "2-fold+T1 trig.";
                trigName = "trig3f1";
            }
            else if (iTrig==2){
                trigTitle = "2-fold+T2 trig.";
                trigName = "trig3f2";
            }
            if (withPT0){
                trigTitle+=", with pre-T0";
                trigName+="pT0";
            }
            hists[iTrig][withPT0]->SetTitle(Form("Event type %s: %.0f (%.1e /trigger) events in entries with %s",runname.Data(),hists[iTrig][withPT0]->GetEntries(),hists[iTrig][withPT0]->GetEntries()/hist_trig->GetBinContent(4,iTrig+1),trigTitle.Data()));
            hists[iTrig][withPT0]->Draw("COLTEXT");
            for (int x = 0; x<4; x++){
                for (int y = 0; y<4; y++){
                    TLatex * text = new TLatex(x+0.2,y+0.3,Form("%.1e /en.",hists[iTrig][withPT0]->GetBinContent(x+1,y+1)/hist_trig->GetBinContent(4,iTrig+1)));
                    text->SetTextFont(42);
                    text->SetTextSize(0.02);
                    text->Draw("SAME");
                }
            }
            canv->SaveAs(Form("pictures/type/evt.%s.%s.png",outputname.Data(),trigName.Data()));
        }
    }

    return 0;
}
