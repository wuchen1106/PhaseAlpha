#include <iostream>
#include <stdlib.h> /* atoi, atof */

#include "TStyle.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TString.h"

double ADCstep[4] = {7.8277888298034667969,7.8277888298034667969,3.9138944149017333984,3.9138944149017333984};
double ADCoffset[4] = {0,0,0,0};

bool loadCalibADC(const TString & filename, int theRun){
    TChain * ichain = new TChain("t");
    ichain->Add(filename);
    if (ichain->GetEntries()==0) return true;
    double step[4];
    double offset[4];
    int run;
    ichain->SetBranchAddress("run",&run);
    for (int ch = 0; ch<4; ch++){
        ichain->SetBranchAddress(Form("step%c",'A'+ch),&step[ch]);
        ichain->SetBranchAddress(Form("offset%c",'A'+ch),&offset[ch]);
    }
    for (int e = 0; e<ichain->GetEntries(); e++){
        ichain->GetEntry(e);
        if (run!=theRun) continue;
        for (int ch = 0; ch<4; ch++){
            ADCstep[ch] = step[ch];
            ADCoffset[ch] = offset[ch];
        }
    }
    return false;
};

int main(int argc, char ** argv){
    TString input = "anaout/run01802-01881.thr20mV20mV10mV10mV.hit.root";
    TString runname = "run01802-01881.thr20mV20mV10mV10mV";
    TString sizeName = "width";
    if (argc > 1) input = argv[1];
    if (argc > 2) runname = argv[2];
    if (argc > 3) sizeName = argv[3];

    int threshold = 20;
    TString sizeUnit = Form("samples over %d mV",threshold);
    int nbins_size = 150;
    double left_size = -0.5;
    double right_size = nbins_size+left_size;
    bool isWidthV = false;
    if (sizeName[0] == 'a'){
        nbins_size = 1000;
        left_size = -0.5;
        right_size = nbins_size+left_size;
        sizeUnit = "pC";
    }
    else if (sizeName[0] == 'w' &&sizeName.Length()>5){
        threshold = atoi(sizeName(5,3).Data());
        sizeUnit = Form("samples over %d mV",threshold);
        isWidthV = true;
    }
    else if (sizeName == "fwhm"){
        sizeUnit = "samples over HM";
    }

    TChain * t = new TChain("t");
    t->Add(input);
    int run = 0;
    t->SetBranchAddress("run",&run);
    t->GetEntry(0);
    loadCalibADC("calib/ADCstep.root",run);

    gStyle->SetOptStat(0);
    TCanvas * canv = new TCanvas("canv","",1000,1000);
    canv->Divide(2,2);
    for (int ch = 0; ch<4; ch++){
        int nbins_height = 250;
        double left_height = -0.5*ADCstep[ch]+ADCoffset[ch];
        double right_height = nbins_height*ADCstep[ch]+left_height;
        canv->cd(ch+1); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(0); gPad->SetLogz(1);
        if (isWidthV){
            t->Draw(Form("width%c%d:height%c>>h%d(%d,%.16e,%.16e,%d,%.16e,%.16e)",'A'+ch,threshold,'A'+ch,ch,nbins_height,left_height,right_height,nbins_size,left_size,right_size),"","COL");
        }
        else{
            t->Draw(Form("%s%c:height%c>>h%d(%d,%.16e,%.16e,%d,%.16e,%.16e)",sizeName.Data(),'A'+ch,'A'+ch,ch,nbins_height,left_height,right_height,nbins_size,left_size,right_size),"","COL");
        }
        TH2D * h2 = (TH2D*) gDirectory->Get(Form("h%d",ch));
        h2->SetTitle(Form("ch%c;height [mV];%s [%s]",'A'+ch,sizeName.Data(),sizeUnit.Data()));
        TGraph * gr1 = new TGraph();
        TGraph * gr2 = new TGraph();
        int step = 3;
        int counter = 0;
        for (int i = 1; i<=h2->GetXaxis()->GetNbins(); i+=step){
            TH1D * h1 = h2->ProjectionY(Form("h%d_%d",ch,i),i,i+step-1);
            if (h1->GetEntries()<10) continue;
            double mean = h1->GetMean();
            double rms = h1->GetRMS();
            double x = 0;
            if (step%2==0){
                x = h2->GetXaxis()->GetBinLowEdge(i+step/2);
            }
            else{
                x = h2->GetXaxis()->GetBinCenter(i+step/2);
            }
            std::cout<<i<<" ~ "<<i+step-1<<" x "<<x<<std::endl;
            gr1->SetPoint(counter,x,mean);
            gr2->SetPoint(counter,x,(mean==0?0:rms/mean)/0.2*(right_size-left_size)+left_size);
            counter++;
        }
        gr1->SetLineColor(kBlack);
        gr2->SetLineColor(kRed);
        gr1->Draw("PLSAME");
        gr2->Draw("PLSAME");
        TGaxis * axis = new TGaxis(right_height,left_size,right_height,right_size,0,20,510,"+L");
        axis->SetTitleFont(42);
        axis->SetLabelFont(42);
        axis->SetTitle("RMS/Mean [%]");
        axis->Draw("SAME");
    }

    canv->SaveAs(Form("%s.%s.png",runname.Data(),sizeName.Data()));

    return 0;
};
