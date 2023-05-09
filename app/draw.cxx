#include <iostream>
#include <sstream>
#include <cstdlib>
#include <stdlib.h> /* atoi, atof */
#include <map>

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

#define NCH 6
int NCHmax = NCH;

void print_usage(char * prog_name){
    fprintf(stderr,"Usage %s [options] runname\n",prog_name);
    fprintf(stderr,"[options]\n");
};

int main(int argc, char ** argv){
    // about the drawing
    TString runname = "test";
    TString runDescription = "";
    TString trigDescription = "";
    TString evtDescription = "";
    TString trigCuts = "foundTrigger";
    TString evtCuts = "evt_n>0";
    TString drawContent = "single"; // coin, single, trig
    bool    drawDetails = false;
    double daqLT = 1; // seconds

    // about ADC
    char sizeChar = 'h'; // h, w, a

    // about time
    bool useLog4Time = false;
    double t_min = -1300;
    double t_max = 9500;
    double t_step = 16;
    double t_min_trig = 300;
    double t_max_trig = 400;
    double t_step_trig = 0.8;

    // options for fitting
    bool fitTime = false;
    bool useSmallPeak = false;
    double tminFit = t_min;
    double tmaxFit = t_max; // TODO: shall consider to make it customizable

    // Load options
    int    opt_result;
    std::stringstream stream;
    while((opt_result=getopt(argc,argv,"Q:PDF:E:S:t:T:L:s:d:r:c:o:h"))!=-1){
        stream.clear();
        if(optarg) stream.str(optarg);
        switch(opt_result){
            /* INPUTS */
            case 'Q':
                daqLT = atof(optarg);
                std::cerr<<"DAQ live time set to "<<daqLT<<" seconds "<<std::endl;
                break;
            case 'P':
                useSmallPeak = true;
                std::cerr<<"Will use a series of Gaussian functions to fit the baseline"<<std::endl;
                break;
            case 'D':
                drawDetails = true;
                std::cerr<<"Will draw distribution plots in addition to the time plot"<<std::endl;
                break;
            case 'F':
                fitTime = true;
                stream>>tminFit>>tmaxFit;
                std::cerr<<"Will fit the time distribution between "<<tminFit<<" and "<<tmaxFit<<" ns"<<std::endl;
                break;
            case 'T':
                stream>>t_min_trig>>t_max_trig>>t_step_trig;
                std::cerr<<"trigger time set to "<<t_min_trig<<" "<<t_max_trig<<" "<<t_step_trig<<" ns"<<std::endl;
                break;
            case 't':
                stream>>t_min>>t_max>>t_step;
                std::cerr<<"time set to "<<t_min<<" "<<t_max<<" "<<t_step<<" ns"<<std::endl;
                break;
            case 'L':
                useLog4Time = atoi(optarg);
                std::cerr<<"Use log scale for time plot? "<<(useLog4Time?"yes":"no")<<std::endl;
                break;
            case 'd':
                runDescription = optarg;
                std::cerr<<"Use runDescription \""<<runDescription<<"\""<<std::endl;
                break;
            case 'S':
                stream>>trigCuts>>evtCuts;
                std::cerr<<"Use cuts for trigger \""<<trigCuts<<"\""<<" and cuts for evtects \""<<evtCuts<<"\""<<std::endl;
                break;
            case 'E':
                drawContent = optarg;
                std::cerr<<"Use drawContent \""<<drawContent<<"\""<<std::endl;
                break;
            case 's':
                sizeChar = optarg[0];
                std::cerr<<"Use size char '"<<sizeChar<<"'"<<std::endl;
                break;
            case 'r':
                runname = optarg;
                std::cerr<<"Run name \""<<runname<<"\""<<std::endl;
                break;
            case 'c':
                trigDescription = optarg;
                std::cerr<<"choose entries with trigger \""<<trigDescription<<"\""<<std::endl;
                break;
            case 'o':
                evtDescription = optarg;
                std::cerr<<"draw evtect \""<<evtDescription<<"\""<<std::endl;
                break;
            case '?':
                fprintf(stderr,"Wrong option! optopt=%c, optarg=%s\n", optopt, optarg);
            case 'h':
            default:
                print_usage(argv[0]);
                return 1;
        }
    }
    int t_nbins = (t_max-t_min)/t_step;
    int t_nbins_trig = (t_max_trig-t_min_trig)/t_step_trig;

    TChain* t = new TChain("t");
    for (int i = optind; i<argc; i++){
        t->Add(argv[i]);
    }
    if(0 == t->GetNtrees()) return 1;

    int size_nbins[NCH];
    double size_min[NCH];
    double size_max[NCH];
    TString sizeName("");
    TString sizeTitle("");
    if (sizeChar == 'h'){
        sizeName = "height";
        sizeTitle = "height [mV]";
        for (int ch = 0; ch<NCHmax; ch++){
            size_nbins[ch] = 127;
            size_min[ch] = 0;
            size_max[ch] = 1000;
        }
    }
    else if (sizeChar == 'w'){
        sizeName = "width";
        sizeTitle = "width [#samples]";
        for (int ch = 0; ch<NCHmax; ch++){
            size_nbins[ch] = 200;
            size_min[ch] = 0;
            size_max[ch] = 200;
        }
    }
    else if (sizeChar == 'a'){
        sizeName = "area";
        sizeTitle = "area [pC]";
        for (int ch = 0; ch<NCHmax; ch++){
            size_nbins[ch] = 250;
            size_min[ch] = 0;
            size_max[ch] = 500;
        }
    }
    int colors[7] = {kRed,kBlue,kMagenta,kCyan,kGreen,kOrange,kGray};

    TString cuts="("+trigCuts+")&&("+evtCuts+")";
    TString cuts4trig = cuts;
    if (drawContent == "trig") cuts4trig = trigCuts;

    TCanvas * canv = 0;
    TPad * pads[NCH];
    double nEntries = t->GetEntries();
    double nGoodEntries = t->GetEntries("foundTrigger");
    double nChosenEntries = t->GetEntries(trigCuts);
    double nChosenEntriesWithEvent = t->GetEntries(cuts);
    TLatex * text = 0;
    TLatex * textTitle = new TLatex(0.03,0.93,"");
    textTitle->SetTextFont(42);
    textTitle->SetTextSize(0.025);
    gStyle->SetOptStat(0);
    TLegend * legend;
    TH1D * h1;
    TH2D * h2;
    double maximum = 0;
    double height;
    double nEvents;

    //=======================================================================
    // first, draw the time distribution
    if (drawContent!="trig"){
        canv = new TCanvas("time");
        canv->SetCanvasSize(1200,800);
        pads[0] = new TPad("pad0","",0,0,1,0.9);
        pads[0]->Draw();
        pads[0]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogy(useLog4Time);
        t->Draw(Form("evt_t-trig_t>>h_t(%d,%.1f,%.1f)",t_nbins,t_min,t_max),cuts,"HIST");
        h1 = (TH1D*)gDirectory->Get("h_t");
        h1->SetTitle(";Time w.r.t. trigger [ns];Count");
        h1->SetLineColor(kBlue);
        h1->SetFillColor(kBlue);
        h1->GetXaxis()->SetRangeUser(0,t_max);
        //h1->GetXaxis()->SetRangeUser(20,110); // for 4-fold, hit D to avoid the huge reflection peak at 120 ns
        //h1->GetXaxis()->SetRangeUser(tminFit,110); // for new runs to show the body
        height = h1->GetMaximum();
        h1->GetXaxis()->UnZoom();
        if (useLog4Time){
            h1->GetYaxis()->SetRangeUser(0.5,height*1.1);
        }
        else{
            h1->GetYaxis()->SetRangeUser(0,height*1.1);
        }
        nEvents = h1->GetEntries();
        if (nEvents==0){
            std::cout<<"Cannot find any entries!"<<std::endl;
            return 0;
        }
        double nEventsInFitWIndow = h1->Integral(h1->FindBin(tminFit),h1->FindBin(tmaxFit));
        textTitle->SetText(0.03,0.93,Form("#splitline{#splitline{%s, %.0f entries, DAQ L.T. %.0f sec}{%.2f%% entries w/ trigger, %.2f%% w/ %s, %.2e/sec}}{%s: %.0f (%.1e/trigger), %.0f (%.1e/trigger) in the fitting window %.0f ns ~ %.0f ns}",runDescription.Data(),nEntries,daqLT,nGoodEntries/nEntries*100,nChosenEntries/nEntries*100,trigDescription.Data(),nChosenEntries/daqLT,evtDescription.Data(),nEvents,nEvents/nChosenEntries,nEventsInFitWIndow,nEventsInFitWIndow/nChosenEntries,tminFit,tmaxFit));

        if (!fitTime){
            canv->cd();
            textTitle->Draw();
            canv->SaveAs(Form("fitT.%s.%c.png",runname.Data(),sizeChar));
        }
        //=======================================================================
        // if we need to fit
        else{
            TString formula="";
            TString coef;
            for (int i = -2;i<=16; i++){
                if (i==0) coef="[5]+[6]+[7]+[8]";
                //if (i==0) coef="[6]";
                else if ((i+9)%9==0) coef="[5]+[6]+[7]+[8]";
                else if ((i+9)%9==2) coef="[6]+[7]+[8]";
                else if ((i+9)%9==3) coef="[5]";
                else if ((i+9)%9==4) coef="[7]+[8]";
                else if ((i+9)%9==5) coef="[5]+[6]";
                else if ((i+9)%9==6) coef="[8]";
                else if ((i+9)%9==7) coef="[5]+[6]+[7]";
                else continue;
                if (i>-2) formula+="+";
                formula+=Form("(%s)*(exp(-0.5*((x-[0]-%.1f)/[1])^2)+exp(-0.5*((x-[2]-%.1f)/[3])^2)*[4])",coef.Data(),i*1170/2.,i*1170/2.);
            }
            TF1 * fbig = new TF1("fbig",formula,t_min,t_max);
            int nParBig = fbig->GetNpar();
            fbig->SetParameters(-20,40,-20,100,0.25,45,45,45,45);
            //fbig->SetParameters(27.5,45,200);
            fbig->SetNpx(t_nbins);
            fbig->SetLineColor(kBlack);
            fbig->SetLineStyle(2);
            TString formulaTotal=formula;
            TString formulaSubtracted=formula;

            formulaTotal += Form("+(x>0)*(expo(%d)+expo(%d))",nParBig,nParBig+2);
            formulaSubtracted += Form("+(x>0)*(expo(%d)+expo(%d))",nParBig,nParBig+2);
            TF1 * fdecay = new TF1("fdecay","(x>0)*(expo(0)+expo(2))",tminFit,tmaxFit);
            int nParDecay = fdecay->GetNpar();
            fdecay->SetParameters(3.44763e+00,-1.68350e-03,6.01586e+00,-3.11848e-02);
            fdecay->SetNpx(t_nbins);
            fdecay->SetLineColor(kRed);
            fdecay->SetLineStyle(2);

            formula="";
            if (useSmallPeak){
                for (int i = -2;i<=16; i++){
                    if (i>-2) formula+="+";
                    formula+=Form("[2]*exp(-0.5*((x-[0]-%.1f)/[1])^2)",i*1170/2.);
                    formulaTotal+=Form("+[%d]*exp(-0.5*((x-[%d]-%.1f)/[%d])^2)",nParBig+nParDecay+2,nParBig+nParDecay+0,i*1170/2.,nParBig+nParDecay+1);
                }
            }
            else{
                formula = "[0]";
                formulaTotal += Form("+[%d]",nParBig+nParDecay);
            }
            TF1 * fsmall = new TF1("fsmall",formula,t_min,t_max);
            if (useSmallPeak){
                fsmall->SetParameters(45,163,15);
            }
            fsmall->SetNpx(t_nbins);
            fsmall->SetLineColor(kGray+2);
            fsmall->SetLineStyle(2);

            TF1 * ftotal = new TF1("ftotal",formulaTotal,t_min,t_max);
            ftotal->SetNpx(t_nbins);
            ftotal->SetLineColor(kBlack);
            for (int i = 0; i<ftotal->GetNpar(); i++){
                if (i<nParBig){
                    ftotal->SetParameter(i,fbig->GetParameter(i));
                    ftotal->SetParLimits(i,0,1000);
                }
                else if (i<nParBig+nParDecay){
                    ftotal->SetParameter(i,fdecay->GetParameter(i-nParBig));
                    if ((i-nParBig)%2==0){ // height
                        ftotal->SetParLimits(i,0,100);
                    }
                    else{ // decay constant
                        ftotal->SetParLimits(i,-1,0);
                    }
                }
                else{
                    ftotal->SetParameter(i,fsmall->GetParameter(i-nParBig-nParDecay));
                    ftotal->SetParLimits(i,0,1000);
                }
                ftotal->SetParLimits(0,-50,50);
                ftotal->SetParLimits(1,0,50);
                ftotal->SetParLimits(2,-100,100);
                ftotal->SetParLimits(3,10,100);
                ftotal->SetParLimits(4,0,10);
                if (useSmallPeak){
                    ftotal->SetParLimits(nParBig+nParDecay,-50,50); // position for small peaks
                    ftotal->SetParLimits(nParBig+nParDecay+1,0,300); // sigma for small peaks
                }
                else{
                    ftotal->SetParLimits(nParBig+nParBig,0,height*0.5);
                }
            }

            h1->Fit(ftotal,"L","",tminFit,tmaxFit);
            double Tau1 = -1/ftotal->GetParameter(nParBig+1);
            double Tau2 = -1/ftotal->GetParameter(nParBig+3);
            double Tau1min = -1/(ftotal->GetParameter(nParBig+1)-ftotal->GetParError(nParBig+1));
            double Tau2min = -1/(ftotal->GetParameter(nParBig+3)-ftotal->GetParError(nParBig+3));
            double Tau1max = -1/(ftotal->GetParameter(nParBig+1)+ftotal->GetParError(nParBig+1));
            double Tau2max = -1/(ftotal->GetParameter(nParBig+3)+ftotal->GetParError(nParBig+3));
            double dTau1 = (Tau1max-Tau1min)/2;
            double dTau2 = (Tau2max-Tau2min)/2;
            bool fixedTau = false;
            if ((Tau1-dTau1*2<2194&&Tau1+dTau1*2>2194&&Tau2-dTau2*2<165 &&Tau2+dTau2*2>165)
              ||(Tau1-dTau1*2<165 &&Tau1+dTau1*2>165 &&Tau2-dTau2*2<2194&&Tau2+dTau2*2>2194)){
                ftotal->FixParameter(nParBig+1,-1/2195.);
                ftotal->FixParameter(nParBig+3,-1/165.);
                h1->Fit(ftotal,"L","",tminFit,tmaxFit);
                fixedTau = true;
            }

            for (int i = 0; i<ftotal->GetNpar(); i++){
                if (i<nParBig){
                    fbig->SetParameter(i-0,ftotal->GetParameter(i));
                    fbig->SetParError(i-0,ftotal->GetParError(i));
                }
                else if (i<nParBig+nParDecay){
                    fdecay->SetParameter(i-nParBig,ftotal->GetParameter(i));
                    fdecay->SetParError(i-nParBig,ftotal->GetParError(i));
                }
                else{
                    fsmall->SetParameter(i-nParBig-nParDecay,ftotal->GetParameter(i));
                    fsmall->SetParError(i-nParBig-nParDecay,ftotal->GetParError(i));
                }
            }

            /*
            for (int i = -2; i<=16; i++){
                double tstart= -250+1170./2*i;
                double tstop = tstart+1170./2;
                TLine * line = new TLine(tstart,0,tstart,height);
                std::cout<<i<<" "<<tstart<<std::endl;
                line->SetLineStyle(2);
                line->SetLineColor(kRed);
                line->Draw();
                double Ndecay = 0;
                if (tstart>0) Ndecay = fdecay->Integral(tstart,tstop)/t_step;
                else if (tstop>0) Ndecay = fdecay->Integral(0,tstop)/t_step;
                double Npeak = fbig->Integral(tstart,tstop)/t_step;
                double Nbase = fsmall->Integral(tstart,tstop)/t_step;
                int ileft = h1->FindBin(tstart);
                int iright = h1->FindBin(tstop);
                double Nhist = h1->Integral(ileft+1,iright);
                text = new TLatex(tstart,height,Form("#splitline{%.0f}{#splitline{%.0f}{%.0f}}",Ndecay,Npeak,Nbase));
                if (useLog4Time) text->SetY(exp(log(height)*0.9));
                text->SetTextColor(kGray+2);
                text->SetTextFont(42);
                text->SetTextSize(0.023);
                text->Draw();
                text = new TLatex(tstart,height*1.05,Form("%.0f",Nhist));
                if (useLog4Time) text->SetY(exp(log(height)*0.95));
                text->SetTextColor(kBlack);
                text->SetTextFont(42);
                text->SetTextSize(0.023);
                text->Draw();
            }
            text = new TLatex(t_max,height,"#splitline{Ndecay}{#splitline{Npeak}{Nbase}}");
            if (useLog4Time) text->SetY(exp(log(height)*0.9));
            text->SetTextColor(kGray+2);
            text->SetTextFont(42);
            text->SetTextSize(0.024);
            text->Draw();
            text = new TLatex(t_max,height*1.05,"Nhist");
            if (useLog4Time) text->SetY(exp(log(height)*0.95));
            text->SetTextColor(kBlack);
            text->SetTextFont(42);
            text->SetTextSize(0.024);
            text->Draw();
            */

            ftotal->Draw("SAME");
            fdecay->Draw("SAME");
            fbig->Draw("SAME");
            fsmall->Draw("SAME");

            canv->cd();

            TF1 * fexpo = new TF1("fexpo","expo",t_min,t_max);
            fexpo->SetNpx(t_nbins);
            fexpo->SetParameter(0,fdecay->GetParameter(0));
            fexpo->SetParameter(1,fdecay->GetParameter(1));
            double nDecay1 = fexpo->Integral(tminFit,tmaxFit)/t_step;
            fexpo->SetParameter(0,fdecay->GetParameter(0)+fdecay->GetParError(0));
            fexpo->SetParameter(1,fdecay->GetParameter(1)+fdecay->GetParError(1));
            double nDecay1Max = fexpo->Integral(tminFit,tmaxFit)/t_step-nDecay1;
            fexpo->SetParameter(0,fdecay->GetParameter(0)-fdecay->GetParError(0));
            fexpo->SetParameter(1,fdecay->GetParameter(1)-fdecay->GetParError(1));
            double nDecay1Min = fexpo->Integral(tminFit,tmaxFit)/t_step-nDecay1;
            double nDecay1Error = (nDecay1Max-nDecay1Min)/2;
            fexpo->SetParameter(0,fdecay->GetParameter(2));
            fexpo->SetParameter(1,fdecay->GetParameter(3));
            double nDecay2 = fexpo->Integral(tminFit,tmaxFit)/t_step;
            fexpo->SetParameter(0,fdecay->GetParameter(2)+fdecay->GetParError(2));
            fexpo->SetParameter(1,fdecay->GetParameter(3)+fdecay->GetParError(3));
            double nDecay2Max = fexpo->Integral(tminFit,tmaxFit)/t_step-nDecay2;
            fexpo->SetParameter(0,fdecay->GetParameter(2)-fdecay->GetParError(2));
            fexpo->SetParameter(1,fdecay->GetParameter(3)-fdecay->GetParError(3));
            double nDecay2Min = fexpo->Integral(tminFit,tmaxFit)/t_step-nDecay2;
            double nDecay2Error = (nDecay2Max-nDecay2Min)/2;

            double nPeak = fbig->Integral(tminFit,tmaxFit)/t_step;
            double nBase = fsmall->Integral(tminFit,tmaxFit)/t_step;
            
            TLatex * text_decay;
            if (fixedTau){
                if (Tau1<Tau2){
                    double temp;
                    temp = Tau2; Tau2=Tau1; Tau1=temp;
                    temp = dTau2; dTau2=dTau1; dTau1=temp;
                }
                text_decay = new TLatex(0.03,0.89,Form("Decay: #tau_{1}=%.0f#pm%.0f ns (set to 2195 ns), n_{1}=%.0f#pm%.0f (%.2e#pm%.2e s^{-1}), #tau_{2}=%.0f#pm%.0f ns (set to 165 ns), n_{2}=%.0f#pm%.0f (%.2e#pm%.2e s^{-1})",Tau1,dTau1,nDecay1,nDecay1Error,nDecay1/daqLT,nDecay1Error/daqLT,Tau2,dTau2,nDecay2,nDecay2Error,nDecay2/daqLT,nDecay2Error/daqLT));
            }
            else
                text_decay = new TLatex(0.03,0.89,Form("Decay: #tau_{1}=%.0f#pm%.0f ns, n_{1}=%.0f#pm%.0f (%.2e#pm%.2e s^{-1}), #tau_{2}=%.0f#pm%.0f ns, n_{2}=%.0f#pm%.0f (%.2e#pm%.2e s^{-1})",Tau1,dTau1,nDecay1,nDecay1Error,nDecay1/daqLT,nDecay1Error/daqLT,Tau2,dTau2,nDecay2,nDecay2Error,nDecay2/daqLT,nDecay2Error/daqLT));
            text_decay->SetTextColor(kRed);
            text_decay->SetTextFont(42);
            text_decay->SetTextSize(0.024);
            text_decay->Draw();
            TLatex * text_big = new TLatex(0.03,0.86,Form("Peaks: double Gaussian, #sigma_{1}=%.1f#pm%.1fns, t_{1}=%.1f#pm%.1fns, #sigma_{2}=%.1f#pm%.1fns, t_{2}=%.1f#pm%.1fns, n=%.0f (%.2e s^{-1})",fbig->GetParameter(1),fbig->GetParError(1),fbig->GetParameter(0),fbig->GetParError(0),fbig->GetParameter(3),fbig->GetParError(3),fbig->GetParameter(2),fbig->GetParError(2),nPeak,nPeak/daqLT));
            text_big->SetTextColor(kBlack);
            text_big->SetTextFont(42);
            text_big->SetTextSize(0.024);
            text_big->Draw();
            TLatex * text_small = 0;
            if (useSmallPeak){
                text_small = new TLatex(0.03,0.83,Form("Baseline: single Gaussian, #sigma=%.1f#pm%.1fns, t=%.1f#pm%.1fns, n=%.0f (%.2e s^{-1})",fsmall->GetParameter(1),fsmall->GetParError(1),fsmall->GetParameter(0),fsmall->GetParError(0),nBase,nBase/daqLT));
            }
            else{
                text_small = new TLatex(0.03,0.83,Form("Baseline: flat bottom, height %.1f#pm%.1f per %.0fns, n = %.0f (%.2e s^{-1})",fsmall->GetParameter(0),fsmall->GetParError(0),t_step,nBase,nBase/daqLT));
            }
            text_small->SetTextColor(kGray+2);
            text_small->SetTextFont(42);
            text_small->SetTextSize(0.024);
            text_small->Draw();
            textTitle->Draw();
            canv->SaveAs(Form("fitT.%s.png",runname.Data()));

            TH1D * h1Subtracted = new TH1D("h1Subtracted",";Time w.r.t. trigger [ns];Count",t_nbins,t_min,t_max);
            for (int i = 1; i<=t_nbins; i++){
                double content = h1->GetBinContent(i);
                double x = h1->GetBinCenter(i);
                //h1Subtracted->SetBinContent(i,content-fdecay->Eval(x));
                //h1Subtracted->SetBinContent(i,content-fsmall->Eval(x)-fdecay->Eval(x));
                double y = content-fsmall->Eval(x)-fbig->Eval(x);
                h1Subtracted->SetBinContent(i,y);
            }
            if (useLog4Time){
                h1Subtracted->GetYaxis()->SetRangeUser(0.5,height*1.1);
            }
            else{
                h1Subtracted->GetYaxis()->SetRangeUser(0,height*1.1);
            }
            h1Subtracted->SetLineColor(kCyan);
            h1Subtracted->SetFillColor(kCyan);

            pads[0]->cd();
            if (!useLog4Time){
                gPad->SetMargin(0.1,0.1,0,0.1);
                h1Subtracted->GetYaxis()->SetRangeUser(-gPad->GetUymax()/8.,gPad->GetUymax());
            }
            h1Subtracted->Draw();
            for (int i = -2; i<=16; i++){
                double tstart= -250+1170./2*i;
                double tstop = tstart+1170./2;
                TLine * line = new TLine(tstart,0,tstart,height);
                std::cout<<i<<" "<<tstart<<std::endl;
                line->SetLineStyle(2);
                line->SetLineColor(kRed);
                line->Draw();
                double Ndecay = 0;
                if (tstart>0) Ndecay = fdecay->Integral(tstart,tstop)/t_step;
                else if (tstop>0) Ndecay = fdecay->Integral(0,tstop)/t_step;
                double Npeak = fbig->Integral(tstart,tstop)/t_step;
                double Nbase = fsmall->Integral(tstart,tstop)/t_step;
                int ileft = h1Subtracted->FindBin(tstart);
                int iright = h1Subtracted->FindBin(tstop);
                double Nhist = h1Subtracted->Integral(ileft+1,iright);
                text = new TLatex(tstart,height,Form("%.0f",Ndecay));
                if (useLog4Time) text->SetY(exp(log(height)*0.9));
                text->SetTextColor(kGray+2);
                text->SetTextFont(42);
                text->SetTextSize(0.023);
                text->Draw();
                text = new TLatex(tstart,height*1.05,Form("%.0f",Nhist));
                if (useLog4Time) text->SetY(exp(log(height)*0.95));
                text->SetTextColor(kBlack);
                text->SetTextFont(42);
                text->SetTextSize(0.023);
                text->Draw();
            }
            text = new TLatex(t_max,height,"Ndecay");
            if (useLog4Time) text->SetY(exp(log(height)*0.9));
            text->SetTextColor(kGray+2);
            text->SetTextFont(42);
            text->SetTextSize(0.024);
            text->Draw();
            text = new TLatex(t_max,height*1.05,"Nhist");
            if (useLog4Time) text->SetY(exp(log(height)*0.95));
            text->SetTextColor(kBlack);
            text->SetTextFont(42);
            text->SetTextSize(0.024);
            text->Draw();

            fdecay->Draw("SAME");

            canv->SaveAs(Form("fitT.%s.subtract.png",runname.Data()));
        }
    }

    if (!drawDetails) return 0;

    //=======================================================================
    // second, draw the trigger distribution plots
    canv = new TCanvas("trigdist");
    canv->SetCanvasSize(1200,1200);
    pads[0] = new TPad("pad0","",0,0,0.5,0.9);
    pads[1] = new TPad("pad1","",0.5,0.45,1,0.9);
    pads[2] = new TPad("pad2","",0.5,0,1,0.45);
    for (int i = 0; i<3; i++){
        pads[i]->Draw(); pads[i]->SetGridx(1); pads[i]->SetGridy(1); pads[i]->SetTickx(1); pads[i]->SetTicky(1);
    }
    pads[0]->cd(); pads[0]->SetLogy(1);
    legend = new TLegend(0.5,0.7,0.9,0.9);
    maximum = 0;
    for (int ch = 0; ch<NCHmax; ch++){
        TString opt = "HISTSAME";
        if (ch == 0) opt = "HIST";
        t->Draw(Form("trig_%c%c>>h_s_%c(%d,%.16e,%.16e",sizeChar,'A'+ch,'A'+ch,size_nbins[0],size_min[0],size_max[0]),cuts4trig+Form("&&trig_h%c>0",'A'+ch),opt);
        h1 = (TH1D*)gDirectory->Get(Form("h_s_%c",'A'+ch));
        h1->SetLineColor(colors[ch]);
        h1->SetTitle(Form(";%s;Count",sizeTitle.Data()));
        if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
        legend->AddEntry(h1,Form("ch%c: %.1f, %.0f",'A'+ch,h1->GetMean(),h1->GetEntries()));
    }
    legend->SetTextSize(0.022);
    h1 = (TH1D*)gDirectory->Get("h_s_A");
    h1->GetYaxis()->SetRangeUser(0.5,maximum*2);
    legend->Draw("SAME");
    // size VS time
    pads[1]->cd(); pads[1]->SetLogz(0);
    t->Draw(Form("trig_%c:trig_t>>h_st(%d,%.1f,%.1f,%d,%.16e,%.16e)",sizeChar,t_nbins_trig,t_min_trig,t_max_trig,size_nbins[0],size_min[0],size_max[0]),cuts4trig,"COLZ");
    h2 = (TH2D*)gDirectory->Get("h_st");
    h2->SetTitle(Form("Average size %.1f, %.0f;Time [ns];%s",h2->GetMean(2),h2->GetEntries(),sizeTitle.Data()));
    // time
    pads[2]->cd(); pads[2]->SetLogy(1);
    t->Draw(Form("trig_t>>h_t(%d,%.1f,%.1f)",t_nbins_trig,t_min_trig,t_max_trig),cuts4trig,"HIST");
    h1 = (TH1D*)gDirectory->Get("h_t");
    h1->SetTitle(";Time [ns];Count");
    canv->cd();
    if (drawContent=="trig"){
        textTitle->SetText(0.03,0.95,Form("#splitline{%s, %.0f entries, %.0f sec}{%.2f%% entries w/ trigger, %.2f%% w/ %s, %.2e/sec}",runDescription.Data(),nEntries,daqLT,nGoodEntries/nEntries*100,nChosenEntries/nEntries*100,trigDescription.Data(),nChosenEntries/daqLT));
    }
    else{
        textTitle->SetText(0.03,0.93,Form("#splitline{#splitline{%s, %.0f entries, %.0f sec}{%.2f%% entries w/ trigger, %.2f%% w/ %s, %.2e/sec}}{%s: %.0f, %.1e/trigger, %.2e/sec}",runDescription.Data(),nEntries,daqLT,nGoodEntries/nEntries*100,nChosenEntries/nEntries*100,trigDescription.Data(),nChosenEntries/daqLT,evtDescription.Data(),nEvents,nEvents/nChosenEntries,nEvents/daqLT));
    }
    textTitle->Draw();
    canv->SaveAs(Form("trigdist.%s.%c.png",runname.Data(),sizeChar));

    //      and the trigger coincidence plots
    canv = new TCanvas("trigcoin");
    canv->SetCanvasSize(1200,1200);
    for (int i = 0; i<2; i++){
        for (int j = 0; j<2; j++){
            int index = i*2+j;
            pads[index] = new TPad(Form("pad%d_%d",i,j),"",0.5*j,0.45*(1-i),0.5*(j+1),0.45*(2-i));
            pads[index]->Draw();
            pads[index]->SetGridx(1); pads[index]->SetGridy(1); pads[index]->SetTickx(1); pads[index]->SetTicky(1);
        }
    }
    // delta t
    pads[0]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogy(1);
    legend = new TLegend(0.6,0.7,1.0,0.9);
    t->Draw("trig_tp-trig_t0>>h_dt_p00(51,-20.4,20.4)",cuts4trig+"&&trig_n0>0&&trig_np>0","HIST");
    h1 = (TH1D*) gDirectory->Get("h_dt_p00");
    h1->SetLineColor(kBlack);
    h1->SetTitle(";#Deltat [ns];Count");
    legend->AddEntry(h1,Form("pre-T_{0}-T_{0} %.1f (%.1f) ns, %.0f",h1->GetMean(),h1->GetRMS(),h1->GetEntries()));
    maximum = 0;
    if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
    t->Draw("trig_tD-trig_tC>>h_dt_21(51,-20.4,20.4)",cuts4trig+"&&trig_hC>0&&trig_hD>0","SAME");
    h1 = (TH1D*) gDirectory->Get("h_dt_21");
    h1->SetLineColor(kBlue);
    legend->AddEntry(h1,Form("T_{2}-T_{1} %.1f (%.1f) ns, %.0f",h1->GetMean(),h1->GetRMS(),h1->GetEntries()));
    if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
    t->Draw("trig_tC-trig_t0>>h_dt_10(51,-20.4,20.4)",cuts4trig+"&&trig_n0>0&&trig_hC>0","SAME");
    h1 = (TH1D*) gDirectory->Get("h_dt_10");
    h1->SetLineColor(kRed);
    legend->AddEntry(h1,Form("T_{1}-T_{0} %.1f (%.1f) ns, %.0f",h1->GetMean(),h1->GetRMS(),h1->GetEntries()));
    if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
    t->Draw("trig_tB-trig_tA>>h_dt_0(51,-20.4,20.4)",cuts4trig+"&&trig_hB>0&&trig_hA>0","SAME");
    h1 = (TH1D*) gDirectory->Get("h_dt_0");
    h1->SetLineColor(kGreen); h1->SetLineStyle(2);
    legend->AddEntry(h1,Form("T_{0}R-T_{0}L %.1f (%.1f) ns, %.0f",h1->GetMean(),h1->GetRMS(),h1->GetEntries()));
    if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
    t->Draw("trig_tF-trig_tE>>h_dt_p0(51,-20.4,20.4)",cuts4trig+"&&trig_hE>0&&trig_hF>0","SAME");
    h1 = (TH1D*) gDirectory->Get("h_dt_p0");
    h1->SetLineColor(kCyan); h1->SetLineStyle(2);
    legend->AddEntry(h1,Form("pre-T_{0}R-pre-T_{0}L %.1f (%.1f) ns, %.0f",h1->GetMean(),h1->GetRMS(),h1->GetEntries()));
    if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
    h1 = (TH1D*) gDirectory->Get("h_dt_p00");
    h1->GetYaxis()->SetRangeUser(0.5,maximum*5);
    legend->Draw("SAME");
    // 2-D size
    pads[1]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogz(1);
    t->Draw(Form("trig_%cC:trig_%c0>>h_s_10(%d,%.16e,%.16e,%d,%.16e,%.16e)",sizeChar,sizeChar,size_nbins[0],size_min[0],size_max[0],size_nbins[2],size_min[2],size_max[2]),cuts4trig+"&&trig_n0>0&&trig_hC>0","COLZ");
    pads[2]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogz(1);
    t->Draw(Form("trig_%cp:trig_%c0>>h_s_p00(%d,%.16e,%.16e,%d,%.16e,%.16e)",sizeChar,sizeChar,size_nbins[0],size_min[0],size_max[0],size_nbins[0],size_min[0],size_max[0]),cuts4trig+"&&trig_n0>0&&trig_np>0","COLZ");
    pads[3]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogz(1);
    t->Draw(Form("trig_%cD:trig_%cC>>h_s_21(%d,%.16e,%.16e,%d,%.16e,%.16e)",sizeChar,sizeChar,size_nbins[2],size_min[2],size_max[2],size_nbins[2],size_min[2],size_max[2]),cuts4trig+"&&trig_hC>0&&trig_hD>0","COLZ");
    h2 = (TH2D*) gDirectory->Get("h_s_10");
    h2->GetYaxis()->SetMaxDigits(2);
    h2->SetTitle(Form("T0 %.1f T1 %.1f, %.0f;T0 %s;T1 %s",h2->GetMean(1),h2->GetMean(2),h2->GetEntries(),sizeTitle.Data(),sizeTitle.Data()));
    h2 = (TH2D*) gDirectory->Get("h_s_p00");
    h2->GetYaxis()->SetMaxDigits(2);
    h2->SetTitle(Form("T0 %.1f pre-T0 %.1f, %.0f;T0 %s;pre-T0 %s",h2->GetMean(1),h2->GetMean(2),h2->GetEntries(),sizeTitle.Data(),sizeTitle.Data()));
    h2 = (TH2D*) gDirectory->Get("h_s_21");
    h2->GetYaxis()->SetMaxDigits(2);
    h2->SetTitle(Form("T1 %.1f T2 %.1f, %.0f;T1 %s;T2 %s",h2->GetMean(1),h2->GetMean(2),h2->GetEntries(),sizeTitle.Data(),sizeTitle.Data()));
    canv->cd();
    textTitle->Draw();
    canv->SaveAs(Form("trigcoin.%s.%c.png",runname.Data(),sizeChar));

    if (drawContent=="trig") return 0;

    //=======================================================================
    // third, draw the event distribution plots
    canv = new TCanvas("dist");
    canv->SetCanvasSize(1200,1200);
    // size
    pads[0] = new TPad("pad0","",0,0.45,0.5,0.9);
    pads[1] = new TPad("pad1","",0.5,0.45,1,0.9);
    pads[2] = new TPad("pad2","",0,0,0.5,0.45);
    pads[3] = new TPad("pad3","",0.5,0,1,0.45);
    for (int i = 0; i<NCHmax; i++){
        pads[i]->Draw(); pads[i]->SetGridx(1); pads[i]->SetGridy(1); pads[i]->SetTickx(1); pads[i]->SetTicky(1);
    }
    pads[0]->cd(); pads[0]->SetLogy(1);
    legend = new TLegend(0.5,0.7,0.9,0.9);
    maximum = 0;
    for (int ch = 0; ch<NCHmax; ch++){
        // here skip channels if needed
        // ...
        TString opt = "HISTSAME";
        if (ch == 0) opt = "HIST";
        t->Draw(Form("evt_%c%c>>h_s_%c(%d,%.16e,%.16e",sizeChar,'A'+ch,'A'+ch,size_nbins[0],size_min[0],size_max[0]),cuts+Form("&&evt_h%c>0",'A'+ch),opt);
        h1 = (TH1D*)gDirectory->Get(Form("h_s_%c",'A'+ch));
        h1->SetLineColor(colors[ch]);
        h1->SetTitle(Form(";%s;Count",sizeTitle.Data()));
        if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
        legend->AddEntry(h1,Form("ch%c: %.1f, %.0f",'A'+ch,h1->GetMean(),h1->GetEntries()));
    }
    legend->SetTextSize(0.03);
    h1 = (TH1D*)gDirectory->Get("h_s_A");
    h1->GetYaxis()->SetRangeUser(0.5,maximum*2);
    legend->Draw("SAME");
    // size VS time
    pads[1]->cd(); pads[1]->SetLogz(0);
    t->Draw(Form("evt_%c:evt_t-trig_t>>h_st(%d,%.1f,%.1f,%d,%.16e,%.16e)",sizeChar,t_nbins,t_min,t_max,size_nbins[0],size_min[0],size_max[0]),cuts,"COLZ");
    h2 = (TH2D*)gDirectory->Get("h_st");
    h2->GetYaxis()->SetMaxDigits(2);
    h2->SetTitle(Form("Average size %.1f, %.0f;Time w.r.t. trigger [ns];%s",h2->GetMean(2),h2->GetEntries(),sizeTitle.Data()));
    // trigger size VS size
    pads[2]->cd(); pads[2]->SetLogz(0);
    t->Draw(Form("trig_%c:evt_%c>>h_tss(%d,%.1f,%.1f,%d,%.16e,%.16e)",sizeChar,sizeChar,size_nbins[0],size_min[0],size_max[0],size_nbins[0],size_min[0],size_max[0]),cuts,"COLZ");
    h2 = (TH2D*)gDirectory->Get("h_tss");
    h2->GetYaxis()->SetMaxDigits(2);
    h2->SetTitle(Form("Average event size %.1f, average trigger size %.1f, %.0f;%s;trigger %s",h2->GetMean(1),h2->GetMean(2),h2->GetEntries(),sizeTitle.Data(),sizeTitle.Data()));
    // trigger size VS time
    pads[3]->cd(); pads[3]->SetLogz(0);
    t->Draw(Form("trig_%c:evt_t-trig_t>>h_tst(%d,%.1f,%.1f,%d,%.16e,%.16e)",sizeChar,t_nbins,t_min,t_max,size_nbins[0],size_min[0],size_max[0]),cuts,"COLZ");
    h2 = (TH2D*)gDirectory->Get("h_tst");
    h2->GetYaxis()->SetMaxDigits(2);
    h2->SetTitle(Form("Average trigger size %.1f, %.0f;Time w.r.t. trigger [ns];trigger %s",h2->GetMean(2),h2->GetEntries(),sizeTitle.Data()));
    canv->cd();
    textTitle->Draw();
    canv->SaveAs(Form("dist.%s.%c.png",runname.Data(),sizeChar));

    //  and the coincidence plots
    if (drawContent!="single"){
        canv = new TCanvas("coin");
        canv->SetCanvasSize(1200,1200);
        for (int i = 0; i<2; i++){
            for (int j = 0; j<2; j++){
                int index = i*2+j;
                pads[index] = new TPad(Form("pad%d_%d",i,j),"",0.5*j,0.45*(1-i),0.5*(j+1),0.45*(2-i));
                pads[index]->Draw();
                pads[index]->SetGridx(1); pads[index]->SetGridy(1); pads[index]->SetTickx(1); pads[index]->SetTicky(1);
            }
        }
        // delta t
        pads[0]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogy(1);
        legend = new TLegend(0.6,0.7,1.0,0.9);
        t->Draw("evt_tp-evt_t0>>h_dt_p00(51,-20.4,20.4)",cuts+"&&evt_n0>0&&evt_np>0","HIST");
        h1 = (TH1D*) gDirectory->Get("h_dt_p00");
        h1->SetLineColor(kBlack);
        h1->SetTitle(";#Deltat [ns];Count");
        legend->AddEntry(h1,Form("pre-T_{0}-T_{0} %.1f (%.1f) ns, %.0f",h1->GetMean(),h1->GetRMS(),h1->GetEntries()));
        maximum = 0;
        if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
        t->Draw("evt_tD-evt_tC>>h_dt_21(51,-20.4,20.4)",cuts+"&&evt_hC>0&&evt_hD>0","SAME");
        h1 = (TH1D*) gDirectory->Get("h_dt_21");
        h1->SetLineColor(kBlue);
        legend->AddEntry(h1,Form("T_{2}-T_{1} %.1f (%.1f) ns, %.0f",h1->GetMean(),h1->GetRMS(),h1->GetEntries()));
        if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
        t->Draw("evt_tC-evt_t0>>h_dt_10(51,-20.4,20.4)",cuts+"&&evt_n0>0&&evt_hC>0","SAME");
        h1 = (TH1D*) gDirectory->Get("h_dt_10");
        h1->SetLineColor(kRed);
        legend->AddEntry(h1,Form("T_{1}-T_{0} %.1f (%.1f) ns, %.0f",h1->GetMean(),h1->GetRMS(),h1->GetEntries()));
        if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
        t->Draw("evt_tB-evt_tA>>h_dt_0(51,-20.4,20.4)",cuts+"&&evt_hB>0&&evt_hA>0","SAME");
        h1 = (TH1D*) gDirectory->Get("h_dt_0");
        h1->SetLineColor(kGreen); h1->SetLineStyle(2);
        legend->AddEntry(h1,Form("T_{0}R-T_{0}L %.1f (%.1f) ns, %.0f",h1->GetMean(),h1->GetRMS(),h1->GetEntries()));
        if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
        t->Draw("evt_tF-evt_tE>>h_dt_p0(51,-20.4,20.4)",cuts+"&&evt_hE>0&&evt_hF>0","SAME");
        h1 = (TH1D*) gDirectory->Get("h_dt_p0");
        h1->SetLineColor(kCyan); h1->SetLineStyle(2);
        legend->AddEntry(h1,Form("pre-T_{0}R-pre-T_{0}L %.1f (%.1f) ns, %.0f",h1->GetMean(),h1->GetRMS(),h1->GetEntries()));
        if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
        h1 = (TH1D*) gDirectory->Get("h_dt_p00");
        h1->GetYaxis()->SetRangeUser(0.5,maximum*5);
        legend->Draw("SAME");
        pads[1]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogz(1);
        t->Draw(Form("evt_%cC:evt_%c0>>h_s_10(%d,%.16e,%.16e,%d,%.16e,%.16e)",sizeChar,sizeChar,size_nbins[0],size_min[0],size_max[0],size_nbins[2],size_min[2],size_max[2]),cuts+"&&evt_n0>0&&evt_hC>0","COLZ");
        pads[2]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogz(1);
        t->Draw(Form("evt_%cp:evt_%c0>>h_s_p00(%d,%.16e,%.16e,%d,%.16e,%.16e)",sizeChar,sizeChar,size_nbins[0],size_min[0],size_max[0],size_nbins[0],size_min[0],size_max[0]),cuts+"&&evt_n0>0&&evt_np>0","COLZ");
        pads[3]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogz(1);
        t->Draw(Form("evt_%cD:evt_%cC>>h_s_21(%d,%.16e,%.16e,%d,%.16e,%.16e)",sizeChar,sizeChar,size_nbins[2],size_min[2],size_max[2],size_nbins[2],size_min[2],size_max[2]),cuts+"&&evt_hC>0&&evt_hD>0","COLZ");
        h2 = (TH2D*) gDirectory->Get("h_s_p00");
        h2->GetYaxis()->SetMaxDigits(2);
        h2->SetTitle(Form("T0 %.1f pre-T0 %.1f, %.0f;T0 %s;pre-T0 %s",h2->GetMean(1),h2->GetMean(2),h2->GetEntries(),sizeTitle.Data(),sizeTitle.Data()));
        h2 = (TH2D*) gDirectory->Get("h_s_10");
        h2->GetYaxis()->SetMaxDigits(2);
        h2->SetTitle(Form("T0 %.1f T1 %.1f, %.0f;T0 %s;T1 %s",h2->GetMean(1),h2->GetMean(2),h2->GetEntries(),sizeTitle.Data(),sizeTitle.Data()));
        h2 = (TH2D*) gDirectory->Get("h_s_21");
        h2->GetYaxis()->SetMaxDigits(2);
        h2->SetTitle(Form("T1 %.1f T2 %.1f, %.0f;T1 %s;T2 %s",h2->GetMean(1),h2->GetMean(2),h2->GetEntries(),sizeTitle.Data(),sizeTitle.Data()));
        canv->cd();
        textTitle->Draw();
        canv->SaveAs(Form("coin.%s.%c.png",runname.Data(),sizeChar));
    }

    return 0;
}
