{
    TChain * t = new TChain("t","t"); 

    bool useSmallPeak = true;
    double tmax = 9500;
    double tmin = -1300;
    int nbins = 675;
    //double tmax = 9600;
    //double tmin = -400;
    //int nbins = 1250;
    double dt = (tmax-tmin)/nbins;
    double tminFit = tmin;
    double tmaxFit = tmax;

    bool useLog = false;
    gStyle->SetOptStat(0);
    TCanvas * canv0 = new TCanvas("canv0","",1024,768); canv0->SetMargin(0.1,0.1,0.1,0.13); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogy(useLog);
    TLatex * text = 0;

    t->Add("anaout/run00370-00505.thrA20mVB20mVC10mVD10mV.evtall.trigABbCbD90mV90mV70mV70mV.event.root");t->Add("anaout/run00584-00624.thrA20mVB20mVC10mVD10mV.evtall.trigABbCbD90mV90mV70mV70mV.event.root");t->Add("anaout/run00773-01000.thrA20mVB20mVC10mVD10mV.evtall.trigABbCbD90mV90mV70mV70mV.event.root"); TString run="run00370-...-01000"; tminFit = 30; tmaxFit = 9600;
    //t->Add("anaout/run01530-01599.thrA20mVB20mVC20mVD20mV.evtall.trigAB25mV25mV.event.root"); t->Add("anaout/run01600-01699.thrA20mVB20mVC20mVD20mV.evtall.trigAB25mV25mV.event.root"); t->Add("anaout/run01700-01798.thrA20mVB20mVC20mVD20mV.evtall.trigAB25mV25mV.event.root"); t->Add("anaout/run01802-01881.thrA20mVB20mVC20mVD20mV.evtall.trigAB25mV25mV.event.root"); TString run="run01530-01881"; tminFit = 60; tmaxFit = 8600; useSmallPeak = false; // use flat bottom instead

    int thrT[4] = {120,120,70,70};
    //TString cutT = Form("foundTrigger&&trig_hA>%d&&trig_hB>%d&&trig_n==2",thrT[0],thrT[1]); TString trigTitle = Form("hA>%d mV hB>%d mV, no hits in T1 T2",thrT[0],thrT[1]); TString trigName = Form("trig-Agt%dmV-Bgt%dmV-noCD",thrT[0],thrT[1]);
    //TString cutT = Form("foundTrigger&&trig_hA>%d&&trig_hB>%d&&trig_hC<%d&&trig_hD<%d",thrT[0],thrT[1],thrT[2],thrT[3]); TString trigTitle = Form("hA>%d mV hB>%d mV hC<%d mV hD<%d mV",thrT[0],thrT[1],thrT[2],thrT[3]); TString trigName = Form("trig-Agt%dmV-Bgt%dmV-Clt%dmV-Dlt%dmV",thrT[0],thrT[1],thrT[2],thrT[3]);
    //TString cutT = Form("foundTrigger&&trig_hA>%d&&trig_hB>%d",thrT[0],thrT[1]); TString trigTitle = Form("hA>%d mV hB>%d mV",thrT[0],thrT[1]); TString trigName = Form("trig-Agt%dmV-Bgt%dmV",thrT[0],thrT[1]);
    //TString cutT = Form("foundTrigger&&(trig_hA<%d||trig_hB<%d)",thrT[0],thrT[1]); TString trigTitle = Form("hA<%d mV or hB<%d mV",thrT[0],thrT[1]); TString trigName = Form("trig-Alt%dmV-or-Blt%dmV",thrT[0],thrT[1]);
    //TString cutT = Form("foundTrigger"); TString trigTitle = "all"; TString trigName = "trig-all";
    TString cutT = ("trig_n==2"); TString trigTitle = "w/ hits in T0L and T0R, w/o hits in T1/T2"; TString trigName = "trig-2f";

    int thrO[4] = {40,40,40,25};
    //int thrO[4] = {0,0,40,25};
    //int thrO[4] = {40,40,100,100};
    //int thrO[4] = {0,0,0,0};

    int ch1 = 3; int ch2 = 1;

    t->Draw(Form("evt_t%c-trig_t>>htime(%d,%.1f,%.1f)",'A'+ch1,nbins,tmin,tmax),cutT+"&&evt_n0==0&&evt_hC>30"); TH1D * htime = (TH1D*) gDirectory->Get("htime"); TString title=Form("#splitline{%s, %.0lld events, %.0lld with trigger: %s}{ %.0f single hits in %c: h>%d mV}",run.Data(),t->GetEntries(),t->GetEntries(cutT),trigTitle.Data(),htime->GetEntries(),'A'+ch1,thrO[ch1]); TString name=Form("single-%cgt%dmV",'A'+ch1,thrO[ch1]); // single-high

    // single hits
    //t->Draw(Form("evt_t%c-(trig_tA+trig_tB)/2>>htime(%d,%.1f,%.1f)",'A'+ch1,nbins,tmin,tmax),cutT+Form("&&evt_t%c!=trig_t%c&&evt_t%c!=peak110_t%c&&evt_t%c!=peak470_t%c&&evt_h%c<%d&&evt_h%c>0&&evt_n==1",'A'+ch1,'A'+ch1,'A'+ch1,'A'+ch1,'A'+ch1,'A'+ch1,'A'+ch1,thrO[ch1],'A'+ch1)); TH1D * htime = (TH1D*) gDirectory->Get("htime"); TString title=Form("#splitline{%s, %.0lld events, %.0lld with trigger: %s}{ %.0f single hits in %c: h<%d mV}",run.Data(),t->GetEntries(),t->GetEntries(cutT),trigTitle.Data(),htime->GetEntries(),'A'+ch1,thrO[ch1]); TString name=Form("single-%clt%dmV",'A'+ch1,thrO[ch1]); // single-low

    // pairs
    //t->Draw(Form("(evt_t%c+evt_t%c)/2-(trig_tA+trig_tB)/2>>htime(%d,%.1f,%.1f)",'A'+ch1,'A'+ch2,nbins,tmin,tmax),cutT+Form("&&evt_t%c!=trig_t%c&&evt_t%c!=peak110_t%c&&evt_t%c!=peak470_t%c&&evt_h%c>0&&evt_h%c>0&&(evt_h%c<%d||evt_h%c<%d)&&evt_n==2",'A'+ch1,'A'+ch1,'A'+ch1,'A'+ch1,'A'+ch1,'A'+ch1,'A'+ch1,'A'+ch2,'A'+ch1,thrO[ch1],'A'+ch2,thrO[ch2])); TH1D * htime = (TH1D*) gDirectory->Get("htime"); TString title=Form("#splitline{%s, %.0lld events, %.0lld with trigger: %s}{ %.0f pair hits in %c%c: h_{%c} <%d mV or h_{%c} <%d mV}",run.Data(),t->GetEntries(),t->GetEntries(cutT),trigTitle.Data(),htime->GetEntries(),'A'+ch1,'A'+ch2,'A'+ch1,thrO[ch1],'A'+ch2,thrO[ch2]); TString name=Form("pair-%clt%dmV-%clt%dmV",'A'+ch1,thrO[ch1],'A'+ch2,thrO[ch2]); // pair-low
    //t->Draw(Form("(evt_t%c+evt_t%c)/2-(trig_tA+trig_tB)/2>>htime(%d,%.1f,%.1f)",'A'+ch1,'A'+ch2,nbins,tmin,tmax),cutT+Form("&&evt_t%c!=trig_t%c&&evt_t%c!=peak110_t%c&&evt_t%c!=peak470_t%c&&evt_h%c>%d&&evt_h%c>%d&&evt_n==2",'A'+ch1,'A'+ch1,'A'+ch1,'A'+ch1,'A'+ch1,'A'+ch1,'A'+ch1,thrO[ch1],'A'+ch2,thrO[ch2])); TH1D * htime = (TH1D*) gDirectory->Get("htime"); TString title=Form("#splitline{%s, %.0lld events, %.0lld with trigger: %s}{ %.0f pair hits in %c%c: h_{%c} >%d mV h_{%c} >%d mV}",run.Data(),t->GetEntries(),t->GetEntries(cutT),trigTitle.Data(),htime->GetEntries(),'A'+ch1,'A'+ch2,'A'+ch1,thrO[ch1],'A'+ch2,thrO[ch2]); TString name=Form("pair-%cgt%dmV-%cgt%dmV",'A'+ch1,thrO[ch1],'A'+ch2,thrO[ch2]); // pair-high

    // 4-fold
    //t->Draw(Form("(evt_tA+evt_tB+evt_tC+evt_tD)/4-(trig_tA+trig_tB)/2>>htime(%d,%.1f,%.1f)",nbins,tmin,tmax),cutT+Form("&&evt_tA!=trig_tA&&evt_tA!=peak110_tA&&evt_tA!=peak470_tA&&(evt_hA<%d||evt_hB<%d||evt_hC<%d||evt_hD<%d)&&evt_n==4",thrO[0],thrO[1],thrO[2],thrO[3])); TH1D * htime = (TH1D*) gDirectory->Get("htime"); TString title=Form("#splitline{%s, %.0lld events, %.0lld with trigger: %s}{ %.0f 4-fold hits: h_{A} <%d mV or h_{B} <%d mV or h_{C} <%d mV or h_{D}<>%d mV}",run.Data(),t->GetEntries(),t->GetEntries(cutT),trigTitle.Data(),htime->GetEntries(),thrO[0],thrO[1],thrO[2],thrO[3]); TString name=Form("4fold-Alt%dmV-or-Blt%dmV-or-Clt%dmV-or-Dlt%dmV",thrO[0],thrO[1],thrO[2],thrO[3]); // 4-fold-low
    //t->Draw(Form("(evt_tA+evt_tB+evt_tC+evt_tD)/4-(trig_tA+trig_tB)/2>>htime(%d,%.1f,%.1f)",nbins,tmin,tmax),cutT+Form("&&evt_tA!=trig_tA&&evt_tA!=peak110_tA&&evt_tA!=peak470_tA&&evt_hA>%d&&evt_hB>%d&&evt_hC>%d&&evt_hD>%d&&evt_n==4",thrO[0],thrO[1],thrO[2],thrO[3])); TH1D * htime = (TH1D*) gDirectory->Get("htime"); TString title=Form("#splitline{%s, %.0lld events, %.0lld with trigger: %s}{ %.0f 4-fold hits: h_{A} >%d mV h_{B} >%d mV h_{C} >%d mV h_{D} >%d mV}",run.Data(),t->GetEntries(),t->GetEntries(cutT),trigTitle.Data(),htime->GetEntries(),thrO[0],thrO[1],thrO[2],thrO[3]); TString name=Form("4fold-Agt%dmV-Bgt%dmV-Cgt%dmV-Dgt%dmV",thrO[0],thrO[1],thrO[2],thrO[3]); // 4-fold-high

    htime->SetTitle(title+";#Deltat [ns];Count");
    htime->SetLineColor(kBlue);
    htime->SetFillColor(kBlue);
    htime->GetXaxis()->SetRangeUser(20,tmax);
    //htime->GetXaxis()->SetRangeUser(20,110); // for 4-fold, hit D to avoid the huge reflection peak at 120 ns
    double height = htime->GetMaximum();
    htime->GetXaxis()->UnZoom();
    if (useLog){
        htime->GetYaxis()->SetRangeUser(0.5,height*1.1);
    }
    else{
        htime->GetYaxis()->SetRangeUser(0,height*1.1);
    }

    TString formula="";
    TString coef;
    for (int i = -2;i<=16; i++){
        //if (i==0) coef="[5]+[6]+[7]+[8]";
        ////if (i==0) coef="[6]";
        //else if ((i+9)%9==0) coef="[5]+[6]+[7]+[8]";
        //else if ((i+9)%9==2) coef="[6]+[7]+[8]";
        //else if ((i+9)%9==3) coef="[5]";
        //else if ((i+9)%9==4) coef="[7]+[8]";
        //else if ((i+9)%9==5) coef="[5]+[6]";
        //else if ((i+9)%9==6) coef="[8]";
        //else if ((i+9)%9==7) coef="[5]+[6]+[7]";
        //else continue;
        if (i==0) coef="4*[5]";
        //if (i==0) coef="[6]";
        else if ((i+9)%9==0) coef="4*[5]";
        else if ((i+9)%9==2) coef="4*[5]";
        else if ((i+9)%9==3) coef="[5]";
        else if ((i+9)%9==4) coef="2*[5]";
        else if ((i+9)%9==5) coef="2*[5]";
        else if ((i+9)%9==6) coef="[5]";
        else if ((i+9)%9==7) coef="3*[5]";
        else continue;
        if (i>-2) formula+="+";
        formula+=Form("(%s)*(exp(-0.5*((x-[0]-%.1f)/[1])^2)+exp(-0.5*((x-[2]-%.1f)/[3])^2)*[4])",coef.Data(),i*1170/2.,i*1170/2.);
    }
    TF1 * fbig = new TF1("fbig",formula,tmin,tmax);
    int nParBig = fbig->GetNpar();
    //fbig->SetParameters(-20,20,-20,50,1,45,45,45,45);
    fbig->SetParameters(-20,20,-20,50,1,10);
    //fbig->SetParameters(27.5,45,200);
    fbig->SetNpx(nbins);
    fbig->SetLineColor(kBlack);
    fbig->SetLineStyle(2);
    TString formulaTotal=formula;
    TString formulaSubtracted=formula;

    formulaTotal += Form("+(x>0)*(expo(%d)+expo(%d))",nParBig,nParBig+2);
    formulaSubtracted += Form("+(x>0)*(expo(%d)+expo(%d))",nParBig,nParBig+2);
    TF1 * fdecay = new TF1("fdecay","(x>0)*(expo(0)+expo(2))",0,tmax);
    int nParDecay = fdecay->GetNpar();
    fdecay->SetParameters(3.44763e+00,-1.68350e-03,6.01586e+00,-3.11848e-02);
    fdecay->SetNpx(nbins);
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
    TF1 * fsmall = new TF1("fsmall",formula,tmin,tmax);
    if (useSmallPeak){
        fsmall->SetParameters(45,163,15);
    }
    fsmall->SetNpx(nbins);
    fsmall->SetLineColor(kGray+2);
    fsmall->SetLineStyle(2);

    TF1 * ftotal = new TF1("ftotal",formulaTotal,tmin,tmax);
    ftotal->SetNpx(nbins);
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
        ftotal->SetParLimits(0,-50,20);
        ftotal->SetParLimits(1,1,50);
        ftotal->SetParLimits(2,-50,20);
        ftotal->SetParLimits(3,1,100);
        ftotal->SetParLimits(4,0,10);
        if (useSmallPeak){
            ftotal->SetParLimits(nParBig+nParDecay,0,100); // position for small peaks
            ftotal->SetParLimits(nParBig+nParDecay+1,0,300); // sigma for small peaks
        }
        else{
            ftotal->SetParLimits(nParBig+nParBig,0,height*0.5);
        }
    }

    htime->Fit(ftotal,"L","",tminFit,tmaxFit);
    htime->SetTitle(title+";#Deltat [ns];Count");

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

    for (int i = -2; i<=16; i++){
        double tstart= -250+1170./2*i;
        double tstop = tstart+1170./2;
        //int ibin = htime->FindBin(tstart);
        //int jbin = htime->FindBin(tstop)-1;
        //double N0 = htime->Integral(ibin,jbin);
        TLine * line = new TLine(tstart,0,tstart,height);
        std::cout<<i<<" "<<tstart<<std::endl;
        line->SetLineStyle(2);
        line->SetLineColor(kRed);
        line->Draw();
        double Ndecay = 0;
        if (tstart>0) Ndecay = fdecay->Integral(tstart,tstop)/dt;
        else if (tstop>0) Ndecay = fdecay->Integral(0,tstop)/dt;
        double Nbig = fbig->Integral(tstart,tstop)/dt;
        double Nsmall = fsmall->Integral(tstart,tstop)/dt;
        int ileft = htime->FindBin(tstart);
        int iright = htime->FindBin(tstop);
        double Nhist = htime->Integral(ileft+1,iright);
        text = new TLatex(tstart,height,Form("#splitline{%.0f}{#splitline{%.0f}{%.0f}}",Ndecay,Nbig,Nsmall));
        if (useLog) text->SetY(exp(log(height)*0.9));
        text->SetTextColor(kGray+2);
        text->SetTextFont(42);
        text->SetTextSize(0.023);
        text->Draw();
        text = new TLatex(tstart,height*1.05,Form("%.0f",Nhist));
        if (useLog) text->SetY(exp(log(height)*0.95));
        text->SetTextColor(kBlack);
        text->SetTextFont(42);
        text->SetTextSize(0.023);
        text->Draw();
    }
    text = new TLatex(tmax,height,"#splitline{Ndecay}{#splitline{Nbig}{Nsmall}}");
    if (useLog) text->SetY(exp(log(height)*0.9));
    text->SetTextColor(kGray+2);
    text->SetTextFont(42);
    text->SetTextSize(0.024);
    text->Draw();
    text = new TLatex(tmax,height*1.05,"Nhist");
    if (useLog) text->SetY(exp(log(height)*0.95));
    text->SetTextColor(kBlack);
    text->SetTextFont(42);
    text->SetTextSize(0.024);
    text->Draw();

    TF1 * fexpo = new TF1("fexpo","expo",tmin,tmax);
    fexpo->SetNpx(nbins);
    fexpo->SetParameter(0,fdecay->GetParameter(0));
    fexpo->SetParameter(1,fdecay->GetParameter(1));
    double nDecay1 = fexpo->Integral(30,tmax)/dt;
    fexpo->SetParameter(0,fdecay->GetParameter(2));
    fexpo->SetParameter(1,fdecay->GetParameter(3));
    double nDecay2 = fexpo->Integral(30,tmax)/dt;
    TLatex * text_decay = new TLatex((tmax+tmin)/2.5,height*0.8,Form("Decay: #tau_{1}=%.0f#pm%.0fns, n_{1}=%.0f, #tau_{2}=%.0f#pm%.0fns, n_{2}=%.0f",-1/fdecay->GetParameter(1),1/fdecay->GetParameter(1)-1/(fdecay->GetParameter(1)+fdecay->GetParError(1)),nDecay1,-1/fdecay->GetParameter(3),1/fdecay->GetParameter(3)-1/(fdecay->GetParameter(3)+fdecay->GetParError(3)),nDecay2));
    if (useLog) text_decay->SetY(exp(log(height)*0.8));
    text_decay->SetTextColor(kRed);
    text_decay->SetTextFont(42);
    text_decay->SetTextSize(0.024);
    text_decay->Draw();
    TLatex * text_big = new TLatex((tmax+tmin)/2.5,height*0.75,Form("Big peaks: #sigma_{1}=%.0f#pm%.1fns, t_{1}=%.0f#pm%.1fns, #sigma_{2}=%.0f#pm%.1fns, t_{2}=%.0f#pm%.1fns",fbig->GetParameter(1),fbig->GetParError(1),fbig->GetParameter(0),fbig->GetParError(0),fbig->GetParameter(3),fbig->GetParError(3),fbig->GetParameter(2),fbig->GetParError(2)));
    if (useLog) text_big->SetY(exp(log(height)*0.75));
    text_big->SetTextColor(kBlack);
    text_big->SetTextFont(42);
    text_big->SetTextSize(0.024);
    text_big->Draw();
    TLatex * text_small = 0;
    if (useSmallPeak){
        text_small = new TLatex((tmax+tmin)/2.5,height*0.7,Form("Small peaks: #sigma=%.0f#pm%.1fns, t=%.0f#pm%.1fns",fsmall->GetParameter(1),fsmall->GetParError(1),fsmall->GetParameter(0),fsmall->GetParError(0)));
    }
    else{
        text_small = new TLatex((tmax+tmin)/2.5,height*0.7,Form("Flat bottom: height %.0f#pm%.1f per %.1fns",fsmall->GetParameter(0),fsmall->GetParError(0),dt));
    }
    if (useLog) text_small->SetY(exp(log(height)*0.7));
    text_small->SetTextColor(kGray+2);
    text_small->SetTextFont(42);
    text_small->SetTextSize(0.024);
    text_small->Draw();
    
    ftotal->Draw("SAME");
    fdecay->Draw("SAME");
    fbig->Draw("SAME");
    fsmall->Draw("SAME");

    TCanvas * canv1 = new TCanvas("canv1","",1024,768); canv1->SetMargin(0.1,0.1,0.1,0.13); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogy(useLog);
    TH1D * htimeSubtracted = new TH1D("htimeSubtracted","",nbins,tmin,tmax);
    for (int i = 1; i<=nbins; i++){
        double content = htime->GetBinContent(i);
        double x = htime->GetBinCenter(i);
        //htimeSubtracted->SetBinContent(i,content-fdecay->Eval(x));
        //htimeSubtracted->SetBinContent(i,content-fsmall->Eval(x)-fdecay->Eval(x));
        double y = content-fsmall->Eval(x)-fbig->Eval(x);
        if (y<0) y = 0;
        htimeSubtracted->SetBinContent(i,y);
    }
    if (useLog){
        htimeSubtracted->GetYaxis()->SetRangeUser(0.5,height*1.1);
    }
    else{
        htimeSubtracted->GetYaxis()->SetRangeUser(0,height*1.1);
    }
    htimeSubtracted->SetTitle(title+";#Deltat [ns];Count");
    htimeSubtracted->SetLineColor(kCyan);
    htimeSubtracted->SetFillColor(kCyan);
    htimeSubtracted->Draw();
    TF1 * fdecayFixTau = new TF1("fdecayFixTau","[0]*exp(-x/2195)+[1]*exp(-x/165)",0,tmax);
    htimeSubtracted->Fit(fdecayFixTau,"L","",tminFit,tmaxFit);
    text_decay->Draw("SAME");

    fdecayFixTau->Draw("SAME");

    canv0->SaveAs(Form("fitT.%s.%s.%s.all.png",run.Data(),trigName.Data(),name.Data()));
    canv1->SaveAs(Form("fitT.%s.%s.%s.subtrackSmall.png",run.Data(),trigName.Data(),name.Data()));

    //htime->GetXaxis()->SetRangeUser(-400,-400+1170*1.5);
    //htimeSubtracted->GetXaxis()->SetRangeUser(-400,-400+1170*1.5);
    //text_decay->SetX(-400+1170/2.);text_decay->Draw();
    //text_big->SetX(-400+1170/2.);text_big->Draw();
    //text_small->SetX(-400+1170/2.);text_small->Draw();
    //canv0->Update();
    //canv1->Update();
    //canv0->SaveAs(Form("fitT.%s.%s.%s.zoomin.all.pdf",run.Data(),trigName.Data(),name.Data()));
    //canv0->SaveAs(Form("fitT.%s.%s.%s.zoomin.all.png",run.Data(),trigName.Data(),name.Data()));
    //canv1->SaveAs(Form("fitT.%s.%s.%s.zoomin.subtrackSmall.png",run.Data(),trigName.Data(),name.Data()));
}
