{
    TString suffix = "";

    int useLog4LastFig = 1;
    //double tzoom_min = 40; useLog4LastFig = 0;
    //double tzoom_max = 400;
    //double tzoom_step = 1.6;
    //double tzoom_min = -400.8;
    //double tzoom_max = 1200.8;
    //double tzoom_step = 1.6;
    double tzoom_min = 1200;
    double tzoom_max = 1200.8;
    double tzoom_step = 80;
    int tzoom_n = (tzoom_max-tzoom_min)/tzoom_step;

    bool ignore_tD = true;
    bool skipABinPeak12 = true; suffix = suffix+"_skipPeak12";
    //bool skipABinPeak12 = false;
    bool skipVeto = true; suffix = suffix+"_skipVeto";
    double   ADCstep = 7.8277888298035;
    double   ADCstep12 = 3.9138944149017;
    int ADCmax = 13;
    int ADCmin = -128;

    double tAB_left = -5;
    double tAB_right = 2;

    double tCD_left = -11;
    double tCD_right = -3;
    //double aC_min = 14;
    //double aD_min = 8;
    //double hminA = 20;
    //double hminB = 20;
    double hminA = 50;
    double hminB = 50;
    //double hminA = 100;
    //double hminB = 100;
    double hminC = 40;
    double hminD = 25;

    double tveto_min = -80;
    double tveto_max = 40;

    TChain * t = new TChain("t");
    TString runname = "run00370-00505_00773-01005"; TString comment = "mu, 12 mm, -90 mV"; t->Add("anaout/run00370-00505.pair.root"); t->Add("anaout/run00773-01005.pair.root");
    //TString runname = "run00700-00772"; TString comment = "mu, 63 mm, -90 mV"; t->Add(Form("anaout/%s.pair.root",runname.Data()));
    //TString runname = "run01006-01040"; TString comment = "mu, 12 mm, -90 mV, SO off, DS off"; t->Add(Form("anaout/%s.pair.root",runname.Data()));
    //TString runname = "run00506-00530"; TString comment = "4-fold, 12 mm, -40/-90 mV"; t->Add(Form("anaout/%s.pair.root",runname.Data())); tveto_min = -10; tveto_max = 10; // hminA = 50; hminB = 50;
    //TString runname = "run00531-00583"; TString comment = "5-fold, 12 mm, -35 mV"; t->Add(Form("anaout/%s.pair.root",runname.Data())); tveto_min = -10; tveto_max = 10;// hminA = 40; hminB = 40;

    TString cuts0 = "foundTrigger";
    TString cuts;
    TString cutsAdd = "";
    if (skipABinPeak12) cutsAdd = cutsAdd+"&&!peak1&&!peak2";
    if (skipVeto) cutsAdd = cutsAdd+"&&!veto";

    TF1 * f2expo = new TF1("f2expo","expo(0)+expo(2)",tveto_max,400); f2expo->SetLineColor(kRed);
    TF1 * f1expo = new TF1("f1expo","expo",tveto_max,400); f1expo->SetLineColor(kBlue);
    TF1 * f1exponeg = new TF1("f1exponeg","expo",-400,tveto_min); f1exponeg->SetLineColor(kCyan);
    TF1 * f2exponeg = new TF1("f2exponeg","expo(0)+expo(2)",-400,tveto_min); f2exponeg->SetLineColor(kMagenta);
    TF1 * fgaus = new TF1("fgaus","gaus",-20,20);

    TCanvas * canvas = 0;
    TH1D * h1 = 0;
    TH2D * h2 = 0;
    TLine * line = 0;
    TLegend * legend = 0;
    TLatex * text = 0;
    TLatex * texttemp = 0;
    gStyle->SetOptStat(0);
    //gStyle->SetOptFit();

    canvas = new TCanvas();
    canvas->SetCanvasSize(1024,768);
    TPad * pads[4];
    pads[0] = new TPad("pads0","",0,0.45,0.5,0.9); pads[0]->Draw();
    pads[1] = new TPad("pads1","",0.5,0.45,1,0.9); pads[1]->Draw();
    pads[2] = new TPad("pads2","",0,0,0.5,0.45); pads[2]->Draw();
    pads[3] = new TPad("pads2","",0.5,0,1,0.45); pads[3]->Draw();
    text = new TLatex(0.03,0.93,"");
    text->SetTextFont(42);
    text->SetTextSize(0.027);
    text->Draw();

    // first of all, draw the trigger pair
    //
    pads[0]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(1);
    if (runname=="run00370-00505_00773-01005"){
        legend = new TLegend(0.1,0.7,0.4,0.9);
        cuts = cuts0+"&&run<500"+cutsAdd;
        t->Draw("(trigger_tA+trigger_tB)/2>>htrig_t0(75,340.4,400.4)",cuts,"HIST");
        h1 = (TH1D*)gDirectory->Get("htrig_t0");
        h1->GetYaxis()->SetRangeUser(0.5,h1->GetMaximum()*10);
        h1->SetLineColor(kRed);
        h1->SetTitle(Form("Trigger time;t_{trig} = (t_{A} + t_{B})/2 [ns];Count"));
        legend->AddEntry(h1,"run00370-00499");
        cuts = cuts0+"&&run>=500&&run<=505"+cutsAdd;
        t->Draw("(trigger_tA+trigger_tB)/2>>htrig_t1(75,340.4,400.4)",cuts,"HISTSAME");
        h1 = (TH1D*)gDirectory->Get("htrig_t1");
        h1->SetLineColor(kGreen);
        legend->AddEntry(h1,"run00500-00505");
        cuts = cuts0+"&&run>505"+cutsAdd;
        t->Draw("(trigger_tA+trigger_tB)/2>>htrig_t2(75,340.4,400.4)",cuts,"HISTSAME");
        h1 = (TH1D*)gDirectory->Get("htrig_t2");
        h1->SetLineColor(kBlue);
        legend->AddEntry(h1,"run00773-01005");
        legend->Draw();
    }
    else{
        t->Draw("(trigger_tA+trigger_tB)/2>>htrig_t(75,340.4,400.4)",cuts0+cutsAdd,"HIST");
        h1 = (TH1D*)gDirectory->Get("htrig_t");
        h1->SetTitle(Form("Trigger time;t_{trig} = (t_{A} + t_{B})/2 [ns];Count"));
    }
    //
    pads[1]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0);
    t->Draw("trigger_tA-trigger_tB>>htrig_dt(51,-20.4,20.4)",cuts0+cutsAdd,"HIST");
    h1 = (TH1D*)gDirectory->Get("htrig_dt");
    h1->SetTitle(Form("Time difference from T0L and T0R in one trigger;#Deltat = t_{A} - t_{B} [ns];Count"));
    h1->Fit(fgaus,"","",-20,20);
    fgaus->Draw("SAME");
    texttemp = new TLatex(5,h1->GetMaximum()*0.8,Form("Mean %.1f ns, #sigma %.1f ns",fgaus->GetParameter(1),fgaus->GetParameter(2)));
    texttemp->SetTextFont(42);
    texttemp->SetTextSize(0.04);
    texttemp->Draw();
    //
    pads[2]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0); gPad->SetLogz(1);
    t->Draw(Form("trigger_hB:trigger_hA>>htrig_hAB(%d,%.16f,%.16f,%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep,(ADCmax-ADCmin+0.5)*ADCstep,ADCmax-ADCmin+1,-0.5*ADCstep,(ADCmax-ADCmin+0.5)*ADCstep),cuts0+cutsAdd,"COLZ");
    h2 = (TH2D*)gDirectory->Get("htrig_hAB");
    h2->SetTitle("Height in T0L and T0R in the trigger;Peak Height in chA [mV];Peak Height in chB [mV]");
    //
    legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->SetFillStyle(0);
    pads[3]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0);
    t->Draw(Form("trigger_hA>>htrig_hA(%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep,(ADCmax-ADCmin+0.5)*ADCstep),cuts0+cutsAdd,"HIST");
    h1 = (TH1D*)gDirectory->Get("htrig_hA");
    h1->SetLineColor(kRed);
    h1->SetTitle("Height in T0L and T0R in the trigger;Peak Height [mV]");
    legend->AddEntry(h1,"T0L (A)");
    t->Draw(Form("trigger_hB>>htrig_hB(%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep,(ADCmax-ADCmin+0.5)*ADCstep),cuts0+cutsAdd,"HISTSAME");
    h1 = (TH1D*)gDirectory->Get("htrig_hB");
    h1->SetLineColor(kBlue);
    legend->AddEntry(h1,"T0R (B)");
    legend->Draw();
    double nGoodEvents = h1->GetEntries();
    text->SetText(0.03,0.93,Form("%s, %s: %lld events, %.0f (%.2f%%) with trigger pair",runname.Data(),comment.Data(),t->GetEntries(),nGoodEvents,nGoodEvents*100./t->GetEntries()));
    canvas->SaveAs(Form("results/trigger.%s%s.png",runname.Data(),suffix.Data()));

    // second, draw the AB pairs (excluding the trigger)
    //
    pads[0]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(1);
    t->Draw("(ptA+ptB)/2-(trigger_tA+trigger_tB)/2>>hp_tAB(500,-500,9500)",cuts0+cutsAdd+Form("&&phA>%.1f&&phB>%.1f",hminA,hminB),"HIST");
    h1 = (TH1D*)gDirectory->Get("hp_tAB");
    h1->SetTitle(Form("Pair time;t_{pair} = (t_{A} + t_{B})/2 [ns];Count"));
    //
    pads[1]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0);
    t->Draw("ptA-ptB>>hp_dtAB(51,-20.4,20.4)",cuts0+cutsAdd+Form("&&phA>%.1f&&phB>%.1f",hminA,hminB),"HIST");
    h1 = (TH1D*)gDirectory->Get("hp_dtAB");
    h1->SetTitle(Form("Time difference from T0L and T0R in one pair;#Deltat = t_{A} - t_{B} [ns];Count"));
    h1->Fit(fgaus,"","",-20,20);
    double nPairsABGood = h1->Integral();
    fgaus->Draw("SAME");
    texttemp = new TLatex(5,h1->GetMaximum()*0.8,Form("Mean %.1f ns, #sigma %.1f ns",fgaus->GetParameter(1),fgaus->GetParameter(2)));
    texttemp->SetTextFont(42);
    texttemp->SetTextSize(0.04);
    texttemp->Draw();
    //
    pads[2]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0); gPad->SetLogz(1);
    t->Draw(Form("(phB+phA)/2:(ptB+ptA)/2-(trigger_tA+trigger_tB)/2>>hp_h2AB(1000,-500,9500,%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep,(ADCmax-ADCmin+0.5)*ADCstep),cuts0+cutsAdd,"COLZ");
    h2 = (TH2D*)gDirectory->Get("hp_h2AB");
    h2->SetTitle("Height VS time of T0 pairs;(t_{A} + t_{B})/2 - (trig_{A} + trig_{B})/2 [ns];Height (h_{A} + h_{B})/2 [mV]");
    line = new TLine(tzoom_min,hminA,tzoom_max,hminA); line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
    line = new TLine(tzoom_min,hminB,tzoom_max,hminB); line->SetLineColor(kBlue); line->SetLineStyle(2); line->Draw();
    //
    legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->SetFillStyle(1);
    pads[3]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(1);
    t->Draw(Form("phA>>hp_hA(%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep,(ADCmax-ADCmin+0.5)*ADCstep),cuts0+cutsAdd,"HIST");
    h1 = (TH1D*)gDirectory->Get("hp_hA");
    h1->SetLineColor(kRed);
    h1->SetTitle("Height in T0L and T0R in the pair;Peak Height [mV]");
    line = new TLine(hminA,0.5,hminA,h1->GetMaximum()); line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
    legend->AddEntry(h1,"T0L (A)");
    t->Draw(Form("phB>>hp_hB(%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep,(ADCmax-ADCmin+0.5)*ADCstep),cuts0+cutsAdd,"HISTSAME");
    h1 = (TH1D*)gDirectory->Get("hp_hB");
    h1->SetLineColor(kBlue);
    line = new TLine(hminB,0.5,hminB,h1->GetMaximum()); line->SetLineColor(kBlue); line->SetLineStyle(2); line->Draw();
    double nPairsAB = h1->Integral();
    legend->AddEntry(h1,"T0R (B)");
    legend->Draw();
    canvas->cd();
    text->SetText(0.03,0.93,Form("%s, %s: %.0f (%.0f) pairs in T0L&R from %.0f good events",runname.Data(),comment.Data(),nPairsAB,nPairsABGood,nGoodEvents));
    text->Draw();
    canvas->SaveAs(Form("results/pairAB.%s%s.hminA%.0fmV.hminB%.0fmV.png",runname.Data(),suffix.Data(),hminA,hminB));

    // second-2, draw the AB pairs (excluding the trigger) in zoom-in view
    //
    pads[0]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(1);
    t->Draw(Form("(ptA+ptB)/2-(trigger_tA+trigger_tB)/2>>hp_tAB(%d,%.1f,%.1f)",tzoom_n,tzoom_min,tzoom_max),cuts0+cutsAdd+Form("&&phA>%.1f&&phB>%.1f",hminA,hminB),"HIST");
    h1 = (TH1D*)gDirectory->Get("hp_tAB");
    h1->SetTitle(Form("Pair time;t_{pair} = (t_{A} + t_{B})/2 [ns];Count"));
    //
    pads[1]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0);
    t->Draw("ptA-ptB>>hp_dtAB(51,-20.4,20.4)",cuts0+cutsAdd+Form("&&(ptA+ptB)/2-(trigger_tA+trigger_tB)/2>%.1f&&(ptA+ptB)/2-(trigger_tA+trigger_tB)/2<=%.1f&&phA>%.1f&&phB>%.1f",tzoom_min,tzoom_max,hminA,hminB),"HIST");
    h1 = (TH1D*)gDirectory->Get("hp_dtAB");
    h1->SetTitle(Form("Time difference from T0L and T0R in one pair;#Deltat = t_{A} - t_{B} [ns];Count"));
    h1->Fit(fgaus,"","",-20,20);
    nPairsABGood = h1->Integral();
    fgaus->Draw("SAME");
    texttemp = new TLatex(5,h1->GetMaximum()*0.8,Form("Mean %.1f ns, #sigma %.1f ns",fgaus->GetParameter(1),fgaus->GetParameter(2)));
    texttemp->SetTextFont(42);
    texttemp->SetTextSize(0.04);
    texttemp->Draw();
    //
    pads[2]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0); gPad->SetLogz(1);
    t->Draw(Form("(phB+phA)/2:(ptB+ptA)/2-(trigger_tA+trigger_tB)/2>>hp_h2AB(%d,%.1f,%.1f,%d,%.16f,%.16f)",tzoom_n,tzoom_min,tzoom_max,ADCmax-ADCmin+1,-0.5*ADCstep,(ADCmax-ADCmin+0.5)*ADCstep),cuts0+cutsAdd,"COLZ");
    h2 = (TH2D*)gDirectory->Get("hp_h2AB");
    h2->SetTitle("Height VS time of T0 pairs;(t_{A} + t_{B})/2 - (trig_{A} + trig_{B})/2 [ns];Height (h_{A} + h_{B})/2 [mV]");
    line = new TLine(tzoom_min,hminA,tzoom_max,hminA); line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
    line = new TLine(tzoom_min,hminB,tzoom_max,hminB); line->SetLineColor(kBlue); line->SetLineStyle(2); line->Draw();
    //
    legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->SetFillStyle(1);
    pads[3]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(1);
    t->Draw(Form("phA>>hp_hA(%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep,(ADCmax-ADCmin+0.5)*ADCstep),cuts0+cutsAdd+Form("&&(ptA+ptB)/2-(trigger_tA+trigger_tB)/2>%.1f&&(ptA+ptB)/2-(trigger_tA+trigger_tB)/2<=%.1f",tzoom_min,tzoom_max),"HIST");
    h1 = (TH1D*)gDirectory->Get("hp_hA");
    h1->SetLineColor(kRed);
    h1->SetTitle("Height in T0L and T0R in the pair;Peak Height [mV]");
    line = new TLine(hminA,0.5,hminA,h1->GetMaximum()); line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
    legend->AddEntry(h1,"T0L (A)");
    t->Draw(Form("phB>>hp_hB(%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep,(ADCmax-ADCmin+0.5)*ADCstep),cuts0+cutsAdd+Form("&&(ptA+ptB)/2-(trigger_tA+trigger_tB)/2>%.1f&&(ptA+ptB)/2-(trigger_tA+trigger_tB)/2<=%.1f",tzoom_min,tzoom_max),"HISTSAME");
    h1 = (TH1D*)gDirectory->Get("hp_hB");
    h1->SetLineColor(kBlue);
    line = new TLine(hminB,0.5,hminB,h1->GetMaximum()); line->SetLineColor(kBlue); line->SetLineStyle(2); line->Draw();
    nPairsAB = h1->Integral();
    legend->AddEntry(h1,"T0R (B)");
    legend->Draw();
    canvas->cd();
    text->SetText(0.03,0.93,Form("%s, %s: %.0f (%.0f) pairs in T0L&R from %.0f good events",runname.Data(),comment.Data(),nPairsAB,nPairsABGood,nGoodEvents));
    text->Draw();
    canvas->SaveAs(Form("results/pairABzoom%.0fns_%.0fns.%s%s.hminA%.0fmV.hminB%.0fmV.png",tzoom_min,tzoom_max,runname.Data(),suffix.Data(),hminA,hminB));

    // third, draw the CD pairs
    //
    pads[0]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(1);
    if (ignore_tD){
        t->Draw("ptC-(trigger_tA+trigger_tB)/2>>hp_tCD(500,-500,9500)",cuts0+cutsAdd+Form("&&phC>%.1f&&phD>%.1f",hminC,hminD),"HIST");
        h1 = (TH1D*)gDirectory->Get("hp_tCD");
        h1->SetTitle(Form("Pair time (h_{C}>%.0f mV h_{D}>%.0f mV);t_{pair} = t_{C} [ns];Count",hminC,hminD));
    }
    else{
        t->Draw("(ptC+ptD)/2-(trigger_tA+trigger_tB)/2>>hp_tCD(500,-500,9500)",cuts0+cutsAdd+Form("&&phC>%.1f&&phD>%.1f",hminC,hminD),"HIST");
        h1 = (TH1D*)gDirectory->Get("hp_tCD");
        h1->SetTitle(Form("Pair time (h_{C}>%.0f mV h_{D}>%.0f mV);t_{pair} = (t_{C} + t_{D})/2 [ns];Count",hminC,hminD));
    }
    //
    pads[1]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0);
    t->Draw("ptC-ptD>>hp_dtCD(51,-20.4,20.4)",cuts0+cutsAdd+Form("&&phC>%.1f&&phD>%.1f",hminC,hminD),"HIST");
    h1 = (TH1D*)gDirectory->Get("hp_dtCD");
    h1->SetTitle(Form("Time difference from T1/2 in one pair (h_{C}>%.0f mV h_{D}>%.0f mV);#Deltat = t_{C} - t_{D} [ns];Count",hminC,hminD));
    h1->Fit(fgaus,"","",-20,20);
    double nPairsCDGood = h1->Integral();
    fgaus->Draw("SAME");
    texttemp = new TLatex(5,h1->GetMaximum()*0.8,Form("Mean %.1f ns, #sigma %.1f ns",fgaus->GetParameter(1),fgaus->GetParameter(2)));
    texttemp->SetTextFont(42);
    texttemp->SetTextSize(0.04);
    texttemp->Draw();
    //
    pads[2]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0); gPad->SetLogz(1);
    if (ignore_tD){
        t->Draw(Form("phC:ptC-(trigger_tA+trigger_tB)/2>>hp_h2CD(1000,-500,9500,%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep12,(ADCmax-ADCmin+0.5)*ADCstep12),cuts0+cutsAdd,"COLZ");
        h2 = (TH2D*)gDirectory->Get("hp_h2CD");
        h2->SetTitle("Height VS time of T1/2 pairs;t_{C} - (trig_{A} + trig_{B})/2 [ns];Height h_{C} [mV]");
    }
    else{
        t->Draw(Form("(phD+phC)/2:(ptD+ptC)/2-(trigger_tA+trigger_tB)/2>>hp_h2CD(1000,-500,9500,%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep12,(ADCmax-ADCmin+0.5)*ADCstep12),cuts0+cutsAdd,"COLZ");
        h2 = (TH2D*)gDirectory->Get("hp_h2CD");
        h2->SetTitle("Height VS time of T1/2 pairs;(t_{C} + t_{D})/2 - (trig_{A} + trig_{B})/2 [ns];Height (h_{C} + h_{D})/2 [mV]");
    }
    line = new TLine(tzoom_min,hminC,tzoom_max,hminC); line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
    //line = new TLine(tzoom_min,hminD,tzoom_max,hminD); line->SetLineColor(kBlue); line->SetLineStyle(2); line->Draw();
    //
    legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->SetFillStyle(1);
    pads[3]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(1);
    t->Draw(Form("phC>>hp_hC(%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep12,(ADCmax-ADCmin+0.5)*ADCstep12),cuts0+cutsAdd,"HIST");
    h1 = (TH1D*)gDirectory->Get("hp_hC");
    h1->SetLineColor(kRed);
    h1->SetTitle("Height in T1 and T2 in the pair;Peak Height [mV]");
    line = new TLine(hminC,0.5,hminC,h1->GetMaximum()); line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
    legend->AddEntry(h1,"T1 (C)");
    t->Draw(Form("phD>>hp_hD(%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep12,(ADCmax-ADCmin+0.5)*ADCstep12),cuts0+cutsAdd,"HISTSAME");
    h1 = (TH1D*)gDirectory->Get("hp_hD");
    h1->SetLineColor(kBlue);
    line = new TLine(hminD,0.5,hminD,h1->GetMaximum()); line->SetLineColor(kBlue); line->SetLineStyle(2); line->Draw();
    double nPairsCD = h1->Integral();
    legend->AddEntry(h1,"T2 (D)");
    legend->Draw();
    canvas->cd();
    text->SetText(0.03,0.93,Form("%s, %s: %.0f (%.0f) pairs in T1/2 from %.0f good events",runname.Data(),comment.Data(),nPairsCD,nPairsCDGood,nGoodEvents));
    text->Draw();
    canvas->SaveAs(Form("results/pairCD.%s%s.hminC%.0fmV.hminD%.0fmV.png",runname.Data(),suffix.Data(),hminC,hminD));

    // third-2, draw the CD pairs in zoom view
    //
    pads[0]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(useLog4LastFig);
    if (ignore_tD){
        t->Draw(Form("ptC-(trigger_tA+trigger_tB)/2>>hp_tCD(%d,%.1f,%.1f)",tzoom_n,tzoom_min,tzoom_max),cuts0+cutsAdd+Form("&&phC>%.1f&&phD>%.1f",hminC,hminD),"HIST");
        h1 = (TH1D*)gDirectory->Get("hp_tCD");
        h1->SetTitle(Form("Pair time (h_{C}>%.0f mV h_{D}>%.0f mV);t_{pair} = t_{C} [ns];Count",hminC,hminD));
    }
    else{
        t->Draw(Form("(ptC+ptD)/2-(trigger_tA+trigger_tB)/2>>hp_tCD(%d,%.1f,%.1f)",tzoom_n,tzoom_min,tzoom_max),cuts0+cutsAdd+Form("&&phC>%.1f&&phD>%.1f",hminC,hminD),"HIST");
        h1 = (TH1D*)gDirectory->Get("hp_tCD");
        h1->SetTitle(Form("Pair time (h_{C}>%.0f mV h_{D}>%.0f mV);t_{pair} = (t_{C} + t_{D})/2 [ns];Count",hminC,hminD));
    }
    h1->Fit(f2expo,"","",tveto_min,150);
    h1->Fit(f2expo,"","",tveto_min,400);
    h1->Fit(f1expo,"","",tveto_min,400);
    f1expo->Draw("SAME");
    f2expo->Draw("SAME");
    if (tzoom_min>-400){
        texttemp = new TLatex(50,h1->GetMaximum(),Form("#splitline{Single-expo: #tau = %.0f ns}{Double-expo: #tau_{1} = %.0f ns, #tau_{2} = %.0f ns}",-1/f1expo->GetParameter(1),-1/f2expo->GetParameter(1),-1/f2expo->GetParameter(3)));
    }
    else{
        //h1->Fit(f2exponeg,"","",-150,-80);
        //h1->Fit(f2exponeg,"","",-400,-80);
        //texttemp = new TLatex(50,h1->GetMaximum()*0.8,Form("#splitline{Single-expo: #tau = %.0f ns}{#splitline{Double-expo: #tau_{1} = %.0f ns, #tau_{2} = %.0f ns}{Double-expo neg. side: #tau_{1} = %.0f ns, #tau_{2} = %.0f ns}}",-1/f1expo->GetParameter(1),-1/f2expo->GetParameter(1),-1/f2expo->GetParameter(3),1/f2exponeg->GetParameter(1),1/f2exponeg->GetParameter(3)));
        //f2exponeg->Draw("SAME");
        h1->Fit(f1exponeg,"","",-150,tveto_min);
        texttemp = new TLatex(50,h1->GetMaximum()*0.8,Form("#splitline{Single-expo: #tau = %.0f ns}{#splitline{Double-expo: #tau_{1} = %.0f ns, #tau_{2} = %.0f ns}{Single-expo neg. side: #tau = %.0f ns}}",-1/f1expo->GetParameter(1),-1/f2expo->GetParameter(1),-1/f2expo->GetParameter(3),1/f1exponeg->GetParameter(1)));
        f1exponeg->Draw("SAME");
    }
    texttemp->SetTextFont(42);
    texttemp->SetTextSize(0.04);
    texttemp->Draw();
    //
    pads[1]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0);
    if (ignore_tD){
        t->Draw("ptC-ptD>>hp_dtCD(51,-20.4,20.4)",cuts0+cutsAdd+Form("&&ptC-(trigger_tA+trigger_tB)/2>%.1f&&ptC-(trigger_tA+trigger_tB)/2<=%.1f&&phC>%.1f&&phD>%.1f",tzoom_min,tzoom_max,hminC,hminD),"HIST");
    }
    else{
        t->Draw("ptC-ptD>>hp_dtCD(51,-20.4,20.4)",cuts0+cutsAdd+Form("&&(ptC+ptD)/2-(trigger_tA+trigger_tB)/2>%.1f&&(ptC+ptD)/2-(trigger_tA+trigger_tB)/2<=%.1f&&phC>%.1f&&phD>%.1f",tzoom_min,tzoom_max,hminC,hminD),"HIST");
    }
    h1 = (TH1D*)gDirectory->Get("hp_dtCD");
    h1->SetTitle(Form("Time difference from T1/2 in one pair (h_{C}>%.0f mV h_{D}>%.0f mV);#Deltat = t_{C} - t_{D} [ns];Count",hminC,hminD));
    h1->Fit(fgaus,"","",-20,20);
    nPairsCDGood = h1->Integral();
    fgaus->Draw("SAME");
    texttemp = new TLatex(5,h1->GetMaximum()*0.8,Form("Mean %.1f ns, #sigma %.1f ns",fgaus->GetParameter(1),fgaus->GetParameter(2)));
    texttemp->SetTextFont(42);
    texttemp->SetTextSize(0.04);
    texttemp->Draw();
    //
    pads[2]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0); gPad->SetLogz(1);
    if (ignore_tD){
        t->Draw(Form("phC:ptC-(trigger_tA+trigger_tB)/2>>hp_h2CD(%d,%.1f,%.1f,%d,%.16f,%.16f)",tzoom_n,tzoom_min,tzoom_max,ADCmax-ADCmin+1,-0.5*ADCstep12,(ADCmax-ADCmin+0.5)*ADCstep12),cuts0+cutsAdd,"COLZ");
        h2 = (TH2D*)gDirectory->Get("hp_h2CD");
        h2->SetTitle("Height VS time of T1/2 pairs;t_{C} - (trig_{A} + trig_{B})/2 [ns];Height h_{C} [mV]");
    }
    else{
        t->Draw(Form("(phD+phC)/2:(ptD+ptC)/2-(trigger_tA+trigger_tB)/2>>hp_h2CD(%d,%.1f,%.1f,%d,%.16f,%.16f)",tzoom_n,tzoom_min,tzoom_max,ADCmax-ADCmin+1,-0.5*ADCstep12,(ADCmax-ADCmin+0.5)*ADCstep12),cuts0+cutsAdd,"COLZ");
        h2 = (TH2D*)gDirectory->Get("hp_h2CD");
        h2->SetTitle("Height VS time of T1/2 pairs;(t_{C} + t_{D})/2 - (trig_{A} + trig_{B})/2 [ns];Height (h_{C} + h_{D})/2 [mV]");
    }
    line = new TLine(tzoom_min,hminC,tzoom_max,hminC); line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
    //line = new TLine(tzoom_min,hminD,tzoom_max,hminD); line->SetLineColor(kBlue); line->SetLineStyle(2); line->Draw();
    //
    legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->SetFillStyle(1);
    pads[3]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(1);
    if (ignore_tD){
        t->Draw(Form("phC>>hp_hC(%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep12,(ADCmax-ADCmin+0.5)*ADCstep12),cuts0+cutsAdd+Form("&&ptC-(trigger_tA+trigger_tB)/2>%.1f&&ptC-(trigger_tA+trigger_tB)/2<=%.1f",tzoom_min,tzoom_max),"HIST");
    }
    else{
        t->Draw(Form("phC>>hp_hC(%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep12,(ADCmax-ADCmin+0.5)*ADCstep12),cuts0+cutsAdd+Form("&&(ptC+ptD)/2-(trigger_tA+trigger_tB)/2>%.1f&&(ptC+ptD)/2-(trigger_tA+trigger_tB)/2<=%.1f",tzoom_min,tzoom_max),"HIST");
    }
    h1 = (TH1D*)gDirectory->Get("hp_hC");
    h1->SetLineColor(kRed);
    h1->SetTitle("Height in T1 and T2 in the pair;Peak Height [mV]");
    line = new TLine(hminC,0.5,hminC,h1->GetMaximum()); line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
    legend->AddEntry(h1,"T1 (C)");
    t->Draw(Form("phD>>hp_hD(%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep12,(ADCmax-ADCmin+0.5)*ADCstep12),cuts0+cutsAdd+Form("&&(ptC+ptD)/2-(trigger_tA+trigger_tB)/2>%.1f&&(ptC+ptD)/2-(trigger_tA+trigger_tB)/2<=%.1f",tzoom_min,tzoom_max),"HISTSAME");
    h1 = (TH1D*)gDirectory->Get("hp_hD");
    h1->SetLineColor(kBlue);
    line = new TLine(hminD,0.5,hminD,h1->GetMaximum()); line->SetLineColor(kBlue); line->SetLineStyle(2); line->Draw();
    nPairsCD = h1->Integral();
    legend->AddEntry(h1,"T2 (D)");
    legend->Draw();
    canvas->cd();
    text->SetText(0.03,0.93,Form("%s, %s: %.0f (%.0f) pairs in T1/2 from %.0f good events",runname.Data(),comment.Data(),nPairsCD,nPairsCDGood,nGoodEvents));
    text->Draw();
    canvas->SaveAs(Form("results/pairCDzoom%.0fns_%.0fns.%s%s.hminC%.0fmV.hminD%.0fmV.png",tzoom_min,tzoom_max,runname.Data(),suffix.Data(),hminC,hminD));
}
