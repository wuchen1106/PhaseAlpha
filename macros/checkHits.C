{
    TString suffix = "";

    TChain * t = new TChain("t");
    TString runname = "run00370-00505_00584-00624_00773-01000"; TString comment = "mu trigger, 12 mm, -90 mV"; t->Add("anaout/run00370-00505.pair.root"); t->Add("anaout/run00584-00624.pair.root"); t->Add("anaout/run00773-01000.pair.root");
    //TString runname = "run00700-00772"; TString comment = "mu trigger, 63 mm, -90 mV"; t->Add(Form("anaout/%s.pair.root",runname.Data()));
    //TString runname = "run01006-01040"; TString comment = "mu trigger, 12 mm, -90 mV, SO off, DS off"; t->Add(Form("anaout/%s.pair.root",runname.Data()));
    //TString runname = "run00506-00530"; TString comment = "4-fold trigger, 12 mm, -40/-90 mV"; t->Add(Form("anaout/%s.pair.root",runname.Data())); tveto_min = -10; tveto_max = 10; // hmin[0] = 50; hmin[1] = 50;
    //TString runname = "run00531-00583"; TString comment = "5-fold trigger, 12 mm, -35 mV"; t->Add(Form("anaout/%s.pair.root",runname.Data())); tveto_min = -10; tveto_max = 10;// hmin[0] = 40; hmin[1] = 40;

    int useLog4Time = 1;
    //double t_min = 40; useLog4Time = 0;
    //double t_max = 400;
    //double t_step = 1.6;
    double t_min = -400.8;
    double t_max = 1200.8;
    double t_step = 1.6;
    int t_n = (t_max-t_min)/t_step;

    bool ignore_tD = true;
    //bool skipABinPeak12 = true; suffix = suffix+"_skipPeak12";
    bool skipABinPeak12 = false;
    bool skipVeto = true; suffix = suffix+"_skipVeto"; comment+=", skip events with T1/T2 hits in veto window";
    int    ADCmax = 13;
    int    ADCmin = -128;
    double ADCstep[4] = {7.8277888298034667969,7.8277888298034667969,3.9138944149017333984,3.9138944149017333984};
    double hmin[4] = {50,50,40,25};

    double tveto_min = -80;
    double tveto_max = 40;

    TString cutsAdd("1");
    TString cutsTrig("&&foundTrigger");
    TString cutsHeight("");
    TString cutsTime("");
    if (skipABinPeak12) cutsAdd = cutsAdd+"&&!peak1&&!peak2";
    if (skipVeto) cutsAdd = cutsAdd+"&&!veto";

    // count events
    double nEvents = t->GetEntries(cutsAdd);
    double nGoodEvents = t->GetEntries(cutsAdd+cutsTrig);

    // for fitting
    TF1 * f2expo = new TF1("f2expo","expo(0)+expo(2)",tveto_max,400); f2expo->SetLineColor(kRed);
    TF1 * f1expo = new TF1("f1expo","expo",tveto_max,400); f1expo->SetLineColor(kBlue);
    TF1 * f1exponeg = new TF1("f1exponeg","expo",-400,tveto_min); f1exponeg->SetLineColor(kCyan);
    TF1 * f2exponeg = new TF1("f2exponeg","expo(0)+expo(2)",-400,tveto_min); f2exponeg->SetLineColor(kMagenta);
    TF1 * fgaus = new TF1("fgaus","gaus",-20,20);

    // for drawing
    TCanvas * canvas = 0;
    TPad * pads[5];
    TLegend * legend = 0;
    TH1D * h1 = 0;
    TH2D * h2 = 0;
    TLine * line = 0;
    TLatex * texttemp = 0;
    TLatex * text = new TLatex(0.03,0.93,"");
    text->SetTextFont(42);
    text->SetTextSize(0.027);
    gStyle->SetOptStat(0);
    //gStyle->SetOptFit();

    // first of all, draw the trigger pair
    canvas = new TCanvas();
    canvas->SetCanvasSize(1024,1152);
    pads[0] = new TPad("padstrig0","",0  ,0.45,0.5,0.9); pads[0]->Draw();
    pads[1] = new TPad("padstrig1","",0.5,0.45,1  ,0.9); pads[1]->Draw();
    pads[2] = new TPad("padstrig2","",0  ,0   ,0.5,0.45); pads[2]->Draw();
    pads[3] = new TPad("padstrig3","",0.5,0   ,1  ,0.45); pads[3]->Draw();
    pads[0]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0), gPad->SetLogz(1);
    int runId_min = t->GetMinimum("run");
    int runId_max = t->GetMaximum("run");
    t->Draw(Form("run:(trigger_tA+trigger_tB)/2>>htrig_rt(75,340.4,400.4,%d,%.1f,%.1f)",runId_max-runId_min+1,runId_min-0.5,runId_max+0.5),cutsAdd+cutsTrig,"COLZ");
    h2 = (TH2D*)gDirectory->Get("htrig_rt");
    h2->SetTitle(Form("Trigger time among the runs;t_{trig} = (t_{A} + t_{B})/2 [ns];Run ID"));
    //
    pads[1]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0);
    t->Draw("trigger_tA-trigger_tB>>htrig_dt(51,-20.4,20.4)",cutsAdd+cutsTrig,"HIST");
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
    t->Draw(Form("trigger_hB:trigger_hA>>htrig_hAB(%d,%.16f,%.16f,%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep[0],(ADCmax-ADCmin+0.5)*ADCstep[0],ADCmax-ADCmin+1,-0.5*ADCstep[0],(ADCmax-ADCmin+0.5)*ADCstep[0]),cutsAdd+cutsTrig,"COLZ");
    h2 = (TH2D*)gDirectory->Get("htrig_hAB");
    h2->SetTitle("Height in T0L and T0R in the trigger;Peak Height in chA [mV];Peak Height in chB [mV]");
    //
    legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->SetFillStyle(0);
    pads[3]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0);
    t->Draw(Form("trigger_hA>>htrig_hA(%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep[0],(ADCmax-ADCmin+0.5)*ADCstep[0]),cutsAdd+cutsTrig,"HIST");
    h1 = (TH1D*)gDirectory->Get("htrig_hA");
    h1->SetLineColor(kRed);
    h1->SetTitle("Height in T0L and T0R in the trigger;Peak Height [mV]");
    legend->AddEntry(h1,"T0L (A)");
    t->Draw(Form("trigger_hB>>htrig_hB(%d,%.16f,%.16f)",ADCmax-ADCmin+1,-0.5*ADCstep[0],(ADCmax-ADCmin+0.5)*ADCstep[0]),cutsAdd+cutsTrig,"HISTSAME");
    h1 = (TH1D*)gDirectory->Get("htrig_hB");
    h1->SetLineColor(kBlue);
    legend->AddEntry(h1,"T0R (B)");
    legend->Draw();
    canvas->cd();
    text->SetText(0.03,0.95,Form("#splitline{%s, %s}{ %.0f events, found trigger in %.0f (%.2f%%) events}",runname.Data(),comment.Data(),nEvents,nGoodEvents,nGoodEvents/nEvents*100));
    text->Draw();
    canvas->SaveAs(Form("results/trigger.%s%s.png",runname.Data(),suffix.Data()));

    // second, draw the pairs
    canvas = new TCanvas();
    canvas->SetCanvasSize(1024,1152);
    pads[0] = new TPad("padspair0","",0  ,0.6 ,0.5,0.9); pads[0]->Draw();
    pads[1] = new TPad("padspair1","",0.5,0.6 ,1  ,0.9); pads[1]->Draw();
    pads[2] = new TPad("padspair2","",0  ,0.3 ,0.5,0.6); pads[2]->Draw();
    pads[3] = new TPad("padspair3","",0.5,0.3 ,1  ,0.6); pads[3]->Draw();
    pads[4] = new TPad("padsPair4","",0  ,0   ,1  ,0.3); pads[4]->Draw();
    // in T0 (excluding the trigger pair)
    for (int iPair = 0; iPair<2; iPair++){
        int ch1 = iPair*2;
        int ch2 = ch1+1;
        cutsHeight = Form("&&ph%c>%.1f&&ph%c>%.1f",'A'+ch1,hmin[ch1],'A'+ch2,hmin[ch2]);
        cutsTime = Form("&&(pt%c+pt%c)/2-(trigger_tA+trigger_tB)/2>%.1f&&(pt%c+pt%c)/2-(trigger_tA+trigger_tB)/2<=%.1f",'A'+ch1,'A'+ch2,t_min,'A'+ch1,'A'+ch2,t_max);
        //
        pads[0]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0); gPad->SetLogz(1);
        t->Draw(Form("ph%c:ph%c>>hp_hh(%d,%.16f,%.16f,%d,%.16f,%.16f)",'A'+ch2,'A'+ch2,ADCmax-ADCmin+1,-0.5*ADCstep[0],(ADCmax-ADCmin+0.5)*ADCstep[0],ADCmax-ADCmin+1,-0.5*ADCstep[1],(ADCmax-ADCmin+0.5)*ADCstep[1]),cutsAdd+cutsTrig+cutsTime,"COLZ");
        h2 = (TH2D*)gDirectory->Get("hp_hh");
        h2->SetTitle(Form(";Height in ch%c [mV];Height in ch%c [mV]",'A'+ch1,'A'+ch2));
        line = new TLine(0,hmin[ch2],(ADCmax-ADCmin+0.5)*ADCstep[0],hmin[ch2]); line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
        line = new TLine(hmin[ch1],0,hmin[ch1],(ADCmax-ADCmin+0.5)*ADCstep[1]); line->SetLineColor(kBlue); line->SetLineStyle(2); line->Draw();
        //
        pads[1]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0);
        t->Draw(Form("pt%c-pt%c>>hp_dt(51,-20.4,20.4)",'A'+ch1,'A'+ch2),cutsAdd+cutsTrig+cutsHeight+cutsTime,"HIST");
        h1 = (TH1D*)gDirectory->Get("hp_dt");
        h1->SetTitle(Form(";#Deltat = t_{%c} - t_{%c} [ns];Count",'A'+ch1,'A'+ch2));
        h1->Fit(fgaus,"","",-20,20);
        double nPairsGood = h1->Integral();
        fgaus->Draw("SAME");
        texttemp = new TLatex(5,h1->GetMaximum()*0.8,Form("Mean %.1f ns, #sigma %.1f ns",fgaus->GetParameter(1),fgaus->GetParameter(2)));
        texttemp->SetTextFont(42);
        texttemp->SetTextSize(0.04);
        texttemp->Draw();
        //
        pads[2]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0); gPad->SetLogz(1);
        t->Draw(Form("(ph%c+ph%c)/2:(pt%c+pt%c)/2-(trigger_tA+trigger_tB)/2>>hp_ht(%d,%.1f,%.1f,%d,%.16f,%.16f)",'A'+ch1,'A'+ch2,'A'+ch1,'A'+ch2,t_n,t_min,t_max,ADCmax-ADCmin+1,-0.5*ADCstep[0],(ADCmax-ADCmin+0.5)*ADCstep[0]),cutsAdd+cutsTrig,"COLZ");
        h2 = (TH2D*)gDirectory->Get("hp_ht");
        h2->SetTitle(Form(";(t_{%c} + t_{%c})/2 - (trig_{A} + trig_{B})/2 [ns];Height (h_{%c} + h_{%c})/2 [mV]",'A'+ch1,'A'+ch2,'A'+ch1,'A'+ch2));
        line = new TLine(t_min,hmin[ch1],t_max,hmin[ch1]); line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
        line = new TLine(t_min,hmin[ch2],t_max,hmin[ch2]); line->SetLineColor(kBlue); line->SetLineStyle(2); line->Draw();
        //
        legend = new TLegend(0.7,0.7,0.9,0.9);
        legend->SetFillStyle(1);
        pads[3]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(1);
        t->Draw(Form("ph%c>>hp_h%c(%d,%.16f,%.16f)",'A'+ch1,'A'+ch1,ADCmax-ADCmin+1,-0.5*ADCstep[0],(ADCmax-ADCmin+0.5)*ADCstep[0]),cutsAdd+cutsTrig+cutsTime,"HIST");
        h1 = (TH1D*)gDirectory->Get(Form("hp_h%c",'A'+ch1));
        h1->SetLineColor(kRed);
        h1->SetTitle("Height [mV]");
        line = new TLine(hmin[ch1],0.5,hmin[ch1],h1->GetMaximum()); line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
        legend->AddEntry(h1,Form("ch%c",'A'+ch1));
        t->Draw(Form("ph%c>>hp_h%c(%d,%.16f,%.16f)",'A'+ch1,'A'+ch1,ADCmax-ADCmin+1,-0.5*ADCstep[1],(ADCmax-ADCmin+0.5)*ADCstep[1]),cutsAdd+cutsTrig+cutsTime,"HISTSAME");
        h1 = (TH1D*)gDirectory->Get(Form("hp_h%c",'A'+ch2));
        h1->SetLineColor(kBlue);
        line = new TLine(hmin[ch2],0.5,hmin[ch2],h1->GetMaximum()); line->SetLineColor(kBlue); line->SetLineStyle(2); line->Draw();
        double nPairs = h1->Integral();
        legend->AddEntry(h1,Form("ch%c",'A'+ch2));
        legend->Draw();
        //
        pads[4]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(1);
        t->Draw(Form("(pt%c+pt%c)/2-(trigger_tA+trigger_tB)/2>>hp_t(%d,%.1f,%.1f)",'A'+ch1,'A'+ch2,t_n,t_min,t_max),cutsAdd+cutsTrig+cutsHeight,"HIST");
        h1 = (TH1D*)gDirectory->Get("hp_t");
        h1->SetTitle(Form(";t_{pair} = (t_{%c} + t_{%c})/2 [ns];Count",'A'+ch1,'A'+ch2));
        //
        canvas->cd();
        text->SetText(0.03,0.95,Form("#splitline{%s, %s}{#splitline{ %.0f events, found trigger in %.0f (%.2f%%) events}{In T0, %.0f (%.1e/evt) pairs within %.0f~%.0f ns, %.0f (%.1e/evt) pairs over %.0f mV (T0L) %.0f mV (T0R)}}",runname.Data(),comment.Data(),nEvents,nGoodEvents,nGoodEvents/nEvents*100,nPairs,nPairs/nGoodEvents,t_min,t_max,nPairsGood,nPairsGood/nGoodEvents,hmin[ch1],hmin[ch2]));
        text->Draw();
        canvas->SaveAs(Form("results/pair%c%c.t%.0fns_%.0fns_%.1fns.%s%s.hmin%c%.0fmV.hmin%c%.0fmV.png",'A'+ch1,'A'+ch2,t_min,t_max,t_step,runname.Data(),suffix.Data(),'A'+ch1,hmin[ch1],'A'+ch2,hmin[ch2]));
    }

    // at last, draw the left over hits
    canvas = new TCanvas();
    canvas->SetCanvasSize(1024,768);
    pads[0] = new TPad("padsHit0","",0,0.45,0.5,0.9); pads[0]->Draw();
    pads[1] = new TPad("padsHit1","",0.5,0.45,1,0.9); pads[1]->Draw();
    pads[2] = new TPad("padsHit2","",0,0,1,0.45); pads[2]->Draw();
    //
    for (int ch = 0; ch<4; ch++){
        cutsHeight = Form("&&height%c>%.1f",'A'+ch,hmin[ch]);
        cutsTime = Form("&&start%c-(trigger_tA+trigger_tB)/2>=%.1f&&start%c-(trigger_tA+trigger_tB)/2<=%.1f",'A'+ch,t_min,'A'+ch,t_max);
        pads[0]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(0); gPad->SetLogz(1);
        t->Draw(Form("height%c:start%c-(trigger_tA+trigger_tB)/2>>hh_ht%c(%d,%.1f,%.1f,%d,%.16f,%.16f)",'A'+ch,'A'+ch,'A'+ch,t_n,t_min,t_max,ADCmax-ADCmin+1,-0.5*ADCstep[ch],(ADCmax-ADCmin+0.5)*ADCstep[ch]),cutsAdd+cutsTrig,"COLZ");
        h2 = (TH2D*)gDirectory->Get(Form("hh_ht%c",'A'+ch));
        h2->SetTitle("Height VS time;t - (trig_{A} + trig_{B})/2 [ns];Height [mV]");
        line = new TLine(t_min,hmin[ch],t_max,hmin[ch]); line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
        //
        pads[1]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(1);
        t->Draw(Form("height%c>>hh_h%c(%d,%.16f,%.16f)",'A'+ch,'A'+ch,ADCmax-ADCmin+1,-0.5*ADCstep[ch],(ADCmax-ADCmin+0.5)*ADCstep[ch]),cutsAdd+cutsTrig+cutsTime,"HIST");
        h1 = (TH1D*)gDirectory->Get(Form("hh_h%c",'A'+ch));
        h1->SetTitle("Peak Height;Peak Height [mV]");
        line = new TLine(hmin[ch],0.5,hmin[ch],h1->GetMaximum()); line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
        double nHits = h1->Integral();
        //
        pads[2]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetLogy(useLog4Time);
        t->Draw(Form("start%c-(trigger_tA+trigger_tB)/2>>hh_t%c(%d,%.1f,%.1f)",'A'+ch,'A'+ch,t_n,t_min,t_max),cutsAdd+cutsTrig+cutsHeight+cutsTime,"HIST");
        h1 = (TH1D*)gDirectory->Get(Form("hh_t%c",'A'+ch));
        h1->SetTitle(Form("Hit time;t [ns];Count"));
        double nHitsGood = h1->Integral();
        //
        canvas->cd();
        text->SetText(0.03,0.95,Form("#splitline{%s, %s}{#splitline{ %.0f events, found trigger in %.0f (%.2f%%) events}{In ch%c, %.0f (%.1e/evt) single hits within %.0f~%.0f ns, %.0f (%.1e/evt) over %.0f mV}}",runname.Data(),comment.Data(),nEvents,nGoodEvents,nGoodEvents/nEvents*100,'A'+ch,nHits,nHits/nGoodEvents,t_min,t_max,nHitsGood,nHitsGood/nGoodEvents,hmin[ch]));
        text->Draw();
        canvas->SaveAs(Form("results/Hit%c.t%.0fns_%.0fns_%.1fns.%s%s.hmin%.0fmV.png",'A'+ch,t_min,t_max,t_step,runname.Data(),suffix.Data(),hmin[ch]));
    }

}
