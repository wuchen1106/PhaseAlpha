{
    double ADCstep[4] = {7.8277888298034667969,7.8277888298034667969,3.9138944149017333984,3.9138944149017333984};

    TChain * t = new TChain("t");
    //t->Add("anaout/run00370-00505.event.root");t->Add("anaout/run00584-00624.event.root");t->Add("anaout/run00773-01000.event.root"); TString runname="12mm_mu";
    t->Add("anaout/run00700-00772.event.root");TString runname="63mm_mu";
    //t->Add("anaout/run00509-00517.event.root");t->Add("anaout/run00518-00530.event.root");TString runname="12mm_4fold";
    //t->Add("anaout/run00531-00583.event.root");TString runname="12mm_5fold_p12_40mV25mV";

    TString suffix = ".all";
    TString comment = "all events with 4-fold coincident";
    TString cuts = "foundTrigger";
    TString eventType = "4-fold"; // 4-fold, 2-fold, 3-fold, single
    //cuts+="&&evt_tA!=trig_tA&&evt_hA>0&&evt_hB>0&&evt_hC>0&&evt_hD>0&&(evt_hA+evt_hB)/2>100&&(evt_hC+evt_hD)/2<100";
    cuts+="&&evt_tA!=trig_tA&&evt_hA>0&&evt_hB>0&&evt_hC>0&&evt_hD>0";
    char sizeChar = 'w'; // h, w, a
    int size_nbins[4];
    double size_min[4];
    double size_max[4];
    TString sizeName("");
    TString sizeTitle("");
    if (sizeChar == 'h'){
        sizeName = "height";
        sizeTitle = "height [mV]";
        for (int ch = 0; ch<4; ch++){
            size_nbins[ch] = 131;
            size_min[ch] = -0.5*ADCstep[ch];
            size_max[ch] = (size_nbins[ch]-1+0.5)*ADCstep[ch];
        }
    }
    else if (sizeChar == 'w'){
        sizeName = "width";
        sizeTitle = "width [#samples]";
        for (int ch = 0; ch<4; ch++){
            size_nbins[ch] = 201;
            size_min[ch] = -0.5;
            size_max[ch] = 200.5;
        }
    }
    double t_min = -400;
    double t_max = 9600;
    double t_step = 16;
    int t_n = (t_max-t_min)/t_step;
    int colors[7] = {kRed,kBlue,kMagenta,kCyan,kGreen,kOrange,kGray};

    TCanvas * canv = 0;
    TPad * pads[4];
    double nEvents = t->GetEntries();
    double nGoodEvents = t->GetEntries("foundTrigger");
    TLatex * text = new TLatex(0.03,0.95,"");
    text->SetTextFont(42);
    text->SetTextSize(0.027);
    gStyle->SetOptStat(0);
    TLegend * legend;
    TH1D * h1;
    TH2D * h2;
    double maximum = 0;

    //=======================================================================
    // first, draw the distribution plots
    canv = new TCanvas("dist");
    canv->SetCanvasSize(1500,1000);
    pads[0] = new TPad("pad0","",0,0,1./3,0.9);
    pads[1] = new TPad("pad1","",1./3,0.45,1,0.9);
    pads[2] = new TPad("pad2","",1./3,0,1,0.45);
    for (int i = 0; i<3; i++){
        pads[i]->Draw(); pads[i]->SetGridx(1); pads[i]->SetGridy(1); pads[i]->SetTickx(1); pads[i]->SetTicky(1);
    }
    // size
    pads[0]->cd(); pads[0]->SetLogy(1);
    legend = new TLegend(0.7,0.7,0.9,0.9);
    maximum = 0;
    for (int ch = 0; ch<4; ch++){
        // here skip channels if needed
        // ...
        TString opt = "HISTSAME";
        if (ch == 0) opt = "HIST";
        t->Draw(Form("evt_%c%c>>h_s_%c(%d,%.16e,%.16e",sizeChar,'A'+ch,'A'+ch,size_nbins[ch],size_min[ch],size_max[ch]),cuts,opt);
        h1 = (TH1D*)gDirectory->Get(Form("h_s_%c",'A'+ch));
        h1->SetLineColor(colors[ch]);
        h1->SetTitle(Form(";%s;Count",sizeTitle.Data()));
        if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
        legend->AddEntry(h1,Form("ch%c",'A'+ch));
    }
    h1 = (TH1D*)gDirectory->Get("h_s_A");
    h1->GetYaxis()->SetRangeUser(0.5,maximum*2);
    double nObjects = h1->GetEntries();
    legend->Draw("SAME");
    // size VS time
    pads[1]->cd(); pads[1]->SetLogz(0);
    if (eventType=="4-fold"){
        t->Draw(Form("(evt_%cA+evt_%cB+evt_%cC+evt_%cD)/4:(evt_tA+evt_tB+evt_tC+evt_tD)/4-(trig_tA+trig_tB)/2>>h_st(%d,%.1f,%.1f,%d,%.16e,%.16e)",sizeChar,sizeChar,sizeChar,sizeChar,t_n,t_min,t_max,size_nbins[0],size_min[0],size_max[0]),cuts,"COLZ");
    }
    // time
    pads[2]->cd(); pads[2]->SetLogy(1);
    if (eventType=="4-fold"){
        t->Draw(Form("(evt_tA+evt_tB+evt_tC+evt_tD)/4-(trig_tA+trig_tB)/2>>h_t(%d,%.1f,%.1f)",t_n,t_min,t_max),cuts,"HIST");
    }
    canv->cd();
    text->SetText(0.03,0.95,Form("#splitline{%s, %.0f events, found trigger in %.0f (%.2f%%) events}{%s: %.0f (%.1e/event)}",runname.Data(),nEvents,nGoodEvents,nGoodEvents/nEvents*100,comment.Data(),nObjects,nObjects/nGoodEvents));
    text->Draw();
    canv->SaveAs(Form("dist.%s.%c.%s%s.png",eventType.Data(),sizeChar,runname.Data(),suffix.Data()));

    //=======================================================================
    // second, draw the coincidence plots
    if (eventType!="single"){
        canv = new TCanvas("coin");
        if (eventType=="4-fold"){
            canv->SetCanvasSize(1200,1200);
            for (int i = 0; i<2; i++){
                for (int j = 0; j<2; j++){
                    int index = i*2+j;
                    pads[index] = new TPad(Form("pad%d_%d",i,j),"",0.5*j,0.45*(1-i),0.5*(j+1),0.45*(2-i));
                    pads[index]->Draw();
                    pads[index]->SetGridx(1); pads[index]->SetGridy(1); pads[index]->SetTickx(1); pads[index]->SetTicky(1);
                }
            }
        }
        // delta t
        pads[0]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogy(1);
        if (eventType=="4-fold"){
            legend = new TLegend(0.7,0.7,0.9,0.9);
            t->Draw("evt_tB-evt_tA>>h_dt_AB(51,-20.4,20.4)",cuts,"HIST");
            h1 = (TH1D*) gDirectory->Get("h_dt_AB");
            h1->SetLineColor(kRed);
            h1->SetTitle(";#Deltat [ns];Count");
            legend->AddEntry(h1,"T_{0R}-T_{0L}");
            maximum = 0;
            if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
            t->Draw("evt_tD-evt_tC>>h_dt_CD(51,-20.4,20.4)",cuts,"SAME");
            h1 = (TH1D*) gDirectory->Get("h_dt_CD");
            h1->SetLineColor(kBlue);
            legend->AddEntry(h1,"T_{2}-T_{1}");
            if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
            t->Draw("(evt_tD+evt_tC)/2-(evt_tB+evt_tA)/2>>h_dt_ABCD(51,-20.4,20.4)",cuts,"SAME");
            h1 = (TH1D*) gDirectory->Get("h_dt_ABCD");
            h1->SetLineColor(kBlack);
            legend->AddEntry(h1,"T_{12}-T_{0}");
            if (maximum<h1->GetMaximum()) maximum = h1->GetMaximum();
            h1 = (TH1D*) gDirectory->Get("h_dt_AB");
            h1->GetYaxis()->SetRangeUser(0.5,maximum*2);
            legend->Draw("SAME");
        }
        // 2-D size: T12 VS T0
        pads[1]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogz(1);
        t->Draw(Form("(evt_%cC+evt_%cD)/2:(evt_%cA+evt_%cB)/2>>h_s_ABCD(%d,%.16e,%.16e,%d,%.16e,%.16e)",sizeChar,sizeChar,sizeChar,sizeChar,size_nbins[0],size_min[0],size_max[0],size_nbins[2],size_min[2],size_max[2]),cuts,"COLZ");
        // 2-D size: T0
        pads[2]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogz(1);
        t->Draw(Form("evt_%cB:evt_%cA>>h_s_AB(%d,%.16e,%.16e,%d,%.16e,%.16e)",sizeChar,sizeChar,size_nbins[0],size_min[0],size_max[0],size_nbins[0],size_min[0],size_max[0]),cuts,"COLZ");
        // 2-D size: T12
        pads[3]->cd(); gPad->SetGridx(1); gPad->SetGridy(1); gPad->SetTickx(1); gPad->SetTicky(1); gPad->SetLogz(1);
        t->Draw(Form("evt_%cD:evt_%cC>>h_s_CD(%d,%.16e,%.16e,%d,%.16e,%.16e)",sizeChar,sizeChar,size_nbins[2],size_min[2],size_max[2],size_nbins[2],size_min[2],size_max[2]),cuts,"COLZ");
        canv->cd();
        text->Draw();
        canv->SaveAs(Form("coin.%s.%c.%s%s.png",eventType.Data(),sizeChar,runname.Data(),suffix.Data()));
    }
}
