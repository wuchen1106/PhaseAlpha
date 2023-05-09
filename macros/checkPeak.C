{
    double ps_dt = 0.8; // ns

    TChain * ichain_raw = new TChain("tr","tr");
    ichain_raw->Add("decdata/rawdata00397.root");
    TChain * ichain_ana = new TChain("t","t");
    ichain_ana->Add("anaout/run00397.ana.root");

    std::vector< std::vector<Float_t>* > pwfs(4, nullptr);
    double ped[4] = {0};
    int    nPeaks[4] = {0};
    std::vector<int> * start[4] = {0};
    std::vector<int> * width[4] = {0};
    std::vector<double> * height[4] = {0};
    std::vector<double> * integral[4] = {0};
    for (int ch = 0; ch<4; ch++){
        ichain_raw->SetBranchAddress(Form("pwf%c", 'A'+ch), &pwfs[ch]);
        ichain_ana->SetBranchAddress(Form("ped%c",'A'+ch),&ped[ch]);
        ichain_ana->SetBranchAddress(Form("nPeaks%c",'A'+ch),&nPeaks[ch]);
        ichain_ana->SetBranchAddress(Form("start%c",'A'+ch),&start[ch]);
        ichain_ana->SetBranchAddress(Form("width%c",'A'+ch),&width[ch]);
        ichain_ana->SetBranchAddress(Form("height%c",'A'+ch),&height[ch]);
        ichain_ana->SetBranchAddress(Form("integral%c",'A'+ch),&integral[ch]);
    }

    int iEntry = 0;
    ichain_ana->GetEntry(iEntry);
    ichain_raw->GetEntry(iEntry);

    gStyle->SetOptStat(0);
    TCanvas * canvas[4];
    for (int ch = 0; ch<4; ch++){
        canvas[ch] = new TCanvas(Form("canv%d",ch),"",1024,768);
        TH1D * hwf = new TH1D(Form("hwf%d",ch),Form("Picoscope ch%c, event %d, pedestal %.1f mV, %d peaks;TDC [%.1f ns];[mV]",'A'+ch,iEntry,ped[ch],nPeaks[ch],ps_dt),pwfs[ch]->size(),-0.5,pwfs[ch]->size()-0.5);
        for (int i = 0; i<pwfs[ch]->size(); i++){
            hwf->SetBinContent(i+1,pwfs[ch]->at(i));
        }
        hwf->Draw();
        for (int i = 0; i<nPeaks[ch]; i++){
            TLatex * text = new TLatex(start[ch]->at(i),ped[ch]-height[ch]->at(i),Form("%d: %d %.1f mV",i,start[ch]->at(i),height[ch]->at(i)));
            text->SetTextFont(42);
            text->SetTextSize(0.023);
            text->Draw();
        }
    }

}
