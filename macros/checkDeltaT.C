void checkDeltaT(TString name=""){
    double hmin[4] = {35,35,25,25};
    TChain * ichain = new TChain("t");
    ichain->Add(name);
    std::cout<<name<<": "<<ichain->GetEntries("foundTrigger")<<"/"<<ichain->GetEntries()<<" = "<<100.*ichain->GetEntries("foundTrigger")/ichain->GetEntries()<<"%"<<std::endl;
    TCanvas * canv = new TCanvas("canv","");
    canv->SetCanvasSize(1600,800);
    canv->Divide(4,2);
    TString otherCuts = "&&evt_hA>0&&evt_hB>0&&evt_hC>0&&evt_hD>0";
    for (int ch = 1; ch<4; ch++){
        canv->cd(ch);
        ichain->Draw(Form("evt_tA-evt_t%c>>htA%c(51,-20.4,20.4)",'A'+ch,'A'+ch),Form("evt_hA>%.1f&&evt_h%c>%.1f",hmin[0],'A'+ch,hmin[ch])+otherCuts,"hist");
    }
    canv->cd(4);
    ichain->Draw(Form("evt_tC-evt_tD>>htCD(51,-20.4,20.4)"),Form("evt_hC>%.1f&&evt_hD>%.1f",hmin[2],hmin[3])+otherCuts,"hist");
    for (int ch = 1; ch<4; ch++){
        canv->cd(4+ch);
        ichain->Draw(Form("evt_h%c:evt_hA>>hhA%c(100,0,1000,100,0,1000)",'A'+ch,'A'+ch),Form("evt_hA>0&&evt_h%c>0",'A'+ch)+otherCuts,"COLZ");
    }
    canv->cd(8);
    ichain->Draw(Form("evt_hC:evt_hD>>hhCD(100,0,1000,100,0,1000)"),Form("evt_hC>0&&evt_hD>0")+otherCuts,"COLZ");
}
