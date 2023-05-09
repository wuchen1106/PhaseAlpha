{
    TChain * t = new TChain("t");
    t->Add("anaout/run0*.pair.root");

    //t->Draw("run:startC-(trigger_tA+trigger_tB)/2>>h(1250,-500,9500,800,300,1100)","foundTrigger&&heightC>35","COLZ");
    //t->Draw("run:ptC-(trigger_tA+trigger_tB)/2>>h(1250,-500,9500,800,300,1100)","foundTrigger&&phC>35","COLZ");

    //t->Draw("run:ptA-(trigger_tA+trigger_tB)/2>>h(1251,-0.4,1000.4,800,300,1100)","foundTrigger","COLZ");
    t->Draw("run:(trigger_tA+trigger_tB)/2>>h(101,319.6,400.4,800,300,1100)","foundTrigger","COLZ");
}
