{
    TFile* f=new TFile("EH2Li_P12E.root");
    TH1F* signalWin=(TH1F*)f->Get("signalWindow");
    TH1F* offWin=(TH1F*)f->Get("offWindow");
    TH1F* LiSpec;
    
    //TCanvas *c1=new TCanvas("c1","c1",800,600);
    signalWin->SetLineColor(kRed);
    signalWin->GetXaxis()->SetTitle("Energy(MeV)");
    signalWin->GetYaxis()->SetTitle("Entries");
    signalWin->SetStats(kFALSE); 
    signalWin->Draw("hist");
    offWin->SetLineColor(kGreen);
    offWin->SetStats(kFALSE); 
    offWin->Draw("histsame");
    TH1F* LiSpec=new TH1F(*signalWin);
    LiSpec->Sumw2();
    LiSpec->Add(signalWin,offWin,1,-1);
    LiSpec->SetLineColor(kBlue);
    LiSpec->SetStats(kFALSE);
    //LiSpec->Draw("EPsame");
    LiSpec->Draw("same");
    TLegend *legend=new TLegend(.6,.65,.79,.89);
    legend->AddEntry(signalWin,"signal window","lp");
    legend->AddEntry(off,"off window","lp");
    legend->AddEntry(LiSpec,"^{9}Li/^{8}He","lp");
    legend->SetFillColor(0);
    legend->Draw("same");
    //signalWin->Write();
    //offWin->Write();
    //LiSpec->Write();
    //c1->Write();
    //f->Close();
}
