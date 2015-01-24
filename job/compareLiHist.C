{
    TCanvas *c =new TCanvas("d","d",800,600);
    TFile *f1 = new TFile("/workfs/dyw/lidj/IBDSel/job/EH1Li_P12A_P12E.root");
    //TFile *f2 = new TFile("/publicfs/dyb/data/userdata/zhangfh/IbdAna/ana/HeLiAna/P12E/EH1/mergedHist.root");
    TFile *f2 = new TFile("/afs/ihep.ac.cn/users/l/lidj/file/LiAna/job/EH1totalHisto_P12A_P12E.root");
    TH1F* h1[6];
    TH1F* h2[6];
    TH1F* h3[6];
    TCanvas *c1=new TCanvas("c1","c1",800,900);
    c1->Divide(2,3);
    for( int i=0 ; i<6 ; i++ )
    {
        TString h1Name="time2lastshowermuon";
        h1Name+=i+1;
        h1[i]=(TH1F*)f1->Get(h1Name);
        std::cout<<" "<<endl;
        std::cout<<"h1 0-40  : "<<h1[i]->Integral(1,h1[i]->FindBin(40))<<endl;
        std::cout<<"h1 40  : "<<h1[i]->FindBin(40)<<endl;
        //TString h2Name="zhangfh/time2LastMuonRed";
        TString h2Name="lidj/time2lastshowermuon";
        h2Name+=i+1;
        h2Name+="4Li";
        h2[i]=(TH1F*)f2->Get(h2Name);
        std::cout<<"h2 40  : "<<h2[i]->FindBin(40)<<endl;
        //std::cout<<"h2 all  : "<<h2[i]->GetEntries()<<endl;
        double h2entries=h2[i]->Integral(1,h2[i]->FindBin(40));
        std::cout<<"h2 0-40 : "<<h2entries<<endl;
        double h2entries3=h2[i]->Integral(h2[i]->FindBin(40),h2[i]->FindBin(250));
        std::cout<<"h2 40-250 : "<<h2entries3<<endl;
        double h2entries2=h2[i]->Integral(1,h2[i]->FindBin(250));
        std::cout<<"h2 0-250 : "<<h2entries2<<endl;
        std::cout<<"h2-h1 num  : "<<h2entries-h1[i]->Integral(1,h1[i]->FindBin(40))<<endl;

        h1[i]->Rebin(10);
        h1[i]->GetXaxis()->SetRange(0,h1[i]->FindBin(100));
        h1[i]->SetMinimum(0);
        h2[i]->Rebin(100);

        std::cout<<"h1 100  : "<<h1[i]->FindBin(100)<<endl;
        std::cout<<"h2 100  : "<<h2[i]->FindBin(100)<<endl;
        TString h3name="h3";
        h3name+=i;
        h3[i]=new TH1F(h3name,"h3",1000,0,100);
        std::cout<<"h3 100  : "<<h3[i]->FindBin(100)<<endl;
        //h3[i]->Add(h1[i],h2[i],-1,1);
        for( int j=0 ; j<=1000 ; j++ )
        {
            h3[i]->SetBinContent(j,(h1[i]->GetBinContent(j)-h2[i]->GetBinContent(j)));
        }
        
        c1->cd(i+1);
        h1[i]->SetLineColor(kRed);
        h1[i]->Draw();
        h2[i]->SetLineColor(kBlue);
        h2[i]->Draw("same");
        h3[i]->SetLineColor(kGreen);
        h3[i]->Draw("same");
    }
    /*
    
    std::cout<<"h1->GetEntries()  : "<<h1->GetEntries()<<endl;
    std::cout<<"h1 binNum  : "<< h1->GetXaxis()->GetNbins()<<endl;
    //h1->Rebin(10);
    std::cout<<"h1 binNum  : "<< h1->GetXaxis()->GetNbins()<<endl;
    h1->SetLineColor(2);
    //std::cout<<"h1 num 12-100  : "<<h1->Integral(155,200,"width")<<endl; 
    //std::cout<<"h1 num 0.7-12  : "<<h1->Integral(108,154,"width")<<endl; 
    //std::cout<<"h1 num  : "<<h1->Integral(1,200,"width")<<endl; 
    //h1->Draw();
    std::cout<<"h2->GetEntries()  : "<<h2->GetEntries()<<endl;
    std::cout<<"h2 binNum  : "<< h2->GetXaxis()->GetNbins()<<endl;
    h2->Rebin(4);
    std::cout<<"h2 binNum  : "<< h2->GetXaxis()->GetNbins()<<endl;
    //h2->Draw("same");
    
    TH1F* h3=new TH1F(*h1);
    h3->Add(h1,h2,-1,1);
    //std::cout<<"h3 num 12-100  : "<<h3->Integral(155,200,"width")<<endl; 
    //std::cout<<"h3 num 0.7-12  : "<<h3->Integral(1,h3->FindBin(12),"width")<<endl; 
    //std::cout<<"h3 num  : "<<h3->Integral(1,200,"width")<<endl; 
    h3->Draw();
    */
    //TH1F h5=0*(*h1);
    /*
    TH1F* h5=new TH1F(*h1);
    h5->Add(h1,h1,0,0);
        
    for( int i=2 ; i<=200 ; i++ )
    {
        h5->SetBinContent(i,(h3->Integral(1,i,"width")-h3->Integral(1,i-1,"width")));
    }
    //h5->Draw("same");
    h5->Draw();
    */
}
