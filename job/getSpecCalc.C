{
    TFile* f = new TFile("Model1.root");
    TGraph* g=(TGraph*)f->Get("electronFull");
    double *x=g->GetX();
    double *y=g->GetY();
    //g->Draw();
    int n=g->GetN();
    std::cout<<"n  : "<<n<<endl;
    //TH1F* h= new TH1F("h","h",300,0,14);
    //for( int i=0 ; i<n ; i++ )
    //{
    //std::cout<<"x:y "<<x[i]<<" " << y[i]<<endl;
    //h->SetBinContent(i+1,y[i]);
    //}
    //h->Draw();

    TFile* f1=new TFile("specCalc.root","update");
    gDirectory->Delete("LiSpec;*");
    gDirectory->Delete("LiSpec_12bin;*");
    TH1F* h1=(TH1F*)f1->Get("hSpecLiFine");
    TH1F* h2=new TH1F("LiSpec","LiSpec after cor",1200,0,12);
    for( int i=1 ; i<=1200 ; i++ )
    {
        double xmid=h1->GetBinCenter(i);
        double cof=1.;
        double xcor=0.;
        double binLowEdge=h1->GetBinLowEdge(i);
        double binWidth=h1->GetBinWidth(i);
        double binContent=h1->GetBinContent(i);
        int corNum=10000;
        double xnow=0.;
        int newbinNum;
        for( int j=0 ; j<300 ; j++ )
        {
            if( xmid>x[j]&&xmid<x[j+1] )
            {
                for( int k=0 ; k<corNum ; k++ )
                {
                    xnow=binLowEdge+(binWidth/corNum)*(1/2+k);
                    cof=y[j]+(y[j+1]-y[j])*(xnow-x[j])/(x[j+1]-x[j]);
                    xcor=(xnow)*cof;
                    newbinNum=h2->FindBin(xcor);
                    h2->AddBinContent(newbinNum,binContent/corNum);
                }
                break;
            }
        }
        //if( i<50 )
        //{
        //std::cout<<"old -> new "<<i<<" -> "<<newbinNum<<endl;
        //}
        
    }
    TH1F* h3=new TH1F("LiSpec_12bin","LiSpec after cor",12,0,12);
    for( int i=0 ; i<12 ; i++ )
    {
        double binContent=0.;
        for( int j=i*100 ; j<(i+1)*100 ; j++ )
        {
            if(j>70) binContent+=h2->GetBinContent(j+1);
        }
        h3->SetBinContent(i+1,binContent);
    }
    h3->Scale(1/h3->Integral());
    h3->SetLineColor(kGreen);
    h3->Draw();
    h3->Write();
    h1->SetLineColor(kRed);
    h1->Draw("same");
    h2->SetLineColor(kBlue);
    h2->Draw("same");
    h2->Write();
    //f1->Close();
    
}
