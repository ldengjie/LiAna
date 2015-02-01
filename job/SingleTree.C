#define SingleTree_cxx

#include "SingleTree.h"
#include <TH2.h>
#include <TStyle.h>
#include  <THStack.h>
#include  <iostream>

void SingleTree::Begin(TTree * /*tree*/)
{
    TString option = GetOption();
    TString IsoModeTmp=option(7,3);
    TString evtNumStr=option(10,(option.Sizeof()-10));
    totalEntries=evtNumStr.Atoi();
    cout<<"totalEntries  : "<<totalEntries<<endl;
    //IsoMode="C9";//Li8,C99,B12
    IsoMode=IsoModeTmp;//Li8,C99,B12
    binNum=80;
    for( int i=0 ; i<9 ; i++ )
    {
        promptT2Muon[i]=0.;
    }

    if( IsoMode=="Li9" )
    {
        LowEdge=3.5;
        //LowEdge=4.0;
        HighEdge=20.0;
        LowEdge4e=LowEdge;
        HighEdge4e=HighEdge;
        signalWinLow=0.002;
        signalWinHigh=1.000;
        offWinLow=1.002;
        offWinHigh=2.000;
        signalWinLowI=0.002;
        signalWinHighI=0.100;
        offWinLowI=-0.098;
        offWinHighI=-0.0;
    }

    TString DataVerTmp=option(3,4);
    dataVer=DataVerTmp;
    option=option(0,3);
    std::cout<<"dataVer  : "<<dataVer<<endl;
    int daqHistNum;
    if( dataVer=="P12E"||dataVer=="P12C" )
    {
        daqHistNum=4;
    } else
    {
        daqHistNum=5;
    }
    std::cout<<"option  : "<<option<<endl;
    ADNum=daqHistNum-1;
    site=2;
    if( option=="EH1")
    {
        site=0;
        ADNum=2;
    }
    if( option=="EH2" )
    {
        site=1;
        ADNum=daqHistNum-3;
    }
    option+="Li9_";
    option+=dataVer;
    option+=".root";
    option=dataVer+"/"+option;
    file = new TFile(option,"update");

    //histogram - binned data for MML or Chi2
    for( int i=0 ; i<4 ; i++ )
    {
        histname=Form("time2lastshowermuon%i_%0.1f_%0.1f",i+1,LowEdge,HighEdge);
        histname2=histname+";*";
        gDirectory->Delete(histname2);
        time2lastshowermuon[i]=new TH1F(histname,"time2lastshowermuon",99999,0.001,100);
    }


    //tree - unbinned data for MML
    for( int i=0 ; i<4 ; i++ )
    {
        histname=Form("slice%i_%0.1f_%0.1f",i+1,LowEdge,HighEdge);
        histname2=histname+";*";
        gDirectory->Delete(histname2);
        time2lastmuon[i]=new TTree(histname,"time2lastmuon");
        time2lastmuon[i]->Branch("xt",&promptT2Muon[4+i],"xt/D");
    }

    //energy[0] spectrum / xy / rz /specVsTime  after muon w/o reduction cut ( there are neutrons after it).
    for( int i=0 ; i<4 ; i++ )
    {
        offTheoNum[i]=1;
        offRealNum[i]=1;
        histname=Form("%ssignalEnergySlice%i_%0.1f_%0.1f",IsoMode.c_str(),i+1,LowEdge,HighEdge);
        histname2=histname+";*";
        gDirectory->Delete(histname2);
        signalWin[i]=new TH1F(histname,histname,binNum,0,20);
        histname=Form("%soffEnergySlice%i_%0.1f_%0.1f",IsoMode.c_str(),i+1,LowEdge,HighEdge);
        histname2=histname+";*";
        gDirectory->Delete(histname2);
        offWin[i]=new TH1F(histname,histname,binNum,0,20);
        histname=Form("%ssignalXYSlice%i_%0.1f_%0.1f",IsoMode.c_str(),i+1,LowEdge,HighEdge);
        histname2=histname+";*";
        gDirectory->Delete(histname2);
        std::cout<<"histname  : "<<histname<<endl;
        signalWinXY[i]=new TH2F(histname,histname,600,-3000.,3000.,600,-3000.,3000.);
        histname=Form("%soffXYSlice%i_%0.1f_%0.1f",IsoMode.c_str(),i+1,LowEdge,HighEdge);
        histname2=histname+";*";
        gDirectory->Delete(histname2);
        offWinXY[i]=new TH2F(histname,histname,600,-3000,3000,600,-3000,3000);
        histname=Form("%ssignalRZSlice%i_%0.1f_%0.1f",IsoMode.c_str(),i+1,LowEdge,HighEdge);
        histname2=histname+";*";
        gDirectory->Delete(histname2);
        signalWinRZ[i]=new TH2F(histname,histname,300,0,3000,600,-3000,3000);
        histname=Form("%soffRZSlice%i_%0.1f_%0.1f",IsoMode.c_str(),i+1,LowEdge,HighEdge);
        histname2=histname+";*";
        gDirectory->Delete(histname2);
        offWinRZ[i]=new TH2F(histname,histname,300,0,3000,600,-3000,3000);
        histname=Form("%ssingleSpecVsTimeSlice%i_%0.1f_%0.1f",IsoMode.c_str(),i+1,LowEdge,HighEdge);
        histname2=histname+";*";
        gDirectory->Delete(histname2);
        singleSpecVsTime[i]=new TH2F(histname,histname,100000,0,100,400,0,20);
    }
    std::cout<<"finished Begin() "<<endl;
}

void SingleTree::SlaveBegin(TTree * /*tree*/)
{

    TString option = GetOption();

}

Bool_t SingleTree::Process(Long64_t entry)
{
    GetEntry(entry);
    //if( entry%10000==0 )
    //{
    //std::cout<<entry<<endl;
    //}

    //cut all events from Z>=900;cut flahser that exists in P12E P13A  from this zone;
    //if(x[0]>-1700&&x[0]<-1200&&y[0]>-1500&&y[0]<-900&&z[0]>-900&&x[0]<-500)
    //{
    //return true;
    //}

    //isotopes spectrum
    if( energy[0]>=LowEdge4e&&energy[0]<=HighEdge4e )
    {
        for( int j=0 ; j<3 ; j++ )
        {
            singleSpecVsTime[j]->Fill(promptT2Muon[j+4],energy[0]);
            if( promptT2Muon[j+4]>=signalWinLow &&promptT2Muon[j+4]<=signalWinHigh )
            {
                signalWinXY[j]->Fill(x[0],y[0]);
                signalWinRZ[j]->Fill(sqrt(x[0]*x[0]+y[0]*y[0]),z[0]);
                signalWin[j]->Fill(energy[0]);
            }
            if( promptT2Muon[j+4]>=offWinLow &&promptT2Muon[j+4]<=offWinHigh )
            {
                offWinXY[j]->Fill(x[0],y[0]);
                offWinRZ[j]->Fill(sqrt(x[0]*x[0]+y[0]*y[0]),z[0]);
                offTheoNum[j]++;
                isRealOff=1; 
                for( int a=0 ; a<5 ; a++ )
                {
                    if( promptT2Muon[a+4]>=signalWinLow&&promptT2Muon[a+4]<=signalWinHigh )
                    {
                        isRealOff=0;    
                        break;
                    }
                }
                //if( isRealOff )
                //{
                offWin[j]->Fill(energy[0]);
                offRealNum[j]++;
                //} 
            }
        }

        if( promptT2Muon[7]>=signalWinLowI &&promptT2Muon[7]<=signalWinHighI )
        {
            signalWin[3]->Fill(energy[0]);
        }
        if(promptT2Muon[8]>=offWinLowI && promptT2Muon[8]<=offWinHighI)
        {
            offWin[3]->Fill(energy[0]);
        }
    }
    if( energy[0]<LowEdge || energy[0]>HighEdge)
    {
        return 1; 
    }

    for( int j=0 ; j<4 ; j++ )
    {
        time2lastshowermuon[j]->Fill(promptT2Muon[j+4]);

    }
    for( int j=0 ; j<4 ; j++ )
    {
        time2lastmuon[j]->Fill();
    }
    return kTRUE;
}

void SingleTree::SlaveTerminate()
{

}

void SingleTree::Terminate()
{
    std::cout<<"Now is in Terminate() "<<endl;
    //histname=Form("%s/EH%itotalHisto_%s.root",dataVer.Data(),site+1,dataVer.Data());
    //TFile f(histname,"read");
    //if( f.IsZombie() )
    //{
    //std::cout<<"Error : can't open "<<histname<<endl;
    //return;
    //}
    for( int i=0 ; i<3 ; i++ )
    {
        /*
        histname=Form("lidj/muonTimeInterval%s%i",i+1);
        TH1F* h=(TH1F*)f.Get(histname);
        if( !h )
        {
            std::cout<<"Error : can't open TH1F "<<histname<<endl;
            return;
        }
        file->cd();
        int offWinNum=(int)((offWinHigh-offWinLow)/(signalWinHigh-signalWinLow));
        if(offWinNum==0) offWinNum=1;

        double offMuonRate=h->Integral(h->FindBin(offWinLow),h->FindBin(9999));
        double offMuonRateH=h->Integral(h->FindBin(offWinHigh),h->FindBin(9999));

        double offFrac=0.;
        double offFracH=0.;
        for( int k=0 ; k<offWinNum ; k++ )
        {
            double offMuonRateTmp=h->Integral(h->FindBin(offWinLow+(signalWinHigh-signalWinLow)*k),h->FindBin(9999));
            double offMuonRateTmpH=h->Integral(h->FindBin(offWinHigh+(signalWinHigh-signalWinLow)*k),h->FindBin(9999));

            offFrac+=offMuonRateTmp/offMuonRate;
            offFracH+=offMuonRateTmpH/offMuonRateH;
        }

        double signalMuonRate=h->Integral(h->FindBin(signalWinLow),h->FindBin(9999));
        double signalMuonRateH=h->Integral(h->FindBin(signalWinHigh),h->FindBin(9999));
        std::cout<<"signalFrac  : "<<signalMuonRate/offMuonRate<<endl;
        std::cout<<"signalFracH  : "<<signalMuonRateH/offMuonRateH<<endl;
        cout<<"offFrac  : "<<offFrac<<endl;
        cout<<"offFracH  : "<<offFracH<<endl;
        histname=Form("%s Slice%i signalMuonRate:%10.0f offMuonRate:%10.0f signal/off:%5.5f offTheo/Real:%5.5f frac:%.5f",IsoMode.c_str(),i+1,signalMuonRate,offMuonRate,signalMuonRate/offMuonRate,(double)(offTheoNum[i])/offRealNum[i],offFrac);
        std::cout<<histname<<endl;
        delete h;
        */
        signalWin[i]->SetLineColor(kRed);
        signalWin[i]->GetXaxis()->SetTitle("Energy(MeV)");
        signalWin[i]->GetYaxis()->SetTitle("Entries");
        signalWin[i]->SetStats(kFALSE); 
        //signalWin[i]->Draw();
        offWin[i]->SetLineColor(kGreen);
        offWin[i]->SetStats(kFALSE); 
        //offWin[i]->Draw("same");
        histname=Form("%sSpecSlice%i_%0.1f_%0.1f",IsoMode.c_str(),i+1,LowEdge4e,HighEdge4e);
        histname2=histname+";*";
        gDirectory->Delete(histname2);
        isoSpec[i]=new TH1F(histname,histname,binNum,0,20);
        isoSpec[i]->Sumw2();
        isoSpec[i]->Add(signalWin[i],offWin[i],1,-1);
        //isoSpec[i]->Add(signalWin[i],offWin[i],1,-(1/(offFrac*offMuonRate/signalMuonRate*offRealNum[i]/offTheoNum[i])));
        isoSpec[i]->SetLineColor(kBlue);
        isoSpec[i]->SetStats(kFALSE);
        isoSpec[i]->SetMarkerStyle(20);
        isoSpec[i]->SetMarkerSize(0.7);
        isoSpec[i]->SetMarkerColor(kBlue);
        isoSpec[i]->SetOption("E");
        signalWin[i]->Write();
        offWin[i]->Write();
        isoSpec[i]->Write();
        singleSpecVsTime[i]->Write();
        file->cd();
        signalWinXY[i]->Write();
        offWinXY[i]->Write();
        signalWinRZ[i]->Write();
        offWinRZ[i]->Write();
        //c[i]->Write();

    }

    //signalWin[3]->Draw();
    //offWin[3]->Draw();
    histname=Form("%sSpecSlice%i_%0.1f_%0.1f",IsoMode.c_str(),4,LowEdge4e,HighEdge4e);
    histname2=histname+";*";
    gDirectory->Delete(histname2);
    isoSpec[3]=new TH1F(histname,histname,binNum,0,20);
    isoSpec[3]->Sumw2();
    isoSpec[3]->Add(signalWin[3],offWin[3],1,-1);
    //isoSpec[3]->Draw();

    offWin[3]->Write();
    signalWin[3]->Write();
    isoSpec[3]->Write();


    for( int i=0 ; i<4 ; i++ )
    {
        time2lastshowermuon[i]->Write();
    }

    //c1->Write();
    file->Close();
}
