#define LiTree_cxx
// The class definition in LiTree.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("LiTree.C")
// Root > T->Process("LiTree.C","some options")
// Root > T->Process("LiTree.C+")
//

#include "LiTree.h"
#include <TH2.h>
#include  <iostream>
#include <TStyle.h>


void LiTree::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();


   file = new TFile(option,"RECREATE");
   siteStr=option(0,3);

    histname="signalWindow";
    signalWin=new TH1F(histname,"signalwidow",24,0,12);
    histname="offWindow";
    offWin=new TH1F(histname,"offwidow",24,0,12);
    histname="signalWindowWithMR";
    signalWinWithMR=new TH1F(histname,"signalwidowWithMR",24,0,12);
    histname="offWindowWithMR";
    offWinWithMR=new TH1F(histname,"offwidowWithMR",24,0,12);
}

void LiTree::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t LiTree::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either LiTree::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
    GetEntry(entry);
    /*
        for( int j=0 ; j<6 ; j++ )
        {
            t2lastshowermuon[4][j]->Fill(t2lastshowermuon);
            t2lastshowermuon[det-1][j]->Fill(t2lastshowermuon);
        }
        */
        if( !(((x[0]-x[1])*(x[0]-x[1])+(y[0]-y[1])*(y[0]-y[1])+(z[0]-z[1])*(z[0]-z[1])<1000*1000) && (timeInterval<100)) )
        {
            return kTRUE;
        }
        if( t2lastshowermuon>1.e-3 &&t2lastshowermuon<1. )
        {
            signalWin->Fill(energy[0]);
        }
        if( t2lastshowermuon>1. &&t2lastshowermuon<2. )
        {
            offWin->Fill(energy[0]);
        }
        if( t2lastshowermuonWithMR>1.e-3 &&t2lastshowermuonWithMR<1. )
        {
            signalWinWithMR->Fill(energy[0]);
        }
        if( t2lastshowermuonWithMR>1. &&t2lastshowermuonWithMR<2. )
        {
            offWinWithMR->Fill(energy[0]);
        }



   return kTRUE;
}

void LiTree::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void LiTree::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

    double scal=0.;
    if( siteStr=="EH1" )
    {
        scal=0.064;
    }
    if( siteStr=="EH2" )
    {
        scal=0.055;
    }
    if( siteStr=="EH3" )
    {
        scal=0.0055;
    }
    double scale=(exp(-scal*0.001)-exp(-scal))/(exp(-scal)-exp(-scal*2));

    signalWinWithMR->Rebin(2);
    offWinWithMR->Rebin(2);
    signalWinWithMR->Sumw2();
    offWinWithMR->Sumw2();
    TH1F* LiSpecWithMR=new TH1F("LiSpecWithMR","LiSpecWithMR",12,0,12);
    LiSpecWithMR->Sumw2();
    LiSpecWithMR->Add(signalWinWithMR,offWinWithMR,1,-scale);
    TH1F* LiSpecScaledWithMR=new TH1F("LiSpecScaledWithMR","LiSpec after scaled",12,0,12);
    LiSpecScaledWithMR->Add(signalWinWithMR,LiSpecWithMR,0,1);
    LiSpecScaledWithMR->Scale(1/LiSpecScaledWithMR->Integral());
    signalWinWithMR->Write();
    offWinWithMR->Write();
    LiSpecWithMR->Write();
    LiSpecScaledWithMR->Write();
    //TCanvas *c1=new TCanvas("c1","c1",1200,400);
    //c1->Divide(2,1);
    //c1->cd(1);
    signalWin->SetLineColor(kRed);
    signalWin->GetXaxis()->SetTitle("Energy(MeV)");
    signalWin->GetYaxis()->SetTitle("Entries");
    signalWin->SetStats(kFALSE); 
    signalWin->SetMarkerStyle(4);
    signalWin->SetMarkerSize(0.5);
    signalWin->SetMarkerColor(kRed);
    signalWin->Rebin(2);
    //signalWin->Draw("e");
    offWin->SetLineColor(kGreen);
    offWin->SetStats(kFALSE); 
    offWin->SetMarkerStyle(4);
    offWin->SetMarkerSize(0.5);
    offWin->SetMarkerColor(kGreen);
    offWin->Rebin(2);
    //offWin->Draw("esame");
    signalWin->Sumw2();
    offWin->Sumw2();
    TH1F* LiSpec=new TH1F("LiSpec","LiSpec",12,0,12);
    LiSpec->Sumw2();
    LiSpec->Add(signalWin,offWin,1,-scale);
    double eventNum=0.;
    for( int i=1 ; i<=12 ; i++ )
    {
        eventNum+=LiSpec->GetBinContent(i);
    }
    std::cout<<"LiSpec eventNumber  : "<<eventNum<<endl;
    LiSpec->SetLineColor(kBlue);
    LiSpec->SetStats(kFALSE);
    LiSpec->SetMarkerStyle(20);
    LiSpec->SetMarkerSize(0.7);
    LiSpec->SetMarkerColor(kBlue);
    //LiSpec->Draw("sameEP");
    TLegend *legend=new TLegend(.6,.65,.79,.89);
    legend->AddEntry(signalWin,"signal window","lp");
    legend->AddEntry(offWin,"off window","lp");
    legend->AddEntry(LiSpec,"^8He/^9Li","lp");
    legend->SetFillColor(0);
    //legend->Draw("same");
    //c1->cd(2);
    TH1F* LiSpecScaled=new TH1F("LiSpecScaled","LiSpec after scaled",12,0,12);
    LiSpecScaled->Add(signalWin,LiSpec,0,1);
    LiSpecScaled->Scale(1/LiSpecScaled->Integral());
    LiSpecScaled->SetLineColor(kRed);
    LiSpecScaled->SetStats(kFALSE);
    LiSpecScaled->SetMarkerStyle(20);
    LiSpecScaled->SetMarkerSize(0.7);
    LiSpecScaled->SetMarkerColor(kRed);
    //LiSpecScaled->Draw("EP");
    TFile* SpecCalFile=new TFile("specCalc_12bin.root");
    TH1F* SpecCal=(TH1F*)SpecCalFile->Get("hSpecLi");
    SpecCal->SetLineColor(kBlue);
    //SpecCal->Draw("same");
    file->cd();
    signalWin->Write();
    offWin->Write();
    LiSpec->Write();
    LiSpecScaled->Write();
    //std::cout<<"signal number  : "<<signalWin->Integral(0,24)<<endl;
    //std::cout<<"off number  : "<<offWin->Integral(0,24)<<endl;
    //std::cout<<"Li number  : "<<LiSpec->Integral(0,24)<<endl;
}
