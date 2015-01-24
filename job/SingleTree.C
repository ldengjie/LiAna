#define SingleTree_cxx
// The class definition in SingleTree.h has been generated automatically
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
// Root > T->Process("SingleTree.C")
// Root > T->Process("SingleTree.C","some options")
// Root > T->Process("SingleTree.C+")
//

#include "SingleTree.h"
#include <TH2.h>
#include <TStyle.h>
#include  <iostream>

void SingleTree::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

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
   option+="iso_";
   option+=dataVer;
   option+=".root";
   file = new TFile(option,"RECREATE");
   for( int i=0 ; i<4 ; i++ )
   {
        histname="singleSpecAD";
        histname+=i+1;
        singleSpec[i]=new TH1F(histname,"spectra of singles",400,0,20);
   }
   
    histname="signalWindow";
    signalWin=new TH1F(histname,"signalwidow",80,0,20);
    histname="offWindow";
    offWin=new TH1F(histname,"offwidow",80,0,20);

}

void SingleTree::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t SingleTree::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either SingleTree::GetEntry() or TBranch::GetEntry()
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
    singleSpec[det-1]->Fill(energy);
    if( energy<3.0 )
    {
        return kTRUE;
    }
        if( t2lastshowermuon>1.e-3 &&t2lastshowermuon<0.1 )
        {
            if( (x*x+y*y<1550*1550)&&(z*z<1550*1550) )
            {
                signalWin->Fill(energy);
            }
        }
        if( t2lastshowermuon>0.302 &&t2lastshowermuon<0.5 )
        {
            if( (x*x+y*y<1550*1550)&&(z*z<1550*1550) )
            {
                offWin->Fill(energy);
            }
        }


   return kTRUE;
}

void SingleTree::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void SingleTree::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

    //TCanvas *c1=new TCanvas("c1","c1",800,600);
    signalWin->SetLineColor(kRed);
    signalWin->GetXaxis()->SetTitle("Energy(MeV)");
    signalWin->GetYaxis()->SetTitle("Entries");
    signalWin->SetStats(kFALSE); 
    //signalWin->Draw();
    offWin->SetLineColor(kGreen);
    offWin->SetStats(kFALSE); 
    //offWin->Draw("same");
    TH1F* B12Spec=new TH1F("B12Spec","B12Spec",80,0,20);
    B12Spec->Sumw2();
    B12Spec->Add(signalWin,offWin,1,-0.5);
    B12Spec->SetLineColor(kBlue);
    B12Spec->SetStats(kFALSE);
    B12Spec->SetMarkerStyle(20);
    B12Spec->SetMarkerSize(0.7);
    B12Spec->SetMarkerColor(kBlue);
    //B12Spec->Draw("sameE");
    TLegend *legend=new TLegend(.6,.65,.79,.89);
    legend->AddEntry(signalWin,"signal window","lp");
    legend->AddEntry(offWin,"off window","lp");
    legend->AddEntry(B12Spec,"^{12}B","lp");
    legend->SetFillColor(0);
    //legend->Draw("same");
    singleSpec[0]->SetOption("E");
    singleSpec[0]->SetLineColor(kBlue);
    singleSpec[0]->Draw();
    singleSpec[1]->SetOption("E");
    singleSpec[1]->SetLineColor(kRed);
    singleSpec[1]->Draw("same");
    for( int i=0 ; i<4 ; i++ )
    {
        singleSpec[i]->Write();
    }
    
    signalWin->Write();
    offWin->Write();
    B12Spec->Write();
    //c1->Write();
    //file->Close();
}
