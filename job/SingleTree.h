//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 10 04:38:04 2013 by ROOT version 5.26/00e
// from TTree SingleTree/SingleTree
// found on file: /afs/ihep.ac.cn/users/l/lidj/largedata/IsotopesAna/P12E/EH1/run21868_IsotopesAna.root
//////////////////////////////////////////////////////////

#ifndef SingleTree_h
#define SingleTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1F.h>
#include <TH2F.h>
#include  <TCanvas.h>
#include  <TLegend.h>
#include  <string>
#include  <TString.h>

class SingleTree : public TSelector {
    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain

        TString histname;
        TString histname2;
        //char nameChar[100];
        int binNum;
        float LowEdge;
        float HighEdge;
        float LowEdge4e;
        float HighEdge4e;
        double signalWinLow;
        double signalWinHigh;
        double offWinLow;
        double offWinHigh;
        double signalWinLowI;
        double signalWinHighI;
        double offWinLowI;
        double offWinHighI;
        string IsoMode;
        TFile* file;
        TString dataVer;
        int ADNum;
        int site;
        TH2F* singleSpecVsTime[4];
        TH1F* signalWin[4];
        TH1F* offWin[4];
        int offTheoNum[4];
        int offRealNum[4];
        bool isRealOff;
        TH2F* signalWinXY[4];
        TH2F* offWinXY[4];
        TH2F* signalWinRZ[4];
        TH2F* offWinRZ[4];
        TH1F* isoSpec[4];
        TH1F* singleUpper[4];
        TH1F* singleLower[4];
        TTree* time2lastmuon[4];
        TH1F* time2lastshowermuon[4];
        int totalEntries;
   // Declaration of leaf types
   Int_t           det;
   Float_t         energy[2];
   Float_t         x[2];
   Float_t         y[2];
   Float_t         z[2];
   Double_t        timeInterval;
   Double_t        promptT2Muon[9];

   // List of branches
   TBranch        *b_det_l;   //!
   TBranch        *b_energy_l;   //!
   TBranch        *b_x_l;   //!
   TBranch        *b_y_l;   //!
   TBranch        *b_z_l;   //!
   TBranch        *b_timeInterval;   //!
   TBranch        *b_promptT2Muon;   //!

   SingleTree(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~SingleTree() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(SingleTree,0);
};

#endif

#ifdef SingleTree_cxx
void SingleTree::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("det", &det, &b_det_l);
   fChain->SetBranchAddress("energy", energy, &b_energy_l);
   fChain->SetBranchAddress("x", x, &b_x_l);
   fChain->SetBranchAddress("y", y, &b_y_l);
   fChain->SetBranchAddress("z", z, &b_z_l);
   fChain->SetBranchAddress("timeInterval", &timeInterval, &b_timeInterval);
   fChain->SetBranchAddress("promptT2Muon", promptT2Muon, &b_promptT2Muon);
}

Bool_t SingleTree::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

#endif // #ifdef SingleTree_cxx
