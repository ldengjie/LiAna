//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 10 04:38:04 2013 by ROOT version 5.26/00e
// from TTree SingleTree/SingleTree
// found on file: /afs/ihep.ac.cn/users/l/lidj/largedata/LiAna/P12E/EH1/run21868_LiAna.root
//////////////////////////////////////////////////////////

#ifndef SingleTree_h
#define SingleTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1F.h>
#include  <TCanvas.h>
#include  <TLegend.h>

class SingleTree : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   TString histname;
   TFile* file;
   TString dataVer;
   int ADNum;
   int site;
    TH1F* singleSpec[4];
    TH1F* signalWin;
    TH1F* offWin;
   // Declaration of leaf types
   Float_t         energy;
   Float_t         x;
   Float_t         y;
   Float_t         z;
   Int_t           det;
   Float_t         t2lastshowermuon;

   // List of branches
   TBranch        *b_energy;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_det;   //!
   TBranch        *b_t2lastshowermuon;   //!

   SingleTree(TTree * /*tree*/ =0) { }
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

   fChain->SetBranchAddress("energy", &energy, &b_energy);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("det", &det, &b_det);
   fChain->SetBranchAddress("t2lastshowermuon", &t2lastshowermuon, &b_t2lastshowermuon);
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
