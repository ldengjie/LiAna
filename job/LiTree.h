//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Dec  8 23:29:50 2013 by ROOT version 5.26/00e
// from TTree LiTree/LiTree
// found on file: /afs/ihep.ac.cn/users/l/lidj/largedata/LiAna/P12E/EH1/run22405_LiAna.root
//////////////////////////////////////////////////////////

#ifndef LiTree_h
#define LiTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1F.h>
#include  <TCanvas.h>
#include  <TLegend.h>

class LiTree : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   TString histname;
   TFile* file;
   TString dataVer;
   TString siteStr;
   int ADNum;
   int site;
   //TH1F* t2lastshowermuon[5][6];
    TH1F* signalWin;
    TH1F* offWin;
    TH1F* signalWinWithMR;
    TH1F* offWinWithMR;
   // Declaration of leaf types
   Int_t           det;
   Float_t         energy[2];
   Float_t         x[2];
   Float_t         y[2];
   Float_t         z[2];
   double         timeInterval;
   double         t2lastshowermuon;
   double         t2lastshowermuonWithMR;
   //float         t2lastshowermuonWithMR;
   //float timeInterval;
   //float         t2lastshowermuon;

   // List of branches
   TBranch        *b_det_l;   //!
   TBranch        *b_energy_l;   //!
   TBranch        *b_x_l;   //!
   TBranch        *b_y_l;   //!
   TBranch        *b_z_l;   //!
   TBranch        *b_timeInterval;   //!
   TBranch        *b_t2lastshowermuon;   //!
   TBranch        *b_t2lastshowermuonWithMR;   //!

   LiTree(TTree * /*tree*/ =0) { }
   virtual ~LiTree() { }
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

   ClassDef(LiTree,0);
};

#endif

#ifdef LiTree_cxx
void LiTree::Init(TTree *tree)
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
   fChain->SetBranchAddress("t2lastshowermuon", &t2lastshowermuon, &b_t2lastshowermuon);
   fChain->SetBranchAddress("t2lastshowermuonWithMR", &t2lastshowermuonWithMR, &b_t2lastshowermuonWithMR);
}

Bool_t LiTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef LiTree_cxx
