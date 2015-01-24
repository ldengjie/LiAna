#ifndef BETASPEC_H
#define BETASPEC_H

#include <string>
#include "TString.h"
#include "TComplex.h"

#define MAXBRANCH 200
#define MAXGAMMA  40

class TH1F;
class TGraph;

class BetaSpec {
public:
    BetaSpec();
    BetaSpec(string name, int z, int a, int charge=-1, string dir=".", bool useHuber=false);
    virtual ~BetaSpec();
    
    void SetIsotope(string name, int z, int a, int charge);
    void LoadSpectra(string filename);
    void MakeBetaSpectra();
    void MakeTree(string filename="", int nEvents=1000000);
    
    // correction functions
    double Fermi(int z, double Te, int charge=-1);
    double FiniteSizeEM(int z, double Te, int charge=-1);
    double FiniteSizeWI(int z, double Te, double W0, int charge=-1);
    double Screening(int z, double Te, int charge=-1);
    double WeakMagnetism(double Te, int charge=-1);
    bool UseForbiddenCorrection(int deltaSpin, int deltaParity);
    double Forbiddenness(int z, double Te, double Q0, int deltaSpin, int deltaParity, int charge=-1);
    double ForbiddennessHuber(double Te, double Q0, int deltaSpin);
    
    
    void Draw();
    TGraph* Draw(TString option);
    void CompareSpectra();
    void CompareSize();
    void Glance();
    void Align(TGraph* g);
    
    static TComplex CGamma(TComplex z);
    static int WeightedChoice(double* weights, int size);
    
    string isotope;
    int Z;
    int A;
    int C;  // charge (-1 for electron, 1 for positron)
    int Zdaughter;
    double Qmax;    // maximum energy level difference
    double R;  // nucleus charge radius
    
    double me;
    double alpha;
    TGraph* gScreeningV;
    bool useForbiddennessHuber;
    string baseDir;
    
    int nBranch;    // number of branches
    double Q[MAXBRANCH];    // Q-value of branches, size nBranch
    double BR[MAXBRANCH];   // branching ratio.
    int dSpin[MAXBRANCH];   // delta spin
    int dParity[MAXBRANCH]; // delta parity
    int nGamma[MAXBRANCH];  // number of correlated gammas
    double eGamma[MAXBRANCH][MAXGAMMA]; //energy of the correlated gammas, size (nBranch x nGamma[br_index])
    TH1F* hBeta[MAXBRANCH]; // the corrected beta spectra, size nBranch, binwidth 1keV
    TH1F* hBetaNoCorrection[MAXBRANCH];  // the un-corrected beta spectra (no Fermi/Screening/etc.) 
    

};

#endif
