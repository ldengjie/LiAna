#include "BetaSpec.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TDirectory.h"
#include "TColor.h"
#include "TMath.h"
#include "TGraph.h"
#include "TComplex.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"

BetaSpec::BetaSpec()
{
    SetIsotope("Unknown", 0, 0, 0);
}

BetaSpec::BetaSpec(string name, int z, int a, int charge, string dir, bool useHuber)
{
    SetIsotope(name, z, a, charge);
    useForbiddennessHuber = useHuber;
    
    baseDir = dir;
    string filename = dir + "/data/beta/";
    if (charge == 1) filename = dir + "/data/EC/";
    filename += isotope + ".txt";
    LoadSpectra(filename);
}

BetaSpec::~BetaSpec()
{
    for (int i=0; i<nBranch; i++) {
        delete hBeta[i];
        delete hBetaNoCorrection[i];
    }
    if (gScreeningV) delete gScreeningV;
    
}

void BetaSpec::SetIsotope(string name, int z, int a, int charge) 
{   isotope=name; 
    Z=z; 
    A=a; 
    nBranch=0; 
    C = charge; 
    Zdaughter = Z-charge; 
    R = 0.0029*pow(A, 1./3) + 0.0063*pow(A, -1./3) - 0.0017/A; // radius in me
    // R = 0.42587*alpha * pow(A,1./3.);

    me = 0.511; // mass of electron
    alpha = 1/137.035999074; // fine structure constant
    
    // tabulated screening potential
    double x[10] = {  1,    8,    13,    16,   23,    27,    29,    49,    84,    92};
    double y[10] = {1.0, 1.42, 1.484, 1.497, 1.52, 1.544, 1.561, 1.637, 1.838, 1.907};
    gScreeningV = new TGraph(10, x, y);
    
}
    
void BetaSpec::LoadSpectra(string filename)
{
    ifstream f(filename.c_str());
    if (!f) {
        cout << filename << " not found. exiting ..." << endl;
        exit(1);
    }
    cout << "loading " << isotope << " spectra info from " << filename << endl;
    cout << "Q\tB.R.\tdeltaSpin\tdeltaParity\tnGamma\teGammas\n";
    cout << "=================================================" << endl;
    string line;
    while (f.good()) {
        getline(f, line);
        if (line.size() > 2) { nBranch++; };
    }
    
    f.clear();
    f.seekg(0, f.beg);    
    for (int i=0; i<nBranch; i++){
        getline(f, line);
        // cout << line << endl;
        
        stringstream input(line);
        input >> Q[i] >> BR[i] >> dSpin[i] >> dParity[i] >> nGamma[i] ;
        cout << Q[i] << "\t " << BR[i] << "\t " 
            << dSpin[i] << "\t " << dParity[i] << "\t " 
            << nGamma[i] << "\t " ;
        
        for (int j=0; j<nGamma[i]; j++) {
            input >> eGamma[i][j];
            cout << eGamma[i][j] << "\t";
        }
        cout << endl;
    }
    f.close();
    cout << "=================================================" << endl;
    
    Qmax = Q[0]; // by default
    MakeBetaSpectra();
}

void BetaSpec::MakeBetaSpectra()
{
    double binWidth = 1e-3; // MeV
    
    for (int i=0; i<nBranch; i++) {
        double Q0 = Q[i];
        TString name = Form("%s_%03i", isotope.c_str(), i);
        TString title = Form("%s_%.1f_MeV", isotope.c_str(), Q0);
        int nBins = int(floor(Q0/binWidth)) + 1;
        hBeta[i] = new TH1F(name.Data(), title.Data(), nBins, 0, nBins*binWidth);
        hBetaNoCorrection[i] = new TH1F((name+"_nc").Data(), title.Data(), nBins, 0, Q0);
        bool useForbiddenCorrection = UseForbiddenCorrection(dSpin[i], dParity[i]);
        for (int j=1; j<nBins; j++) {
            double Te = (j-0.5)*binWidth;
            double Ee = Te + me;
            double dNdE = sqrt(Ee*Ee - me*me) * (Q0-Te)*(Q0-Te) * Ee;
            hBetaNoCorrection[i]->SetBinContent(j, dNdE);
            
            dNdE *= Fermi(Zdaughter, Te, C)
                  * Screening(Zdaughter, Te, C);
            
            if (useForbiddenCorrection) {
                dNdE *= Forbiddenness(Zdaughter, Te, Q0, dSpin[i], dParity[i], C);
            }            
            if (!(dParity[i]==1 && (dSpin[i]==0 || dSpin[i]==1))) {  // not first non-unique
                dNdE *= FiniteSizeEM(Zdaughter, Te, C)
                      * FiniteSizeWI(Zdaughter, Te, Q0, C)
                      * WeakMagnetism(Te, C);
            }

            // cout << Fermi(Z, Te) << endl;
            hBeta[i]->SetBinContent(j, dNdE);
        }
        hBeta[i]->Scale(1./hBeta[i]->Integral());
        hBetaNoCorrection[i]->Scale(1./hBetaNoCorrection[i]->Integral());
    }
}


double BetaSpec::Fermi(int z, double Te, int charge)
{
    double Ee = Te + me;
    double pe = sqrt(Ee*Ee-me*me);
    
    double gamma0 = sqrt(1- pow(alpha*z, 2));
    double nu = -charge * alpha * z * Ee / pe;
    double result = 2 * (1+gamma0) * exp(TMath::Pi()*nu)
                  * pow(2*pe*R/me, 2*(gamma0-1)) 
                  / pow(TMath::Gamma(2*gamma0+1), 2);;
    TComplex a(gamma0, nu);
    TComplex b = CGamma(a);
    result *= b.Rho2();
    return result;
}

double BetaSpec::FiniteSizeEM(int z, double Te, int charge)
{
    double W = Te/me + 1;
    double gamma0 = sqrt(1- pow(alpha*z, 2));
    
    double result = 1 
        + 13. * pow(alpha*z, 2) / 60.
        + charge * W*R*alpha*z*(41-26*gamma0)/15./(2*gamma0-1) 
        + charge * alpha*z*R*gamma0*(17-2*gamma0)/30./W/(2*gamma0-1);
    return result;
}

double BetaSpec::FiniteSizeWI(int z, double Te, double W0, int charge)
{
    double W = Te/me + 1;
    W0 /= me;
    
    double C0, C1, C2;
    C0 = -233/630.*pow(alpha*z, 2) - pow(W0*R, 2)/5. + charge * 2./35.*W0*R*alpha*z;
    C1 = charge * 21./35.*R*alpha*z + 4./9.*W0*R*R;
    C2 = -4./9.*R*R;
    
    double result = 1 + C0 + C1*W + C2*W*W;
    return result;
}

double BetaSpec::Screening(int z, double Te, int charge)
{
    double gamma0 = sqrt(1- pow(alpha*z, 2));
    double W = Te/me + 1;
    double p = sqrt(W*W-1);
    double y = - charge * alpha * z * W / p;
    
    double Ztilt = z - 1;
    double V0 = - charge * alpha * alpha * pow(Ztilt,4./3.) * gScreeningV->Eval(Ztilt);
    if (W<V0) return 1;
    
    double Wtilt = W - V0 > 1 ? W-V0 : 1.00001;
    double ptilt = sqrt(Wtilt*Wtilt-1);
    double ytilt = - charge * alpha * z * Wtilt / ptilt;
    
    TComplex a(gamma0, ytilt);
    TComplex b = CGamma(a);
    TComplex c(gamma0, y);
    TComplex d = CGamma(c);
      
    double result = Wtilt/W 
        * pow(ptilt/p, 2*gamma0-1) 
        * exp(TMath::Pi() * (ytilt-y))
        * b.Rho2() / d.Rho2();
    
    return result;
    
}

double BetaSpec::WeakMagnetism(double Te, int charge)
{
    return 1 -charge * 0.0048 * (Te+me);  // use a universal WM slope
}


bool BetaSpec::UseForbiddenCorrection(int deltaSpin, int deltaParity)
{
    bool doCorrection = false;
    
    if ((deltaSpin == 0 || deltaSpin == 1) && deltaParity == 0) { // allowed decay
        doCorrection = false;      
    } else if ((deltaSpin == 0 || deltaSpin == 1) && deltaParity == 1) { // Non-unique First Forbidden Decay
        doCorrection = false;                                                                                                                       
    } else if (deltaSpin == 2 && deltaParity == 1) { // Unique First Forbidden Decay
        doCorrection = true;
    } else if (deltaSpin == 3 && deltaParity == 0) { // Unique Second Forbidden Decay
        doCorrection = true;
    } else if (deltaSpin == 4 && deltaParity == 1) { // Unique Third Forbidden Decay
        doCorrection = true;
    } else if (deltaSpin == 5 && deltaParity == 0) { // Unique Fourth Forbidden Decay
        doCorrection = true;
    } else if (deltaSpin == 2 && deltaParity == 0) { // Non-unique Second Forbidden Decay
        doCorrection = true;
    } else if (deltaSpin == 3 && deltaParity == 1) { // Non-unique Third Forbidden Decay
        doCorrection = false; //crude approximation
    } else if (deltaSpin == 4 && deltaParity == 0) { // Non-unique Fourth Forbidden Decay
        doCorrection = false; //crude approximation
    } else {  // higher orders
        doCorrection = false; //crude approximation allowed
    }
    
    if (doCorrection) {
        cout << "   Forbidden Shape Correction: (" << deltaSpin << ", " << deltaParity << ")" << endl; 
    }
    else {
        cout << "No Forbidden Shape Correction: (" << deltaSpin << ", " << deltaParity << ")" << endl;
    }
    return doCorrection;
}

double BetaSpec::Forbiddenness(int z, double Te, double Q0, int deltaSpin, int deltaParity, int charge)
{    
    double tmp = deltaParity; tmp = 0; // currently don't use deltaParity in this fuction
    // shape correction to fobidden decays
    double Ee = Te + me;
    double pe = sqrt(Ee*Ee-me*me);
    double gamma0 = sqrt(1- pow(alpha*z, 2));
    double nu = -charge * alpha * z * Ee / pe;

    double result = 1.0;
    double corr = 0;
    
    /* Get the factorials in the normalization of the forbiddenness correction factor */
    double factorial_1 = 1;
    double dubl_fact_1 = 1;
    for(int k=1; k<=(deltaSpin-1); k++) factorial_1 = factorial_1*k;
    for(int k=1; k<=(2*deltaSpin-1); k++) {dubl_fact_1 = dubl_fact_1*k; k++;} 
    
    for (int j=1; j<=deltaSpin; j++) {
        /* Get Double Factorials over spin */
        double dubl_fact_num = 1;
        for(int k=1; k<=j; k++) { dubl_fact_num = dubl_fact_num*(2*k-1); }
        double dubl_fact_denom = 1;
        for(int k=1; k<=(2*deltaSpin-2*j+1); k++){ dubl_fact_denom = dubl_fact_denom*(k); k++; }

        /* Get Single Factorials over spin */
        double factorial_denom_1 = 1;
        for(int k=1; k<=(j-1); k++){ factorial_denom_1=factorial_denom_1*(k);}
        double factorial_denom_2=1;
        for(int k=1; k<=(deltaSpin-j); k++) { factorial_denom_2=factorial_denom_2*(k); }

        double gamma1 = sqrt(j*j - alpha*z*alpha*z);

        TComplex a(gamma0, nu);
        TComplex b = CGamma(a);
        TComplex c(gamma1, nu);
        TComplex d = CGamma(c);

        corr += 25.1327 
              * factorial_1/dubl_fact_1
              * pow(R, 2.*((double)(deltaSpin-1)))/ (1.+gamma0)
              * dubl_fact_num / (dubl_fact_denom*factorial_denom_1*factorial_denom_2) 
              * pow((Q0-Te)/me, 2.*(deltaSpin-(2*j-1)/2.)-1.)
              * 0.5 * pow(2.*pe/me, 2.*(gamma1-gamma0)) 
              * pow(R, 2.*(gamma1-gamma0)-(2.*j-1.))
              * j * (j+gamma1)
              * pow(TMath::Gamma(2*gamma0+1)/TMath::Gamma(2*gamma1+1), 2)
              * d.Rho2()/b.Rho2();
    }
    result *= corr;
    return result;
}

double BetaSpec::ForbiddennessHuber(double Te, double Q0, int deltaSpin)
{
    double Ee = Te + me;
    double pe = sqrt(Ee*Ee-me*me);
    double pnu = Q0 - Te;
    
    double result = 1.0;

    if (deltaSpin == 0 || deltaSpin==1) {
        result = 1.0;
    } else if (deltaSpin == 2) {
        result = pow(pnu, 2) + pe*pe;
    } else if (deltaSpin == 3) {
        result = pow(pnu, 4) + 10./3.*pow(pnu, 2)*pow(pe, 2) + pow(pe, 4);
    } else if (deltaSpin == 4) {
        result = pow(pnu,6) + 7*pow(pnu, 4)*pow(pe, 2) + 7*pow(pnu, 2)*pow(pe, 4) + pow(pe, 6);
    }

    return result;
}

// double BetaSpec::Fermi(int z, double Te, int charge)
// {
//     double me = 0.511; // mass of electron
//     double E = Te + me;
//     double p = sqrt(E*E - me*me); // momentum
//     double alpha = 1/137.035999074;
//     double rho = 0.0029*pow(A, 1./3) + 0.0063*pow(A, -1./3) - 0.0017/A;
//     double S = sqrt(1-alpha*alpha*z*z);
//     double eta = charge *  alpha * z * E / p;
//     
//     double phi = atan(S/eta);
//     double F=1;
//     F *= 4 * (1+S) / pow(TMath::Gamma(2*S+1),2);
//     F *= pow(2*p*rho, 2*S-2);
//     F *= pow(S*S + eta*eta, S-0.5);
//     F *= exp(2*phi*eta - 2*S + S/6./(S*S+eta*eta));
//     F *= TMath::Pi();
//     return F;   
//     // return 1 - charge * TMath::Pi() * alpha * z * E / p;
// }

void BetaSpec::MakeTree(string filename, int nEvents)
{
    const int N = int(nEvents);
    if (filename == "") filename = isotope + "tree.root";
    TFile f(filename.c_str(), "recreate");
    
    TTree t("T","T");
    double eElectron_; 
    int nGamma_; 
    double eGamma_[MAXGAMMA];
    
    t.Branch("eElectron", &eElectron_, "eElectron/D");
    t.Branch("nGamma", &nGamma_, "nGamma/I");
    t.Branch("eGamma", &eGamma_, "eGamma[nGamma]/D");
    
    for (int i=0; i<N; i++) {
        // eElectron_ = 0;
        // nGamma_ = 0;
        // for (int j=0; j<MAXGAMMA; j++) eGamma[j] = 0;
        int br_i = WeightedChoice(BR, nBranch);
        eElectron_ = hBeta[br_i]->GetRandom();
        nGamma_ = nGamma[br_i];
        for (int j=0; j<nGamma_; j++) {
            eGamma_[j] = eGamma[br_i][j];
        }
        t.Fill();
    }
    t.Write();
    f.Close();
    
}

// Returns a random index according to the weight array
int BetaSpec::WeightedChoice(double* weights, int size)
{
    double totals[MAXBRANCH];
    double running_total = 0;
    
    for (int i=0; i<size; i++) {
        running_total += weights[i];
        totals[i] = running_total;
    }
    double rnd = gRandom->Uniform() * running_total;
    for (int i=0; i<size; i++) {
        if (rnd < totals[i]) {
            return i;
        }
    }
    return size-1;
}

void BetaSpec::Draw()
{
    double binWidth = 1e-3; // MeV
    int maxBins = int(floor(Qmax/binWidth))+1;
    TH1F *hTotal = new TH1F("h", (isotope + " Spectra").c_str(), maxBins, 0, maxBins*binWidth);
    
    TH1F *h;
    for (int i=0; i<nBranch; i++) {
        TString name = Form("h%03i", i);
        h = new TH1F(name.Data(), name.Data(), maxBins, 0, maxBins*binWidth);
        double eGammaTotal = 0;
        for (int j=0; j<nGamma[i]; j++) {
            eGammaTotal += eGamma[i][j];
        }
        // cout << eGammaTotal << endl;
        int nBins = hBeta[i]->GetNbinsX();
        for (int bin=1; bin<nBins; bin++) {
            double content = BR[i] * hBeta[i]->GetBinContent(bin);;
            double e = (bin-0.5)*binWidth + eGammaTotal;
            int currentBin = hTotal->FindBin(e);
            double currentContent = hTotal->GetBinContent(currentBin);
            // cout << bin << " " << content << " " << currentContent << endl;
            h->SetBinContent(currentBin, content);  
            hTotal->SetBinContent(currentBin, currentContent + content);
                
        }
    }
    hTotal->SetLineWidth(2);
    hTotal->SetLineColor(kBlack);
    hTotal->GetXaxis()->SetTitle("E (MeV)");
    hTotal->Draw();
    for (int i=0; i<nBranch; i++) {
        TString name = Form("h%03i", i);
        h = (TH1F*)gDirectory->Get(name.Data());
        h->SetLineWidth(1);
        h->SetLineColor(kBlue);
        // h->SetLineStyle(3);
        h->Draw("same");
    }
}

TGraph* BetaSpec::Draw(TString option)
{
    const int N = 1000;
    double x[N], y[N];
    for (int i=0; i<N; i++) {
        x[i] = 0.01 * (i+5); // energy
        if (option == "Fermi") {
            y[i] = Fermi(Zdaughter, x[i], C);
        }
        else if (option == "FiniteSizeEM") {
            y[i] = FiniteSizeEM(Zdaughter, x[i], C);
        }
        else if (option == "FiniteSizeWI") {
            y[i] = FiniteSizeWI(Zdaughter, x[i], Q[0], C); // use 1st branch as a demo.
        }
        else if (option == "Screening") {
            y[i] = Screening(Zdaughter, x[i], C);
        }
        else if (option == "WeakMagnetism") {
            y[i] = WeakMagnetism(x[i], C); 
        }
        else if (option == "Forbiddenness") {
            y[i] = Forbiddenness(Zdaughter, x[i], Q[0], dSpin[0], dParity[0], C); // use 1st branch as a demo.
        }
        else if (option == "ForbiddennessHuber") {
            y[i] = ForbiddennessHuber(x[i], Q[0], dSpin[0]);
        }
        else {
            cout << "option " << option << " unknown" << endl;
            return 0;
        }
        // cout << x[i] << " " << y[i] << endl;
    }
    TGraph *g = new TGraph(N, x, y);
    Align(g);
    g->Draw("AL");
    g->SetTitle((option + " Correction").Data());
    g->GetXaxis()->SetTitle("E (MeV)");
    return g;
    
}

void BetaSpec::CompareSpectra()
{
    hBetaNoCorrection[0]->SetTitle("Demo of Correction (1st Branch)");
    hBetaNoCorrection[0]->Draw();
    hBeta[0]->Draw("same");
    hBeta[0]->SetLineColor(kRed);
    
    TLegend *leg = new TLegend(0.63, 0.74, 0.89, 0.89);
    leg->AddEntry(hBetaNoCorrection[0],  " Not Corrected", "l");
    leg->AddEntry(hBeta[0],  " Corrected", "l");
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0);
    leg->Draw();
}

void BetaSpec::CompareSize()
{
    TGraph *g1 = Draw("Screening");     Align(g1);
    TGraph *g2 = Draw("FiniteSizeEM");  Align(g2);
    TGraph *g3 = Draw("FiniteSizeWI");  Align(g3);
    TGraph *g4 = Draw("WeakMagnetism"); Align(g4);
    
    TH2F *hc = new TH2F("hc", "S/L/C/WM Corrections", 100, 0, 10, 100, 0.9, 1.1);
    hc->Draw();
    hc->GetXaxis()->SetTitle("E [MeV]");
    hc->GetYaxis()->SetTitle("Correction");
    // g1->SetLineColor(kBlack);   g1->Draw("Lsame");
    g1->SetLineColor(kBlack);   g1->Draw("Lsame");
    g2->SetLineColor(kRed);     g2->Draw("Lsame");
    g3->SetLineColor(kBlue);    g3->Draw("Lsame");
    g4->SetLineColor(kMagenta); g4->Draw("Lsame");
    
    TLegend *leg = new TLegend(0.55, 0.67, 0.90, 0.90);
    leg->AddEntry(g1,  " Screening"     , "l");
    leg->AddEntry(g2,  " Finite Size EM", "l");
    leg->AddEntry(g3,  " Finite Size WI", "l");
    leg->AddEntry(g4,  " Weak Magnetism", "l");
    leg->SetFillColor(kWhite);
    leg->Draw();
}

void BetaSpec::Glance()
{
    TCanvas *c1 = new TCanvas("cBeta", "Glance At Beta Spectra", 1000, 750);
    c1->Divide(2, 2);
    
    c1->cd(1);
    Draw();
    
    c1->cd(2);
    CompareSpectra();
    
    c1->cd(3);
    Draw("Fermi");
    
    c1->cd(4);
    CompareSize();
}

void BetaSpec::Align(TGraph *g)
{
    int nPoints = g->GetN();
    double middle = g->GetY()[nPoints/2-1];
    for (int i=0; i<nPoints; i++) {
        g->GetY()[i] /= middle;
    }
}

// Complex Gamma function
TComplex BetaSpec::CGamma(TComplex z)
{
  // complex gamma function
  // modified f2c of CERNLIB cgamma64.F
  // taken from Knat and converted to compile without ROOT libs

  Double_t zr=z.Re();

  if (z.Im()==0.0) {
    if (zr==0.0) {
      TComplex w(0.0,0.0);
      return w;
    } else if (-fabs(zr)==floor(zr)) {
      TComplex w(0.0,0.0);
      return w;
    }
  }

  TComplex u;
  TComplex v;

  if (zr>=1.0) {
    u=1.0;
    v=z;
  } else if (zr>=0.0) {
    u=1.0/z;
    v=1.0+z;
  } else {
    u=1.0;
    v=1.0-z;
  }

  const Double_t c1=2.50662827463100050;
  const Double_t c2[16]={
     41.624436916439068,-51.224241022374774,
     11.338755813488977, -0.747732687772388,
      0.008782877493061, -0.000001899030264,
      0.000000001946335, -0.000000000199345,
      0.000000000008433,  0.000000000001486,
     -0.000000000000806,  0.000000000000293,
     -0.000000000000102,  0.000000000000037,
     -0.000000000000014,  0.000000000000006
  };

  TComplex w(1.0,0.0);
  TComplex s(c2[0],0.0);
  TComplex cpi(3.14159265358979312,0.0);

  for (int i=1;i<16;++i) {
    w*=(v-((Double_t)i))/(v+((Double_t)(i-1)));
    s+=c2[i]*w;
  }

  w=v+4.5;
  TComplex logw = TComplex::Log(w);
  w=c1*s*TComplex::Exp((v-0.5)*TComplex::Log(w)-w);
  if (zr<0.0) w=cpi/(TComplex::Sin(cpi*z)*w);

  return w*u;
}
