{
    gROOT->ProcessLine(".x load.C");
    gStyle->SetOptStat(0);
    gStyle->SetMarkerStyle(24);
    gStyle->SetMarkerSize(0.7);
    

    // // Initialize BetaSpec(string name, int Z, int A, int charge=-1, string dir=".", bool useHuber=false)
    BetaSpec b("B12", 5, 12);
    // BetaSpec b("N12", 7, 12, 1);
    // BetaSpec b("Bi212", 83, 212);
    // BetaSpec b("Bi214", 83, 214);
    // BetaSpec b("Bi210", 83, 210);
    // BetaSpec b("Tl208", 81, 208); b.Qmax = 4.994;
    // BetaSpec b("Kr85", 36, 85);
    // BetaSpec b("Ar39", 18, 39);
    // BetaSpec b("K40", 19, 40);
    // BetaSpec b("C9", 6, 9, 1);
    // BetaSpec b("Li9", 3, 9);
    // BetaSpec b("Li9_beta_only", 3, 9);
    // BetaSpec b("Li9_beta_neutron", 3, 9);
    // BetaSpec b("He8", 2, 8); b.Qmax = 10.652;
    // BetaSpec b("He8_beta_only", 2, 8); b.Qmax = 10.652;
    // BetaSpec b("He8_beta_neutron", 2, 8);
    // BetaSpec b("Li8", 3, 8);
    // BetaSpec b("B8", 5, 8, 1);
    // BetaSpec b("C10", 6, 10, 1); b.Qmax = 2.651;
    // BetaSpec b("C11", 6, 11, 1);
    // BetaSpec b("X117", 46, 117); // hypothetical Z=46, A=117, Q=10
    // BetaSpec b("X117", 46, 117, 1); // hypothetical Z=46, A=117, Q=10, beta+
    
    
    // // Glance at the spectra and the size of the corrections
    b.Glance();
    
    // // Randomly generate event and save to TTree
    // b.MakeTree();
    
    
    // // Or, one can directly access the histograms and gamma/br/etc from the object members:    
        // int nBranch;    // number of branches
        // double Q[MAXBRANCH];    // Q-value of branches, size nBranch
        // double BR[MAXBRANCH];   // branching ratio.
        // int dSpin[MAXBRANCH];   // delta spin
        // int dParity[MAXBRANCH]; // delta parity
        // int nGamma[MAXBRANCH];  // number of correlated gammas
        // double eGamma[MAXBRANCH][MAXGAMMA]; //energy of the correlated gammas, size (nBranch x nGamma[br_index])
        // TH1F* hBeta[MAXBRANCH]; // the corrected beta spectra, size nBranch, binwidth 1keV
        // TH1F* hBetaNoCorrection[MAXBRANCH];  // the un-corrected beta spectra (no Fermi/Screening/etc.)
    
    // // Others Visulization Examples
    // b.Draw("Fermi");    
    // b.Draw("FiniteSizeEM");
    // b.Draw("FiniteSizeWI");
    // b.Draw("Screening");
    // b.Draw("WeakMagnetism");
    // b.Draw("ForbiddennessHuber");
    // b.Draw();    
    // b.CompareSpectra();    
    // b.CompareSize();

}