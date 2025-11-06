// makeDforunfold.cpp 
//----------------------------------
// Create 25 xb–theta 2D histograms per binning (DAT, BIS) from hadron tree,
// organized in 2 directories (one wfor weights in pt2, and one for weightis in pt2*pt2)
// No top level of only_e, we dont need that here. 

// Contingency plan, if there  is not enough phih data, then proceed to N-1 dimensions (only three) 
// by restricting phih bins, i.e., make 5 histograms per binning (DAT, BIS); aka only z 

// Compile: g++ makeDforunfold.cpp -o makeDforunfold $(root-config --cflags --libs)
// how 2 run: 
// Original behavior (with phih directories): 5*25 hists + top-level only-e (DAT & BIS)
// ./makeDforunfold --with-phi
// CONTINGENCY (NO phih binning): 25 hists at top level + top-level only-e 
// ./makeDforunfold --no-phi

//# With custom inputs (all flags are optional)
//./makeDforunfold --no-phi \
  --input build/Trees_CxC_RGD.root \
  --target CxC \
  --tree-had tEv_had_CxC_RGD \
  --tree-ele tEv_onlye_CxC_RGD
//--------------------------------------

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cstdio>
#include <cstdlib>
#include <cmath>

// =================== inputs and whatnot =============================================

// Input ROOT file and TTree name
const std::string INPUT_FILE = "~/dump/Trees_Sn_RGD.root";
const std::string TREE_NAME  = "tEv_had_Sn_RGD";
const std::string TARGET_TAG = "Sn";
// ALSO DO C1 and C2
// Mode: if true → include phi bins (25 histos total)
//       if false → only z bins (5 histos total)
const bool USE_PHI_BINS = false; // toggle manually for now

// =================== Binning definitions=============================================
// 
const double Z_EDGES[6]   = {0.30, 0.38, 0.46, 0.54, 0.62, 0.70};
const double PHI_EDGES[6] = {0.0, 72.0, 144.0, 216.0, 288.0, 360.0};

// DAT and BIS xB–theta grids
struct Grid {
    std::vector<double> xb;
    std::vector<double> th;
};

const Grid GRID_DAT = {
    {0.10, 0.11, 0.13, 0.16, 0.20, 0.25, 0.36, 1.00},
    {7.6, 8.4, 10.3, 23.0}
};

const Grid GRID_BIS = {
    {0.10, 0.11, 0.15, 0.19, 0.29, 1.00},
    {7.6, 8.8, 11.0, 23.0}
};

// ================================================================
//book TH2F with given binning
TH2F* BookHist(const char* name, const Grid& G) {
    TH2F* h = new TH2F(name, name,
                       (int)G.xb.size() - 1, G.xb.data(),
                       (int)G.th.size() - 1, G.th.data());
    h->GetXaxis()->SetTitle("x_{B}");
    h->GetYaxis()->SetTitle("#theta_{e} [deg]");
    h->Sumw2(); // ensure proper weighted errors
    return h;
}

void WriteScheme(TTree* T, const Grid& G, const std::string& tag, bool usePhi) {
    std::string outFile = "pt2_" + tag + "_" + TARGET_TAG + ".root";
    TFile* fOut = new TFile(outFile.c_str(), "RECREATE");

    // Directories: weights and counts
    TDirectory* d_pt2   = fOut->mkdir("w_pt2");
    TDirectory* d_pt2sq = fOut->mkdir("w_pt2sq");
    TDirectory* d_nevt  = fOut->mkdir("n_evt");

    // Tree branches
    double xb, th, z, pt2, phih;
    T->SetBranchAddress("xb",   &xb);
    T->SetBranchAddress("th",   &th);
    T->SetBranchAddress("z",    &z);
    T->SetBranchAddress("pt2",  &pt2);
    T->SetBranchAddress("phih", &phih);

    // Book histograms for each dir
    auto bookVec = [&](TDirectory* d, const char* suffix, bool withPhi) {
        std::vector<TH2F*> v;
        d->cd();
        if (withPhi) {
            for (int iz=0; iz<5; ++iz) {
                for (int iph=0; iph<5; ++iph) {
                    char hname[96];
                    std::snprintf(hname, sizeof(hname),
                                  "h%s_%s_z%d_phi%d%s",
                                  tag.c_str(), TARGET_TAG.c_str(),
                                  iz+1, iph+1, suffix);
                    v.push_back(BookHist(hname, G));
                }
            }
        } else {
            for (int iz=0; iz<5; ++iz) {
                char hname[96];
                std::snprintf(hname, sizeof(hname),
                              "h%s_%s_z%d%s",
                              tag.c_str(), TARGET_TAG.c_str(),
                              iz+1, suffix);
                v.push_back(BookHist(hname, G));
            }
        }
        return v;
    };

    std::vector<TH2F*> H_pt2   = bookVec(d_pt2,   /*suffix*/"",      usePhi);
    std::vector<TH2F*> H_pt2sq = bookVec(d_pt2sq, /*suffix*/"_sq",   usePhi);
    std::vector<TH2F*> H_nevt  = bookVec(d_nevt,  /*suffix*/"_count",usePhi);

    auto indexFor = [&](int iz, int iph)->int {
        return usePhi ? (iz*5 + iph) : iz;
    };

    // Fill
    Long64_t nEntries = T->GetEntries();
    std::cout << "[INFO] Filling " << nEntries << " entries for " << tag << " ..." << std::endl;

    for (Long64_t i=0; i<nEntries; ++i) {
        T->GetEntry(i);

        // z bin
        int iz = -1;
        for (int j=0; j<5; ++j) {
            if (z >= Z_EDGES[j] && z < Z_EDGES[j+1]) { iz = j; break; }
        }
        if (iz == -1) continue;

        if (usePhi) {
            // phi bin (all open on upper edge; adjust last bin to include 360 if needed)
            int iph = -1;
            for (int j=0; j<5; ++j) {
                if (phih >= PHI_EDGES[j] && phih < PHI_EDGES[j+1]) { iph = j; break; }
            }
            if (iph == -1) continue;

            int idx = indexFor(iz, iph);
            H_pt2[idx]->Fill(xb, th, pt2);
            H_pt2sq[idx]->Fill(xb, th, pt2*pt2);
            H_nevt[idx]->Fill(xb, th, 1.0);
        } else {
            int idx = indexFor(iz, 0);
            H_pt2[idx]->Fill(xb, th, pt2);
            H_pt2sq[idx]->Fill(xb, th, pt2*pt2);
            H_nevt[idx]->Fill(xb, th, 1.0);
        }
    }

    // Write out
    d_pt2->cd();   for (auto* h : H_pt2)   h->Write();
    d_pt2sq->cd(); for (auto* h : H_pt2sq) h->Write();
    d_nevt->cd();  for (auto* h : H_nevt)  h->Write();

    fOut->Write();
    fOut->Close();

    std::cout << "[INFO] Wrote file: " << outFile << std::endl;
}

int main() {
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);
    TH1::SetDefaultSumw2(true);

    TFile* fIn = new TFile(INPUT_FILE.c_str(), "READ");
    TTree* T = (TTree*)fIn->Get(TREE_NAME.c_str());
    if (!T) {
        std::cerr << "[ERROR] Could not find tree: " << TREE_NAME << std::endl;
        return 1;
    }

    WriteScheme(T, GRID_DAT, "DAT", USE_PHI_BINS);
    WriteScheme(T, GRID_BIS, "BIS", USE_PHI_BINS);

    std::cout << "[DONE] All output written." << std::endl;
    return 0;
}