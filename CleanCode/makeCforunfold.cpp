// makeCforunfold.cpp 
//----------------------------------
// Create 25 xb–theta 2D histograms per binning (DAT, BIS) from hadron tree,
// organized in 3 directories (one wfor weights in cphih, and one for weightis in cphih*cphih)
// a third one is for normal number of counts (distribution no weights ).
// No top level of only_e, we dont need that here. 

// NO CONTINGENCY PLAN we use phih either way, the analysis is restricted only to 4D  
// How to compile:
// g++ makeCforunfold.cpp -o makeCforunfold $(root-config --cflags --libs)
//
// Default: phih in degrees (converted internally)
// ./makeCforunfold
//
//# If your phih branch is in radians:
//./makeCforunfold --phi-rad
//--------------------------------------

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TStyle.h>
#include <cstdio>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

// edit paths and tags
const std::string INPUT_FILE = "~/dump/Trees_LD2_RGD.root";
const std::string TREE_NAME  = "tEv_had_LD2_RGD";
const std::string TARGET_TAG = "LD2";

// option rad default is  degrees
bool PHI_IS_RAD = false; // set by --phi-rad

// bins
const double Z_EDGES[6]   = {0.30, 0.38, 0.46, 0.54, 0.62, 0.70};
const double PT2_EDGES[6] = {0.00, 0.24, 0.48, 0.72, 0.96, 1.20};

struct Grid { std::vector<double> xb; std::vector<double> th; };

const Grid GRID_DAT = {
  {0.10, 0.11, 0.13, 0.16, 0.20, 0.25, 0.36, 1.00},
  {7.6, 8.4, 10.3, 23.0}
};

const Grid GRID_BIS = {
  {0.10, 0.11, 0.15, 0.19, 0.29, 1.00},
  {7.6, 8.8, 11.0, 23.0}
};

// more
TH2F* BookHist(const char* name, const Grid& G) {
  TH2F* h = new TH2F(name, name,
                     (int)G.xb.size()-1, G.xb.data(),
                     (int)G.th.size()-1, G.th.data());
  h->GetXaxis()->SetTitle("x_{B}");
  h->GetYaxis()->SetTitle("#theta_{e} [deg]");
  h->Sumw2();
  return h;
}

void WriteScheme(TTree* T, const Grid& G, const std::string& tag) {
  const std::string outFile = "cos_" + tag + "_" + TARGET_TAG + ".root";
  TFile* fOut = new TFile(outFile.c_str(), "RECREATE");

  TDirectory* d_dist  = fOut->mkdir("dist");   // counts
  TDirectory* d_wcos  = fOut->mkdir("w_cos");  // sum cos(phi)
  TDirectory* d_wcos2 = fOut->mkdir("w_cos2"); // sum cos^2(phi)

  double xb, th, z, pt2, phih;
  T->SetBranchAddress("xb",   &xb);
  T->SetBranchAddress("th",   &th);
  T->SetBranchAddress("z",    &z);
  T->SetBranchAddress("pt2",  &pt2);
  T->SetBranchAddress("phih", &phih); // degrees by default

  auto bookVec = [&](TDirectory* d, const char* suffix){
    std::vector<TH2F*> v; d->cd();
    for (int iz=0; iz<5; ++iz)
      for (int ip=0; ip<5; ++ip) {
        char hname[96];
        std::snprintf(hname, sizeof(hname), "h%s_%s_z%d_pt%d%s",
                      tag.c_str(), TARGET_TAG.c_str(), iz+1, ip+1, suffix);
        v.push_back(BookHist(hname, G));
      }
    return v;
  };

  std::vector<TH2F*> H_dist  = bookVec(d_dist,  "");
  std::vector<TH2F*> H_wcos  = bookVec(d_wcos,  "_cos");
  std::vector<TH2F*> H_wcos2 = bookVec(d_wcos2, "_cos2");

  auto idx = [](int iz, int ip){ return iz*5 + ip; };

  const double DEG2RAD = M_PI/180.0;
  Long64_t nEntries = T->GetEntries();
  std::cout << "[INFO] Filling " << nEntries << " entries for " << tag
            << "  (phi units: " << (PHI_IS_RAD ? "radians" : "degrees") << ")\n";

  for (Long64_t i=0; i<nEntries; ++i) {
    T->GetEntry(i);

    int iz=-1; for (int j=0;j<5;++j) if (z   >= Z_EDGES[j]   && z   < Z_EDGES[j+1]) { iz=j; break; }
    if (iz==-1) continue;
    int ip=-1; for (int j=0;j<5;++j) if (pt2 >= PT2_EDGES[j] && pt2 < PT2_EDGES[j+1]) { ip=j; break; }
    if (ip==-1) continue;

    double phi_for_cos = PHI_IS_RAD ? phih : (phih * DEG2RAD);
    double c  = std::cos(phi_for_cos);
    double c2 = c*c;

    H_dist [idx(iz,ip)]->Fill(xb, th, 1.0);
    H_wcos [idx(iz,ip)]->Fill(xb, th, c);
    H_wcos2[idx(iz,ip)]->Fill(xb, th, c2);
  }

  d_dist->cd();  for (auto* h: H_dist ) h->Write();
  d_wcos->cd();  for (auto* h: H_wcos ) h->Write();
  d_wcos2->cd(); for (auto* h: H_wcos2) h->Write();

  fOut->Write(); fOut->Close();
  std::cout << "[INFO] Wrote file: " << outFile << "\n";
}

// this for rad option ###############  USELESS
void usage(const char* prog){
  std::cerr << "Usage: " << prog << " [--phi-rad]\n"
            << "  Default: phih in degrees. Use --phi-rad if phih is stored in radians.\n";
}



// ###############################

int main(int argc, char** argv) {
  for (int i=1;i<argc;++i){
    std::string a = argv[i];
    if (a=="--phi-rad") PHI_IS_RAD = true;
    else if (a=="-h" || a=="--help") { usage(argv[0]); return 0; }
    else { std::cerr << "[WARN] Unknown option: " << a << "\n"; usage(argv[0]); return 1; }
  }

  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(true);

  TFile* fIn = new TFile(INPUT_FILE.c_str(), "READ");
  TTree* T = (TTree*)fIn->Get(TREE_NAME.c_str());
  if (!T){ std::cerr << "[ERR] Could not find tree: " << TREE_NAME << "\n"; return 1; }

  WriteScheme(T, GRID_DAT, "DAT");
  WriteScheme(T, GRID_BIS, "BIS");

  std::cout << "[DONE] All output written.\n";
  return 0;
}
