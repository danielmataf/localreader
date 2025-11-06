// makeRforunfold.cpp 
//----------------------------------
// Create 125 xb–theta 2D histograms per binning (DAT, BIS) from hadron tree,
// organized in 5 directories (phi_1..phi_5). Also add ONE top-level 2D histo
// for electrons-only per binning. Target hardcoded below.

// Contingency plan, if there  is not enough phih data, then proceed to N-1 dimensions
// by restricting phih bins, i.e., make 25 histograms per binning (DAT, BIS); aka only z and pt2

// Compile: g++ makeRforunfold.cpp -o makeRforunfold $(root-config --cflags --libs)

// how 2 run: 
// Original behavior (with phih directories): 5*25 hists + top-level only-e (DAT & BIS)
// ./makeRforunfold --with-phi
// CONTINGENCY (NO phih binning): 25 hists at top level + top-level only-e 
// ./makeRforunfold --no-phi

//# With custom inputs (all flags are optional)
//./makeRforunfold --no-phi \
  --input build/Trees_CxC_RGD.root \
  --target CxC \
  --tree-had tEv_had_CxC_RGD \
  --tree-ele tEv_onlye_CxC_RGD
//--------------------------------------
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TString.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TStyle.h>

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <memory>

// Configuration tutles and files and mor
struct Config {
  std::string inFile  = "~/dump/Trees_C1_RGD.root"; //example "~/dump/Trees_CxC_RGD.root"
  std::string target  = "C1";                     // e.g. CxC, Cu, LD2, Sn...    ALSO DO C1 and C2 
  std::string treeHad = "tEv_had_C1_RGD";          // e.g. tEv_had_CxC_RGD
  std::string treeEle = "tEv_onlye_C1_RGD";    // e.g. tEv_onlye_CxC_RGD
  bool usePhiBins     = true; // --with-phi / --no-phi
} cfg;

//Hadron binning (uniflrm in both )
static const double Z_EDGES[6]   = {0.30, 0.38, 0.46, 0.54, 0.62, 0.70};  // 5 bins
static const double PT2_EDGES[6] = {0.00, 0.24, 0.48, 0.72, 0.96, 1.20};  // 5 bins
static const double PHI_EDGES[6] = {0.0, 72.0, 144.0, 216.0, 288.0, 360.0}; // 5 bins (deg)

// xB–theta grids ( DATA & BIS)
struct Grid { std::vector<double> x; std::vector<double> th; };

Grid gridDAT{
  /*x*/ {0.10, 0.11, 0.13, 0.16, 0.20, 0.25, 0.36, 1.00},
  /*θ*/ {7.6, 8.4, 10.3, 23.0}
};

Grid gridBIS{
  /*x*/ {0.10, 0.11, 0.15, 0.19, 0.29, 1.00},
  /*θ*/ {7.6, 8.8, 11.0, 23.0}
};

// Helpers
std::string joinEdges(const std::vector<double>& v) {
  std::string s;
  char buf[64];
  for (size_t i = 0; i < v.size(); ++i) {
    std::snprintf(buf, sizeof(buf), "%.3g", v[i]);
    s += buf;
    if (i + 1 < v.size()) s += ", ";
  }
  return s;
}

TString makeSelZPtPhi(int iz, int ip, int iph, bool includePhi) {
  // z in [z_i, z_{i+1})
  // pt2 in [pt_i, pt_{i+1})
  // phi in [phi_i, phi_{i+1}) except last: <=
  char buf[512];
  TString sel;

  std::snprintf(buf, sizeof(buf),
    "(z>=%g && z<%g) && (pt2>=%g && pt2<%g)",
    Z_EDGES[iz], Z_EDGES[iz+1], PT2_EDGES[ip], PT2_EDGES[ip+1]);
  sel = buf;

  if (includePhi) {
    if (iph < 4) {
      std::snprintf(buf, sizeof(buf), " && (phih>=%g && phih<%g)",
                    PHI_EDGES[iph], PHI_EDGES[iph+1]);
    } else {
      std::snprintf(buf, sizeof(buf), " && (phih>=%g && phih<=%g)",
                    PHI_EDGES[iph], PHI_EDGES[iph+1]); // include 360
    }
    sel += buf;
  }
  return sel;
}

TH2F* bookHist(const char* hname, const Grid& G) {
  TH2F* h = new TH2F(hname, hname,
                     (int)G.x.size()-1, G.x.data(),
                     (int)G.th.size()-1, G.th.data());
  h->GetXaxis()->SetTitle("x_{B}");
  h->GetYaxis()->SetTitle("#theta_{e} [deg]");
  h->SetDirectory(nullptr);
  return h;
}

void fillAndWriteHists(TTree* Thad, TTree* Tele,
                       const Grid& G,
                       const std::string& schemeTag,  // "DAT" or "BIS"
                       const std::string& outPath,
                       bool usePhiBins,
                       const std::string& targetTag)
{
  std::unique_ptr<TFile> fout(TFile::Open(outPath.c_str(), "RECREATE"));
  if (!fout || fout->IsZombie()) {
    std::cerr << "[FATAL] Cannot create output file: " << outPath << "\n";
    std::exit(1);
  }

  // ---- Hadron histograms ----
  if (usePhiBins) {
    for (int iph = 0; iph < 5; ++iph) {
      // Create phi directory
      char dname[32];
      std::snprintf(dname, sizeof(dname), "phi_%d", iph+1);
      TDirectory* dphi = fout->mkdir(dname);
      dphi->cd();

      for (int ip = 0; ip < 5; ++ip) {
        for (int iz = 0; iz < 5; ++iz) {
          char hname[128];
          std::snprintf(hname, sizeof(hname),
                        "h%s_%s_phi%d_z%d_pt%d",
                        schemeTag.c_str(), targetTag.c_str(),
                        iph+1, iz+1, ip+1);

          TH2F* h = bookHist(hname, G);
          h->SetDirectory(gDirectory);                  
          TString drawExpr = "th:xb>>"; drawExpr += hname;
          TString sel = makeSelZPtPhi(iz, ip, iph, /*includePhi*/true);

          // Fill silently
          Thad->Draw(drawExpr.Data(), sel.Data(), "goff");
          h->Write();
          delete h;
        }
      }
      fout->cd(); // back to file root
    }
  } else {
    // No phih binning: write 25 histograms at the top level
    for (int ip = 0; ip < 5; ++ip) {
      for (int iz = 0; iz < 5; ++iz) {
        char hname[128];
        std::snprintf(hname, sizeof(hname),
                      "h%s_%s_z%d_pt%d",
                      schemeTag.c_str(), targetTag.c_str(),
                      iz+1, ip+1);

        TH2F* h = bookHist(hname, G);
        h->SetDirectory(gDirectory);
        TString drawExpr = "th:xb>>"; drawExpr += hname;
        TString sel = makeSelZPtPhi(iz, ip, /*iph*/0, /*includePhi*/false);

        Thad->Draw(drawExpr.Data(), sel.Data(), "goff");
        h->Write();
        delete h;
      }
    }
  }

  // ---- Only-electron histogram at top level ----
  {
    char hname[128];
    std::snprintf(hname, sizeof(hname), "h%s_onlye_%s",
                  schemeTag.c_str(), targetTag.c_str());
    TH2F* hOnly = bookHist(hname, G);
    hOnly->SetDirectory(gDirectory);  
    TString drawExpr = "th:xb>>"; drawExpr += hname;
    Tele->Draw(drawExpr.Data(), "", "goff");
    hOnly->Write();
    delete hOnly;
  }

  fout->Write();
  fout->Close();
}

// ------------------------------
// CLI
// ------------------------------
void usage(const char* prog) {
  std::cerr <<
    "Usage: " << prog << " [--with-phi | --no-phi]\n"
    "                 [--input FILE.root]\n"
    "                 [--target TAG]\n"
    "                 [--tree-had NAME] [--tree-ele NAME]\n\n"
    "Writes two output files in the current directory:\n"
    "  unfold_DAT_<TARGET>.root  and  unfold_BIS_<TARGET>.root\n"
    "With --with-phi: each has 5 phi_* directories, each containing 25 (z×pt2) TH2F(xB,theta)\n"
    "With --no-phi:   each has 25 TH2F(xB,theta) at top-level (z×pt2), no phi directories\n"
    "Both always include a top-level only-electron TH2F.\n\n"
    "Default input: build/Trees_targer_RGD.root; trees: tEv_had_target_RGD and tEv_onlye_target_RGD.\n";
}

void parseArgs(int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--with-phi") {
      cfg.usePhiBins = true;
    } else if (a == "--no-phi") {
      cfg.usePhiBins = false;
    } else if (a == "--input" && i+1 < argc) {
      cfg.inFile = argv[++i];
    } else if (a == "--target" && i+1 < argc) {
      cfg.target = argv[++i];
    } else if (a == "--tree-had" && i+1 < argc) {
      cfg.treeHad = argv[++i];
    } else if (a == "--tree-ele" && i+1 < argc) {
      cfg.treeEle = argv[++i];
    } else if (a == "-h" || a == "--help") {
      usage(argv[0]);
      std::exit(0);
    } else {
      std::cerr << "[WARN] Unknown option: " << a << "\n";
      usage(argv[0]);
      std::exit(1);
    }
  }
}

// ##################################
int main(int argc, char** argv) {
  parseArgs(argc, argv);

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);

  // Open input
  std::unique_ptr<TFile> fin(TFile::Open(cfg.inFile.c_str(), "READ"));
  if (!fin || fin->IsZombie()) {
    std::cerr << "[FATAL] Cannot open input file: " << cfg.inFile << "\n";
    return 2;
  }

  // Get trees
  TTree* Thad = dynamic_cast<TTree*>(fin->Get(cfg.treeHad.c_str()));
  TTree* Tele = dynamic_cast<TTree*>(fin->Get(cfg.treeEle.c_str()));
  if (!Thad) {
    std::cerr << "[FATAL] Hadron tree not found: " << cfg.treeHad << "\n";
    return 3;
  }
  if (!Tele) {
    std::cerr << "[FATAL] Only-electron tree not found: " << cfg.treeEle << "\n";
    return 4;
  }

  // Compose output names
  std::string outDAT = "R_unfold_DAT_" + cfg.target + ".root";
  std::string outBIS = "R_unfold_BIS_" + cfg.target + ".root";

  // Log summary
  std::cout << "Input: " << cfg.inFile << "\n"
            << " Trees: had='" << cfg.treeHad << "', onlye='" << cfg.treeEle << "'\n"
            << " Target: " << cfg.target << "\n"
            << " Mode:   " << (cfg.usePhiBins ? "WITH phih bins (5 dirs × 25 hists)" :
                                           "NO phih bins (25 hists at top level)") << "\n"
            << " DAT grid: x={" << joinEdges(gridDAT.x)
            << "}, theta={" << joinEdges(gridDAT.th) << "}\n"
            << " BIS grid: x={" << joinEdges(gridBIS.x)
            << "}, theta={" << joinEdges(gridBIS.th) << "}\n"
            << " Outputs: " << outDAT << "  and  " << outBIS << "\n";

  // Write both schemes
  fillAndWriteHists(Thad, Tele, gridDAT, "DAT", outDAT, cfg.usePhiBins, cfg.target);
  fillAndWriteHists(Thad, Tele, gridBIS, "BIS", outBIS, cfg.usePhiBins, cfg.target);

  std::cout << "Done.\n";
  return 0;
}
