// makeUNFtable.cpp
// ----------------------------------------------
// Reads unfolded 2D electron histograms from:
//   /home/matamoros/Downloads/RooUnfold/unf2D_<TARGET>outw4.root
//
// Each file must contain a TH2 named:
//   h_xB_thetael_UNF
//
// It writes out a text file with bin indices and counts
// using the same convention as Ecounts_* from makeRcountsTable.cpp:
//
//   unfEcounts_BIS_<TARGET>.txt
//   # xb_bin  theta_bin  N
//   ix  iy  N
//
// IMPORTANT: h_xB_thetael_UNF is based on the BIS xB–theta binning scheme.
// ----------------------------------------------
//
// Compile:
//   g++ makeUNFtable.cpp -o makeUNFtable $(root-config --cflags --libs)
//
// Examples:
//   ./makeUNFtable --target CxC
//   ./makeUNFtable --target Cu
//   ./makeUNFtable --target LD2
//   ./makeUNFtable --target Sn
//   ./makeUNFtable --target C1
//   ./makeUNFtable --target C2
// ----------------------------------------------

#include <TFile.h>
#include <TH2.h>
#include <TROOT.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cctype>

struct Cfg {
  std::string target;  // CxC, Cu, LD2, Sn, C1, C2, ...
  std::string indir = "/home/matamoros/Downloads/RooUnfold/";
} cfg;

static void usage(const char* prog) {
  std::cerr
    << "Usage: " << prog << " --target <TARGET> [--indir DIR]\n"
    << "\n"
    << "Reads:  <indir>/unf2D_<TARGET>outw4.root\n"
    << "  (expecting histogram h_xB_thetael_UNF)\n"
    << "Writes: unfEcounts_BIS_<TARGET>.txt   # xb_bin theta_bin N\n"
    << "\n"
    << "Note: the UNF histogram uses the BIS xB-theta binning scheme.\n";
}

static void parseArgs(int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--target" && i+1 < argc) {
      cfg.target = argv[++i];
    } else if (a == "--indir" && i+1 < argc) {
      cfg.indir = argv[++i];
      // Ensure indir ends with a '/'
      if (!cfg.indir.empty() && cfg.indir.back() != '/')
        cfg.indir += "/";
    } else if (a == "-h" || a == "--help") {
      usage(argv[0]);
      std::exit(0);
    } else {
      std::cerr << "[WARN] Unknown option: " << a << "\n";
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (cfg.target.empty()) {
    std::cerr << "[FATAL] --target is required\n";
    std::exit(2);
  }
}

int main(int argc, char** argv) {
  parseArgs(argc, argv);

  gROOT->SetBatch(true);

  // Build input ROOT filename
  const std::string inroot = cfg.indir + "unf2D_" + cfg.target + "outw4.root";
  TFile* fin = TFile::Open(inroot.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "[FATAL] Cannot open input file: " << inroot << "\n";
    return 10;
  }

  // Get the UNF 2D histogram (xB, theta)
  const std::string hname = "h_xB_thetael_UNF";
  TH2* h2 = dynamic_cast<TH2*>(fin->Get(hname.c_str()));
  if (!h2) {
    std::cerr << "[FATAL] Missing histogram: " << hname
              << " in file " << inroot << "\n";
    return 11;
  }

  const int Nx = h2->GetXaxis()->GetNbins(); // xB bins (BIS)
  const int Ny = h2->GetYaxis()->GetNbins(); // theta bins (BIS)

  // Output text file: unfEcounts_BIS_<TARGET>.txt
  const std::string outName = "unfEcounts_BIS_" + cfg.target + ".txt";
  std::ofstream out(outName);
  if (!out.is_open()) {
    std::cerr << "[FATAL] Cannot create output file: " << outName << "\n";
    return 12;
  }

  out << "# xb_bin  theta_bin  N   (BIS scheme)\n";
  for (int ix = 1; ix <= Nx; ++ix) {
    for (int iy = 1; iy <= Ny; ++iy) {
      const double N = h2->GetBinContent(ix, iy);
      out << ix << " " << iy << " " << N << "\n";
    }
  }
  out.close();

  std::cout << "[INFO] Input : " << inroot  << "\n"
            << "[INFO] Output: " << outName << "\n"
            << "[INFO] Bins : Nx = " << Nx << ", Ny = " << Ny
            << " (BIS scheme)\n";

  fin->Close();
  return 0;
}
