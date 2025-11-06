// makeCcountstable.cpp
// Extract (xb,theta,z,pt2) bin contents from cos_<SCHEME>_<TARGET>.root
// and write a flat table: xb  theta  z  pt2  N  Nw  Nw2
//
// Build:
//   g++ makeCcountstable.cpp -o makeCcountstable $(root-config --cflags --libs)
// Run examples:
//   ./makeCcountstable --target Sn
//   ./makeCcountstable --scheme BIS --target Cu




#include <TFile.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TROOT.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

struct Grid {
  std::vector<double> xb;  // bin edges
  std::vector<double> th;  // bin edges
};

static const Grid GRID_DAT{
  /*x*/ {0.10, 0.11, 0.13, 0.16, 0.20, 0.25, 0.36, 1.00},
  /*θ*/ {7.6,  8.4,  10.3, 23.0}
};
static const Grid GRID_BIS{
  /*x*/ {0.10, 0.11, 0.15, 0.19, 0.29, 1.00},
  /*θ*/ {7.6,  8.8,  11.0, 23.0}
};

struct Cfg {
  std::string scheme = "DAT";   // DAT | BIS
  std::string target;           // CxC, Cu, Sn, LD2, C1, C2, ...
} cfg;

static void usage(const char* p){
  std::cerr
    << "Usage: " << p << " --target <TARGET> [--scheme DAT|BIS]\n"
    << "Reads  : cos_<SCHEME>_<TARGET>.root  (expects dirs: dist, w_cos, w_cos2)\n"
    << "Writes : Ccounts_<SCHEME>_<TARGET>.txt\n";
}

static void parseArgs(int argc, char** argv){
  for (int i=1;i<argc;++i){
    std::string a = argv[i];
    if (a=="--scheme" && i+1<argc) cfg.scheme = argv[++i];
    else if (a=="--target" && i+1<argc) cfg.target = argv[++i];
    else if (a=="-h" || a=="--help"){ usage(argv[0]); std::exit(0); }
    else { std::cerr << "[WARN] Unknown option: " << a << "\n"; usage(argv[0]); std::exit(1); }
  }
  for (auto& c: cfg.scheme) c = std::toupper(c);
  if (cfg.scheme!="DAT" && cfg.scheme!="BIS"){
    std::cerr << "[FATAL] --scheme must be DAT or BIS\n"; std::exit(2);
  }
  if (cfg.target.empty()){
    std::cerr << "[FATAL] --target is required\n"; std::exit(3);
  }
}

static const Grid& pickGrid(const std::string& scheme){
  return (scheme=="DAT") ? GRID_DAT : GRID_BIS;
}

int main(int argc, char** argv){
  parseArgs(argc, argv);
  gROOT->SetBatch(true);

  const Grid& G = pickGrid(cfg.scheme);
  const int Nx  = (int)G.xb.size()-1;   // DAT: 7  | BIS: 5
  const int Nth = (int)G.th.size()-1;   // 3

  const std::string inRoot = "cos_" + cfg.scheme + "_" + cfg.target + ".root";
  TFile* f = TFile::Open(inRoot.c_str(), "READ");
  if (!f || f->IsZombie()){
    std::cerr << "[FATAL] Cannot open input file: " << inRoot << "\n";
    return 10;
  }

  // Directories as created by makeCforunfold.cpp
  TDirectory* dDist  = dynamic_cast<TDirectory*>(f->Get("dist"));
  TDirectory* dWcos  = dynamic_cast<TDirectory*>(f->Get("w_cos"));
  TDirectory* dWcos2 = dynamic_cast<TDirectory*>(f->Get("w_cos2"));
  if (!dDist || !dWcos || !dWcos2){
    std::cerr << "[FATAL] Missing one or more directories (dist, w_cos, w_cos2) in " << inRoot << "\n";
    return 11;
  }

  const std::string outTxt = "Ccounts_" + cfg.scheme + "_" + cfg.target + ".txt";
  std::ofstream out(outTxt);
  if (!out){
    std::cerr << "[FATAL] Cannot create output file: " << outTxt << "\n";
    return 12;
  }

  // Header
  out << "# Source: " << inRoot << "\n"
      << "# Scheme: " << cfg.scheme << "   Target: " << cfg.target << "\n"
      << "# Columns: xb  theta  z  pt2  N  Nw  Nw2\n";

  // Loop over z,pt2 (both 1..5), then over xb/theta bins inside each histogram
  for (int iz=1; iz<=5; ++iz){
    for (int ipt=1; ipt<=5; ++ipt){
      // Hist names as produced by makeCforunfold:
      //   dist  : h<DAT|BIS>_<TARGET>_z%02d_pt%02d
      //   w_cos : h<DAT|BIS>_<TARGET>_z%02d_pt%02d_cos
      //   w_cos2: h<DAT|BIS>_<TARGET>_z%02d_pt%02d_cos2
      // (Note: original code used no zero-padding; keep it simple w/out padding.)
      const std::string base = "h" + cfg.scheme + "_" + cfg.target
                             + "_z" + std::to_string(iz) + "_pt" + std::to_string(ipt);

      TH2* hN   = dynamic_cast<TH2*>( dDist ->Get(base.c_str()) );
      TH2* hNw  = dynamic_cast<TH2*>( dWcos ->Get( (base + "_cos").c_str() ) );
      TH2* hNw2 = dynamic_cast<TH2*>( dWcos2->Get( (base + "_cos2").c_str() ) );

      if (!hN || !hNw || !hNw2){
        std::cerr << "[WARN] Missing one or more histograms for z="<<iz<<" pt2="<<ipt
                  << "  (base='"<<base<<"'). Skipping this pair.\n";
        continue;
      }

      // Sanity: bin numbers should match our x/theta counts
      if (hN->GetNbinsX()!=Nx || hN->GetNbinsY()!=Nth){
        std::cerr << "[WARN] Unexpected binning in '"<<base<<"' "
                  << "(Nx,Ny)=("<<hN->GetNbinsX()<<","<<hN->GetNbinsY()
                  << ") expected ("<<Nx<<","<<Nth<<"). Continuing anyway.\n";
      }

      for (int ix=1; ix<=Nx;  ++ix){
        for (int iy=1; iy<=Nth; ++iy){
          const double N   = hN  ->GetBinContent(ix,iy);
          const double Nw  = hNw ->GetBinContent(ix,iy);
          const double Nw2 = hNw2->GetBinContent(ix,iy);
          // Write coded indices
          out << ix << " " << iy << " " << iz << " " << ipt << " "
              << N  << " " << Nw << " " << Nw2 << "\n";
        }
      }
    }
  }

  out.close();
  f->Close();
  std::cout << "[OK] Wrote " << outTxt << "\n";
  return 0;
}
