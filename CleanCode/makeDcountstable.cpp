// makeDcountstable.cpp
// creates (xB,theta,z) tables: xb theta z N Nw Nw2
// from pt2_<SCHEME>_<TARGET>.root produced in the no-phi contingency mode.
//
// directories:   n_evt/   w_pt2/   w_pt2sq/
// Histogram format :  h<SCH>_<TGT>_zK_count   h<SCH>_<TGT>_zK   h<SCH>_<TGT>_zK_sq
//
// Build: g++ -std=c++17 makeDcountstable.cpp -o makeDcountstable $(root-config --cflags --libs)

//g++ -std=c++17 makeDcountstable.cpp -o makeDcountstable $(root-config --cflags --libs)
//
//# Default: DAT, reads pt2_DAT_<TARGET>.root and writes Dcounts4D_DAT_<TARGET>.txt
//./makeDcountstable --target CxC
//
//# BIS:
//./makeDcountstable --scheme BIS --target Sn
//
//# (Optional) custom input path:
//./makeDcountstable --target LD2 --in mydir/pt2_DAT_LD2.root


#include <TFile.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TROOT.h>
#include <TStyle.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

struct Cfg {
  std::string scheme = "DAT";        // DAT (default) or BIS
  std::string target;                // CxC, LD2, Cu, Sn, C1, C2...
  std::string inpath;                // optional override
} cfg;

// Binnings (must match what you used to make the files)
struct Grid { std::vector<double> x, th; };

static const Grid GRID_DAT{
  {0.10, 0.11, 0.13, 0.16, 0.20, 0.25, 0.36, 1.00}, // 7 bins
  {7.6, 8.4, 10.3, 23.0}                            // 3 bins
};
static const Grid GRID_BIS{
  {0.10, 0.11, 0.15, 0.19, 0.29, 1.00},             // 5 bins
  {7.6, 8.8, 11.0, 23.0}                            // 3 bins
};

static void usage(const char* p){
  std::cerr
    << "Usage: " << p << " --target <TARGET> [--scheme DAT|BIS] [--in <pt2_file.root>]\n"
    << "Reads:  pt2_<SCHEME>_<TARGET>.root (if --in is not given)\n"
    << "Writes: Dcounts4D_<SCHEME>_<TARGET>.txt  (columns: xb theta z N Nw Nw2)\n";
}

static void parseArgs(int argc, char** argv){
  for (int i=1;i<argc;++i){
    std::string a = argv[i];
    if (a=="--scheme" && i+1<argc) cfg.scheme = argv[++i];
    else if (a=="--target" && i+1<argc) cfg.target = argv[++i];
    else if (a=="--in" && i+1<argc) cfg.inpath = argv[++i];
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
  if (cfg.inpath.empty()){
    cfg.inpath = "pt2_" + cfg.scheme + "_" + cfg.target + ".root";
  }
}

static TH2* get2D(TDirectory* d, const std::string& name){
  if (!d) return nullptr;
  TH2* h = dynamic_cast<TH2*>(d->Get(name.c_str()));
  return h;
}

int main(int argc, char** argv){
  parseArgs(argc, argv);
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  const Grid& G = (cfg.scheme=="DAT") ? GRID_DAT : GRID_BIS;
  const int Nx  = (int)G.x.size()-1;  // 7 for DAT, 5 for BIS
  const int Nth = (int)G.th.size()-1; // 3
  const int Nz  = 5;                  // z bins

  const std::string tag = cfg.scheme; // "DAT" or "BIS"

  // Open input file
  TFile* fin = TFile::Open(cfg.inpath.c_str(), "READ");
  if (!fin || fin->IsZombie()){
    std::cerr << "[FATAL] Cannot open input: " << cfg.inpath << "\n";
    return 10;
  }

  // Directories
  TDirectory* dN    = dynamic_cast<TDirectory*>(fin->Get("n_evt"));
  TDirectory* dWpt2 = dynamic_cast<TDirectory*>(fin->Get("w_pt2"));
  TDirectory* dWpt4 = dynamic_cast<TDirectory*>(fin->Get("w_pt2sq"));

  if (!dN || !dWpt2 || !dWpt4){
    std::cerr << "[FATAL] Missing expected directories (n_evt, w_pt2, w_pt2sq)\n";
    return 11;
  }

  // Prepare output
  const std::string outTxt = "Dcounts4D_" + cfg.scheme + "_" + cfg.target + ".txt";
  std::ofstream out(outTxt);
  if (!out){
    std::cerr << "[FATAL] Cannot write: " << outTxt << "\n";
    return 12;
  }

  out << "# target=" << cfg.target << " scheme=" << cfg.scheme
      << " file=" << cfg.inpath << "\n";
  out << "# columns: xb  theta  z  N  Nw  Nw2\n";

  // Loop over z bins -> grab the three histos; then loop (xb,theta) bins
  for (int iz=1; iz<=Nz; ++iz){
    // Build names
    char n_count[128], n_w[128], n_w2[128];
    std::snprintf(n_count, sizeof(n_count), "h%s_%s_z%d_count", tag.c_str(), cfg.target.c_str(), iz);
    std::snprintf(n_w,     sizeof(n_w),     "h%s_%s_z%d",       tag.c_str(), cfg.target.c_str(), iz);
    std::snprintf(n_w2,    sizeof(n_w2),    "h%s_%s_z%d_sq",    tag.c_str(), cfg.target.c_str(), iz);

    TH2* hN   = get2D(dN,    n_count);
    TH2* hW   = get2D(dWpt2, n_w);
    TH2* hW2  = get2D(dWpt4, n_w2);

    if (!hN || !hW || !hW2){
      std::cerr << "[WARN] Missing one or more histos for z=" << iz
                << "  ("
                << (hN?"":"N ") << (hW?"":"W ") << (hW2?"":"W2 ")
                << ") — skipping this z bin\n";
      continue;
    }

    // Sanity: bin counts match grid?
    if (hN->GetNbinsX()!=Nx || hN->GetNbinsY()!=Nth){
      std::cerr << "[WARN] Unexpected binning in " << n_count
                << "  got (" << hN->GetNbinsX() << "," << hN->GetNbinsY()
                << ") expected (" << Nx << "," << Nth << ")\n";
    }

    for (int ith=1; ith<=Nth; ++ith){
      for (int ixb=1; ixb<=Nx; ++ixb){
        const double N   = hN ->GetBinContent(ixb, ith);
        const double Nw  = hW ->GetBinContent(ixb, ith);
        const double Nw2 = hW2->GetBinContent(ixb, ith);

        out << ixb << "  " << ith << "  " << iz << "  "
            << N   << "  " << Nw  << "  " << Nw2 << "\n";
      }
    }
  }

  out.close();
  fin->Close();
  std::cout << "[OK] Wrote " << outTxt << "\n";
  return 0;
}
