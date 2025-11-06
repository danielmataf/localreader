//makeRcountsTable.cpp
// ----------------------------------------------
//gets the root files produced by the from the makeRforunfold.cpp
//gets N (nber of events per bin of (xb,theta,pt2,phih,z)) and writes them to text files
//text files are organizzed this way 
//for hadron counts :
//xb  theta  z  pt2  phih  N
//for electron counts only:
//xb  theta  N
//the nbr of bin is an index, from 1 to 5 for hadrons, following the name of each histo, 
//then according to the binning defined in xb theta, we will call those given bins
//
//how 2 compil and run :
//g++ makeRcountsTable.cpp -o makeRcountsTable $(root-config --cflags --libs)
//
//# Default scheme DAT, file inferred: R_unfold_DAT_CxC.root
//./makeRcountsTable --target CxC
//
//# Explicit scheme + file inferred
//./makeRcountsTable --scheme BIS --target Sn
//
// CONITINGENCY PLAN /!\ DROP PHIH
// if there is no phi dirs in the root file, then it will proceed to 4D histos only
 // ./makeRcountsTable --scheme DAT --target CxC --no-phi
//
//# Custom input filename (overrides pattern)
//./makeRcountsTable --in myfile.root --scheme DAT --target LD2
// ----------------------------------------------


#include <TFile.h>
#include <TDirectory.h>
#include <TH2.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>

struct Cfg {
  std::string scheme = "DAT";     // DAT or BIS
  std::string target;             // CxC, Cu, Sn, LD2, C1, C2, ...
  bool no_phi = false;            // contingency flag
} cfg;

static void usage(const char* p){
  std::cerr
    << "Usage: " << p << " --target <TARGET> [--scheme DAT|BIS] [--no-phi]\n"
    << "Reads:  R_unfold_<SCHEME>_<TARGET>.root\n"
    << "Writes: Ecounts_<SCHEME>_<TARGET>.txt                 (xb theta N)\n"
    << "        Hcounts_<SCHEME>_<TARGET>.txt                 (xb theta z pt2 phih N)   [default]\n"
    << "        Hcounts4D_<SCHEME>_<TARGET>.txt               (xb theta z pt2 N)        [--no-phi]\n";
}

static void parseArgs(int argc, char** argv){
  for (int i=1;i<argc;++i){
    std::string a = argv[i];
    if (a=="--scheme" && i+1<argc) cfg.scheme = argv[++i];
    else if (a=="--target" && i+1<argc) cfg.target = argv[++i];
    else if (a=="--no-phi") cfg.no_phi = true;
    else if (a=="-h" || a=="--help"){ usage(argv[0]); std::exit(0); }
    else { std::cerr << "[WARN] Unknown option: " << a << "\n"; usage(argv[0]); std::exit(1); }
  }
  for (auto& c : cfg.scheme) c = std::toupper(c);
  if (cfg.scheme!="DAT" && cfg.scheme!="BIS"){
    std::cerr << "[FATAL] --scheme must be DAT or BIS\n"; std::exit(2);
  }
  if (cfg.target.empty()){
    std::cerr << "[FATAL] --target is required\n"; std::exit(3);
  }
}

static std::string inputRootPath(){
  return "R_unfold_" + cfg.scheme + "_" + cfg.target + ".root";
}

// Electron-only 2D (xB,theta) histogram
static std::string onlyEName(){
  return "h" + cfg.scheme + "_onlye_" + cfg.target;
}

// Hadron histo names:
// With phi_h directories:  phi_K / h<SCHEME>_<TARGET>_phiK_zI_ptJ
// Without phi_h dirs:      top-level h<SCHEME>_<TARGET>_zI_ptJ
static std::string hadronName_withPhi(int iz, int ip, int iph){
  std::ostringstream ss;
  ss << "h" << cfg.scheme << "_" << cfg.target << "_phi" << iph
     << "_z" << iz << "_pt" << ip;
  return ss.str();
}
static std::string hadronName_noPhi(int iz, int ip){
  std::ostringstream ss;
  ss << "h" << cfg.scheme << "_" << cfg.target
     << "_z" << iz << "_pt" << ip;
  return ss.str();
}

static bool hasPhiDir(TFile* f){
  return f && f->GetDirectory("phi_1");
}

int main(int argc, char** argv){
  parseArgs(argc, argv);

  const std::string inroot = inputRootPath();
  TFile* fin = TFile::Open(inroot.c_str(), "READ");
  if (!fin || fin->IsZombie()){
    std::cerr << "[FATAL] Cannot open " << inroot << "\n"; return 10;
  }

  // --- Electron counts (always the same) ---
  const std::string eName = onlyEName();
  TH2F* hE = dynamic_cast<TH2F*>(fin->Get(eName.c_str()));
  if (!hE){
    std::cerr << "[FATAL] Missing electron histo: " << eName << "\n"; return 11;
  }
  const int Nx = hE->GetXaxis()->GetNbins(); // xB bins
  const int Ny = hE->GetYaxis()->GetNbins(); // theta bins

  {
    const std::string outE = "Ecounts_" + cfg.scheme + "_" + cfg.target + ".txt";
    std::ofstream out(outE);
    out << "# xb  theta  N\n";
    for (int ix=1; ix<=Nx; ++ix){
      for (int iy=1; iy<=Ny; ++iy){
        const double N = hE->GetBinContent(ix,iy);
        out << ix << " " << iy << " " << N << "\n";
      }
    }
    out.close();
    std::cout << "[INFO] Wrote " << outE << "\n";
  }

  // --- Hadron counts ---
  const int Nz = 5, Np = 5, Nphi = 5;
  const bool phiDirsPresent = hasPhiDir(fin);

  if (!cfg.no_phi){
    // 5D mode: expect phi_1..phi_5 (if not present, we still try to read top-level names per iph)
    const std::string outH = "Hcounts_" + cfg.scheme + "_" + cfg.target + ".txt";
    std::ofstream out(outH);
    out << "# xb  theta  z  pt2  phih  N\n";

    for (int iph=1; iph<=Nphi; ++iph){
      TDirectory* dphi = fin->GetDirectory(Form("phi_%d", iph));
      // If phi dirs exist, change into dir; else, we’ll try top-level names (defensive)
      for (int ip=1; ip<=Np; ++ip){
        for (int iz=1; iz<=Nz; ++iz){
          TH2F* h2 = nullptr;
          if (dphi){
            h2 = dynamic_cast<TH2F*>(dphi->Get(hadronName_withPhi(iz,ip,iph).c_str()));
          } else {
            // fallback attempt (rare): a file lacking dirs but user still running 5D mode
            h2 = dynamic_cast<TH2F*>(fin->Get(hadronName_withPhi(iz,ip,iph).c_str()));
          }
          for (int ix=1; ix<=Nx; ++ix){
            for (int iy=1; iy<=Ny; ++iy){
              const double N = h2 ? h2->GetBinContent(ix,iy) : 0.0;
              out << ix << " " << iy << " " << iz << " " << ip << " " << iph << " " << N << "\n";
            }
          }
        }
      }
    }
    out.close();
    std::cout << "[INFO] Wrote " << outH << "\n";
  } else {
    // 4D contingency: no phi dimension — just read top-level histos by (z, pt2)
    if (phiDirsPresent){
      std::cout << "[WARN] File has phi_* directories, but --no-phi was requested. "
                   "Proceeding to read top-level histos only (ignoring phi_*).\n";
    }
    const std::string outH = "Hcounts4D_" + cfg.scheme + "_" + cfg.target + ".txt";
    std::ofstream out(outH);
    out << "# xb  theta  z  pt2  N\n";

    for (int ip=1; ip<=Np; ++ip){
      for (int iz=1; iz<=Nz; ++iz){
        TH2F* h2 = dynamic_cast<TH2F*>(fin->Get(hadronName_noPhi(iz,ip).c_str()));
        for (int ix=1; ix<=Nx; ++ix){
          for (int iy=1; iy<=Ny; ++iy){
            const double N = h2 ? h2->GetBinContent(ix,iy) : 0.0;
            out << ix << " " << iy << " " << iz << " " << ip << " " << N << "\n";
          }
        }
      }
    }
    out.close();
    std::cout << "[INFO] Wrote " << outH << " (no phi_h)\n";
  }

  fin->Close();
  return 0;
}
