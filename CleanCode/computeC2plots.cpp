// computeC2plots.cpp
//
// Cos(2phi_h) ratio plots.
// - Reads C2counts_<SCHEME>_<TARGET>.txt (columns: xb th z pt N Nw Nw2)
// - Produces 5 PDFs (pt=1..5) for A/LD2 with targets CxC, Cu, Sn
// - Produces 5 PDFs (pt=1..5) for C1/C2 self-ratio
//
// Compile:
//   g++ computeC2plots.cpp -o computeC2plots $(root-config --cflags --libs)
// Run:
//                  # scheme=DAT (default)
//   ./computeC2plots
//   ./computeC2plots --scheme BIS
//
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TH1.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TColor.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

// ---------------- CLI cfg ----------------
struct Cfg {
  std::string scheme = "DAT"; // DAT or BIS
} cfg;

void usage(const char* p){
  std::cerr << "Usage: " << p << " [--scheme DAT|BIS]\n";
}

void parseArgs(int argc, char** argv){
  for (int i=1;i<argc;++i){
    std::string a = argv[i];
    if (a=="--scheme" && i+1<argc) cfg.scheme = argv[++i];
    else if (a=="-h" || a=="--help"){ usage(argv[0]); std::exit(0); }
    else { std::cerr << "[WARN] Unknown option: " << a << "\n"; usage(argv[0]); std::exit(1); }
  }
  for (auto& c: cfg.scheme) c = std::toupper(c);
  if (cfg.scheme!="DAT" && cfg.scheme!="BIS"){
    std::cerr << "[FATAL] --scheme must be DAT or BIS\n"; std::exit(2);
  }
}

// ------------- physics binning -------------
static const double Z_EDGES [6] = {0.30, 0.38, 0.46, 0.54, 0.62, 0.70};
static const double PT2_EDGES[6] = {0.00, 0.24, 0.48, 0.72, 0.96, 1.20};

struct Grid { std::vector<double> x, th; };
static const Grid gridDAT{
  {0.10, 0.11, 0.13, 0.16, 0.20, 0.25, 0.36, 1.00}, // 7 bins
  {7.6, 8.4, 10.3, 23.0}                            // 3 bins
};
static const Grid gridBIS{
  {0.10, 0.11, 0.15, 0.19, 0.29, 1.00}, // 5 bins (not used here for X count, but kept for completeness)
  {7.6, 8.8, 11.0, 23.0}
};

static inline double mid(const double* e, int i){ return 0.5*(e[i-1]+e[i]); } // i=1..5
static inline double mid(const std::vector<double>& v, int i){ return 0.5*(v[i-1]+v[i]); }

// ------------- tuple key hash for maps -------------
// key = (xb, th, z, pt)
using Key = std::tuple<int,int,int,int>;

struct KeyHash {
  std::size_t operator()(const Key& k) const noexcept {
    auto [a,b,c,d] = k;
    std::size_t h1 = std::hash<int>{}(a);
    std::size_t h2 = std::hash<int>{}(b);
    std::size_t h3 = std::hash<int>{}(c);
    std::size_t h4 = std::hash<int>{}(d);
    // simple hash combine
    std::size_t h = h1;
    h ^= (h2 + 0x9e3779b9 + (h<<6) + (h>>2));
    h ^= (h3 + 0x9e3779b9 + (h<<6) + (h>>2));
    h ^= (h4 + 0x9e3779b9 + (h<<6) + (h>>2));
    return h;
  }
};
struct KeyEq {
  bool operator()(const Key& a, const Key& b) const noexcept { return a==b; }
};

// ------------- row content -------------
struct Row {
  double N=0.0;   // counts
  double Nw=0.0;  // sum cos
  double Nw2=0.0; // sum cos^2
};

// ------------- reader -------------
bool readCcounts(const std::string& path, std::unordered_map<Key,Row,KeyHash,KeyEq>& M){
  std::ifstream f(path);
  if(!f){ std::cerr << "[ERR] Can't open " << path << "\n"; return false; }
  std::string line;
  int xb, th, z, pt;
  double N, Nw, Nw2;
  long ln=0, ok=0;
  while (std::getline(f,line)){
    ++ln;
    if (line.empty() || line[0]=='#') continue;
    std::istringstream iss(line);
    if (!(iss>>xb>>th>>z>>pt>>N>>Nw>>Nw2)) continue;
    M[Key{xb,th,z,pt}] = Row{N,Nw,Nw2};
    ++ok;
  }
  std::cout << "[OK] " << path << " rows=" << ok << "\n";
  return true;
}

// ------------- mean cos and error per target -------------
struct MeanC { double c=0.0, ec=0.0; bool ok=false; };

static MeanC meanCos_from_Row(const Row& r){
  MeanC m;
  if (r.N <= 0) return m;
  const double mean = r.Nw / r.N;                 // <cos>
  const double mean2 = r.Nw2 / r.N;               // <cos^2>
  const double var = std::max(0.0, mean2 - mean*mean);
  const double err_mean = std::sqrt(var / r.N);   // σ(<cos>)
  if (std::isfinite(mean) && std::isfinite(err_mean) && err_mean>=0){
    m.c = mean; m.ec = err_mean; m.ok = true;
  }
  return m;
}

// ------------- ratio and error -------------
struct Rpt { double x=0.0, y=0.0, ey=0.0; bool ok=false; };

static Rpt ratio_two_targets(const MeanC& A, const MeanC& D, double x){
  Rpt p; p.x=x;
  if (!A.ok || !D.ok || D.c==0) return p;
  const double R  = A.c / D.c;
  const double eA = A.ec;
  const double eD = D.ec;
  const double rel2 = (eA*A.c>0 ? (eA*eA)/(A.c*A.c) : 0.0)
                    + (eD*eD)/(D.c*D.c);
  p.y = R;
  p.ey = R * std::sqrt(std::max(0.0, rel2));
  p.ok = std::isfinite(p.y) && std::isfinite(p.ey);
  return p;
}

// ------------- plot helpers -------------
static void styleGraph(TGraphErrors* g, Color_t col, Style_t mstyle=20){
  g->SetMarkerStyle(mstyle);
  g->SetMarkerSize(1.0);
  g->SetMarkerColor(col);
  g->SetLineColor(col);
}

static void drawFrameAndTitle(TCanvas& c,
                              double xmin, double xmax,
                              double ymin, double ymax,
                              const std::string& titleLine){
  TH1F* frame = (TH1F*)c.DrawFrame(xmin, ymin, xmax, ymax);
  frame->GetXaxis()->SetTitle("z");
  frame->GetYaxis()->SetTitle("<cos2#phi_{h}>_{A}/<cos2#phi_{h}>_{D}");
  frame->SetTitle("");
  TLatex L; L.SetNDC(); L.SetTextSize(0.035);
  L.DrawLatex(0.12, 0.94, titleLine.c_str());
  // horizontal line at 1
  TLine l; l.SetLineColor(kGray+2); l.SetLineStyle(7);
  l.DrawLine(xmin, 1.0, xmax, 1.0);
}

// ------------- main -------------
int main(int argc, char** argv){
  parseArgs(argc, argv);
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);

  // pick xB–theta grid
  const Grid& G = (cfg.scheme=="DAT") ? gridDAT : gridBIS;
  const int Nx  = (int)G.x.size()-1;  // 7 for DAT
  const int Nth = (int)G.th.size()-1; // 3
  const int Nz  = 5;
  const int Npt = 5;

  // read all targets
  std::unordered_map<Key,Row,KeyHash,KeyEq> mCxC, mCu, mSn, mLD2, mC1, mC2;
  auto path = [&](const std::string& T){ return "C2counts_" + cfg.scheme + "_" + T + ".txt"; };
  if (!readCcounts(path("CxC"), mCxC)) return 10;
  if (!readCcounts(path("Cu"),  mCu )) return 11;
  if (!readCcounts(path("Sn"),  mSn )) return 12;
  if (!readCcounts(path("LD2"), mLD2)) return 13;
  if (!readCcounts(path("C1"),  mC1 )) return 14;
  if (!readCcounts(path("C2"),  mC2 )) return 15;

  // ---------- A/LD2 PDFs ----------
  for (int ipt=1; ipt<=Npt; ++ipt){
    const std::string pdf = "C2ratio_" + cfg.scheme + "_A_over_LD2_pt" + std::to_string(ipt) + ".pdf";
    int page=0;
    const int Npages = Nx*Nth; // 21

    for (int ith=1; ith<=Nth; ++ith){
      for (int ixb=1; ixb<=Nx; ++ixb){
        // build three graphs (CxC, Cu, Sn), each across z=1..5
        auto* gCxC = new TGraphErrors(); gCxC->SetName("gCxC");
        auto* gCu  = new TGraphErrors(); gCu ->SetName("gCu");
        auto* gSn  = new TGraphErrors(); gSn ->SetName("gSn");
        styleGraph(gCxC, kBlack);
        styleGraph(gCu,  kGreen+2);
        styleGraph(gSn,  kRed+1);

        int pc=0, pu=0, ps=0;
        for (int iz=1; iz<=Nz; ++iz){
          const Key k{ixb,ith,iz,ipt};
          // LD2 mean
          MeanC mD; auto itD = mLD2.find(k);
          if (itD!=mLD2.end()) mD = meanCos_from_Row(itD->second);

          auto putPoint = [&](const std::unordered_map<Key,Row,KeyHash,KeyEq>& M,
                              TGraphErrors* g, int& idx, double xoff){
            auto it = M.find(k);
            if (it==M.end()) return;
            MeanC mA = meanCos_from_Row(it->second);
            Rpt r = ratio_two_targets(mA, mD, mid(Z_EDGES,iz)+xoff);
            if (!r.ok) return;
            g->SetPoint(idx, r.x, r.y);
            g->SetPointError(idx, 0.0, r.ey);
            ++idx;
          };

          putPoint(mCxC, gCxC, pc, 0.00);
          putPoint(mCu,  gCu,  pu, 0.01);
          putPoint(mSn,  gSn,  ps, 0.02);
        }

        // draw one page
        TCanvas c("c","",900,650);
        const double x1 = Z_EDGES[0], x2 = Z_EDGES[5];
        drawFrameAndTitle(
          c, x1, x2, 0.0, 2.0,
          Form("<cos2#phi_{h}> ratio  A/LD2  |  scheme=%s  |  x_{B} [%g,%g]  #theta_{e} [%g,%g] deg  p_{T}^{2} [%g,%g] GeV^{2}",
               cfg.scheme.c_str(),
               G.x[ixb-1], G.x[ixb], G.th[ith-1], G.th[ith],
               PT2_EDGES[ipt-1], PT2_EDGES[ipt])
        );

        //auto* leg = new TLegend(0.72,0.82,0.92,0.93);
        auto* leg = new TLegend(0.76, 0.70, 0.96, 0.86);

        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
        leg->AddEntry(gCxC, "CxC", "pe");
        leg->AddEntry(gCu , "Cu" , "pe");
        leg->AddEntry(gSn , "Sn" , "pe");

        if (gCxC->GetN()>0) gCxC->Draw("PE SAME");
        if (gCu ->GetN()>0) gCu ->Draw("PE SAME");
        if (gSn ->GetN()>0) gSn ->Draw("PE SAME");
        leg->Draw();

        if (page==0) c.Print((pdf+"(").c_str());
        else if (page==Npages-1) c.Print((pdf+")").c_str());
        else c.Print(pdf.c_str());
        ++page;

        delete gCxC; delete gCu; delete gSn; delete leg;
      }
    }
    std::cout << "[OK] wrote " << pdf << "\n";
  }

  // ---------- C1/C2 PDFs ----------
  for (int ipt=1; ipt<=Npt; ++ipt){
    const std::string pdf = "C2ratio_" + cfg.scheme + "_C1_over_C2_pt" + std::to_string(ipt) + ".pdf";
    int page=0;
    const int Npages = Nx*Nth;

    for (int ith=1; ith<=Nth; ++ith){
      for (int ixb=1; ixb<=Nx; ++ixb){
        auto* gC = new TGraphErrors(); styleGraph(gC, kBlue+1);

        int pi=0;
        for (int iz=1; iz<=Nz; ++iz){
          const Key k{ixb,ith,iz,ipt};
          auto it1 = mC1.find(k);
          auto it2 = mC2.find(k);
          if (it1==mC1.end() || it2==mC2.end()) continue;

          MeanC m1 = meanCos_from_Row(it1->second);
          MeanC m2 = meanCos_from_Row(it2->second);
          Rpt r = ratio_two_targets(m1, m2, mid(Z_EDGES,iz));
          if (!r.ok) continue;

          gC->SetPoint(pi, r.x, r.y);
          gC->SetPointError(pi, 0.0, r.ey);
          ++pi;
        }

        TCanvas c("c","",900,650);
        const double x1 = Z_EDGES[0], x2 = Z_EDGES[5];
        drawFrameAndTitle(
          c, x1, x2, 0.0, 2.0,
          Form("<cos2#phi_{h}> self-ratio  C1/C2  |  scheme=%s  |  x_{B} [%g,%g]  #theta_{e} [%g,%g] deg  p_{T}^{2} [%g,%g] GeV^{2}",
               cfg.scheme.c_str(),
               G.x[ixb-1], G.x[ixb], G.th[ith-1], G.th[ith],
               PT2_EDGES[ipt-1], PT2_EDGES[ipt])
        );

        if (gC->GetN()>0) gC->Draw("PE SAME");
        // no legend by request

        if (page==0) c.Print((pdf+"(").c_str());
        else if (page==Npages-1) c.Print((pdf+")").c_str());
        else c.Print(pdf.c_str());
        ++page;

        delete gC;
      }
    }
    std::cout << "[OK] wrote " << pdf << "\n";
  }

  return 0;
}
