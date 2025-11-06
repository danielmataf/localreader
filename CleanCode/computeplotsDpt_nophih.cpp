// computeplotsDpt_nophih.cpp
// Build: 
// g++ -std=c++17 computeplotsDpt_nophih.cpp -o computeplotsDpt_nophih $(root-config --cflags --libs)

// How 2 Run:
//   ./computeplotsDpt_nophih --scheme DAT
//   ./computeplotsDpt_nophih --scheme BIS

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLegend.h>

#include <TColor.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

struct Cfg {
  std::string scheme = "DAT";   // DAT (default) or BIS
} cfg;

static void usage(const char* prog){
  std::cerr << "Usage: " << prog << " [--scheme DAT|BIS]\n";
}

static void parseArgs(int argc, char** argv){
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

/* ---------- binning (must match producers) ---------- */
static const double Z_EDGES[6]   = {0.30, 0.38, 0.46, 0.54, 0.62, 0.70};

struct Grid { std::vector<double> x, th; };
static const Grid gridDAT{
  {0.10, 0.11, 0.13, 0.16, 0.20, 0.25, 0.36, 1.00}, // 7
  {7.6, 8.4, 10.3, 23.0}                            // 3
};
static const Grid gridBIS{
  {0.10, 0.11, 0.15, 0.19, 0.29, 1.00},             // 5
  {7.6, 8.8, 11.0, 23.0}                            // 3
};

static inline double mid(const double* e, int i/*1..5*/){ return 0.5*(e[i-1]+e[i]); }
static inline std::string key3(int xb,int th,int z){ return std::to_string(xb)+" "+std::to_string(th)+" "+std::to_string(z); }

/* ---------- tables ---------- */
struct Row { double N=0, Nw=0, Nw2=0; };
using Table = std::unordered_map<std::string, Row>;

static bool readTable(const std::string& path, Table& out){
  std::ifstream f(path);
  if(!f) return false;
  std::string line;
  while (std::getline(f,line)){
    if (line.empty() || line[0]=='#') continue;
    std::istringstream iss(line);
    int xb,th,z; double N,Nw,Nw2;
    if (!(iss>>xb>>th>>z>>N>>Nw>>Nw2)) continue;
    out[key3(xb,th,z)] = Row{N,Nw,Nw2};
  }
  return true;
}

/* ---------- math helpers ---------- */
struct AvgVar {
  double avg=0, var=0, err=0; bool ok=false;
};
static AvgVar avg_from(const Row& r){
  AvgVar a;
  if (r.N>0){
    a.avg = r.Nw / r.N;
    double v = (r.Nw2 / r.N) - a.avg*a.avg;
    if (v < 0) v = 0; // guard small negatives
    a.var = v;
    a.err = std::sqrt(a.var / r.N);
    a.ok  = true;
  }
  return a;
}
struct Point { double x=0, y=0, ey=0; bool ok=false; };

/* ---------- plot page ---------- */
static void draw_page_delta(
  const std::vector<Point>& cx, const std::vector<Point>& cu, const std::vector<Point>& sn,
  int xb, int th, const Grid& G, const std::string& titlePrefix,
  bool first, bool last, const std::string& outPdf)
{
  TCanvas c("c","",900,650);
  gStyle->SetOptStat(0);

  const double xmin = Z_EDGES[0], xmax = Z_EDGES[5];
  const double ymin = -0.05, ymax =  0.2;
  TH1F* frame = (TH1F*)c.DrawFrame(xmin, ymin, xmax, ymax);
  frame->GetXaxis()->SetTitle("z");
  frame->GetYaxis()->SetTitle("#Delta#LT p_{T}^{2} #GT ");
  frame->SetTitle("");

  // Title
  const double xL = G.x[xb-1], xH = G.x[xb];
  const double tL = G.th[th-1], tH = G.th[th];
  TLatex L; L.SetNDC(); L.SetTextSize(0.035);
  L.DrawLatex(0.12, 0.94,
    Form("%s  |  x_{B} [%g,%g]   #theta_{e} [%g,%g%s",
         titlePrefix.c_str(),
         xL, xH, tL, tH, (th==(int)G.th.size()-1?"] deg":"] deg"))
  );

  // horizontal zero line
  TLine zline(xmin,0,xmax,0); zline.SetLineStyle(2); zline.Draw();

  // Helper to draw a series with a small z-shift (dx)
  auto draw_shifted = [&](const std::vector<Point>& v, Color_t col, double dx, const char* legLabel) -> TGraphErrors* {
    int m=0; for (auto& p: v) if (p.ok) ++m;
    if (!m) return nullptr;
    auto* g = new TGraphErrors(m);
    int i=0;
    for (auto& p: v){
      if (!p.ok) continue;
      double x = p.x + dx;
      // keep inside frame just in case
      x = std::min(std::max(x, xmin+1e-6), xmax-1e-6);
      g->SetPoint(i, x, p.y);
      g->SetPointError(i, 0.0, p.ey);
      ++i;
    }
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.0);
    g->SetLineColor(col);
    g->SetMarkerColor(col);
    g->Draw("PE SAME");
    return g;
  };

  // Draw with small per-target offsets in z
  TGraphErrors* gCx = draw_shifted(cx, kBlack,    0.00, "CxC");
  TGraphErrors* gCu = draw_shifted(cu, kGreen+2,  0.01, "Cu");
  TGraphErrors* gSn = draw_shifted(sn, kRed,      0.02, "Sn");

  // Legend (only add entries that exist)
  //TLegend leg(0.70, 0.78, 0.90, 0.92);
  TLegend leg(0.76, 0.70, 0.96, 0.86);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  if (gCx) leg.AddEntry(gCx, "CxC", "pe");
  if (gCu) leg.AddEntry(gCu, "Cu",  "pe");
  if (gSn) leg.AddEntry(gSn, "Sn",  "pe");
  if (gCx || gCu || gSn) leg.Draw();

  // Multi-page output
  if (first) c.Print((outPdf+"(").c_str());
  else if (last) c.Print((outPdf+")").c_str());
  else c.Print(outPdf.c_str());
}

/* ---------- build points for D⟨pT²⟩ vs LD2 ---------- */
static std::vector<Point> build_delta_series(const Table& A, const Table& D, int xb, int th){
  std::vector<Point> v; v.reserve(5);
  for (int iz=1; iz<=5; ++iz){
    Point p; p.x = mid(Z_EDGES, iz);
    const auto& ka = key3(xb,th,iz);
    const auto& kd = ka; // same indices
    auto ia = A.find(ka), id = D.find(kd);
    if (ia!=A.end() && id!=D.end()){
      AvgVar aA = avg_from(ia->second);
      AvgVar aD = avg_from(id->second);
      if (aA.ok && aD.ok){
        p.y  = aA.avg - aD.avg;
        p.ey = std::sqrt(aA.err*aA.err + aD.err*aD.err);
        p.ok = true;
      }
    }
    v.push_back(p);
  }
  return v;
}

/* ---------- build points for D⟨pT²⟩ C1−C2 ---------- */
static std::vector<Point> build_delta_series_C1C2(const Table& C1, const Table& C2, int xb, int th){
  std::vector<Point> v; v.reserve(5);
  for (int iz=1; iz<=5; ++iz){
    Point p; p.x = mid(Z_EDGES, iz);
    const auto& k = key3(xb,th,iz);
    auto i1 = C1.find(k), i2 = C2.find(k);
    if (i1!=C1.end() && i2!=C2.end()){
      AvgVar a1 = avg_from(i1->second);
      AvgVar a2 = avg_from(i2->second);
      if (a1.ok && a2.ok){
        p.y  = a1.avg - a2.avg;
        p.ey = std::sqrt(a1.err*a1.err + a2.err*a2.err);
        p.ok = true;
      }
    }
    v.push_back(p);
  }
  return v;
}

int main(int argc, char** argv){
  parseArgs(argc, argv);
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  const Grid& G = (cfg.scheme=="DAT") ? gridDAT : gridBIS;
  const int Nx = (int)G.x.size()-1;
  const int Nth= (int)G.th.size()-1; // 3
  const int Npages = Nx * Nth;       // DAT: 7*3=21; BIS: 5*3=15

  // Load tables
  auto path = [&](const std::string& tgt){
    return "Dcounts4D_" + cfg.scheme + "_" + tgt + ".txt";
  };

  Table T_LD2, T_CxC, T_Cu, T_Sn, T_C1, T_C2;

  if (!readTable(path("LD2"), T_LD2)){ std::cerr << "[FATAL] Missing " << path("LD2") << "\n"; return 10; }

  // For A vs LD2
  if (!readTable(path("CxC"), T_CxC)) std::cerr << "[WARN] Missing " << path("CxC") << "\n";
  if (!readTable(path("Cu"),  T_Cu )) std::cerr << "[WARN] Missing " << path("Cu" ) << "\n";
  if (!readTable(path("Sn"),  T_Sn )) std::cerr << "[WARN] Missing " << path("Sn" ) << "\n";

  // For C1 vs C2
  if (!readTable(path("C1"),  T_C1 )) std::cerr << "[WARN] Missing " << path("C1" ) << "\n";
  if (!readTable(path("C2"),  T_C2 )) std::cerr << "[WARN] Missing " << path("C2" ) << "\n";

  /* ===== PDF 1: A − LD2 (CxC, Cu, Sn) ===== */
  {
    const std::string outPdf = "Dpt_" + cfg.scheme + "_A_vs_LD2.pdf";
    int page=0;
    for (int th=1; th<=Nth; ++th){
      for (int xb=1; xb<=Nx; ++xb){
        const bool first = (page==0);
        const bool last  = (page==Npages-1);

        auto Pcx = T_CxC.empty()? std::vector<Point>{} : build_delta_series(T_CxC, T_LD2, xb, th);
        auto Pcu = T_Cu .empty()? std::vector<Point>{} : build_delta_series(T_Cu , T_LD2, xb, th);
        auto Psn = T_Sn .empty()? std::vector<Point>{} : build_delta_series(T_Sn , T_LD2, xb, th);

        draw_page_delta( Pcx, Pcu, Psn,  xb, th, G, " #Delta#LT p_{T}^{2} #GT (A - LD2)", first, last, outPdf);
        ++page;
      }
    }
    std::cout << "[OK] Wrote " << outPdf << " (" << Npages << " pages)\n";
  }

  /* ===== PDF 2: C1 − C2 ===== */
  if (!T_C1.empty() && !T_C2.empty()){
    const std::string outPdf = "Dpt_" + cfg.scheme + "_C1_minus_C2.pdf";
    int page=0;
    for (int th=1; th<=Nth; ++th){
      for (int xb=1; xb<=Nx; ++xb){
        const bool first = (page==0);
        const bool last  = (page==Npages-1);
        auto Pc12 = build_delta_series_C1C2(T_C1, T_C2, xb, th);
        // Draw single series in black for C1−C2
        std::vector<Point> empty;
        draw_page_delta( Pc12, empty, empty,xb, th, G," #Delta#LT p_{T}^{2} #GT (C1 - C2)",first, last, outPdf);
        ++page;
      }
    }
    std::cout << "[OK] Wrote " << outPdf << " (" << Npages << " pages)\n";
  } else {
    std::cout << "[INFO] Skipping C1−C2 PDF (missing C1 or C2 table)\n";
  }

  return 0;
}
