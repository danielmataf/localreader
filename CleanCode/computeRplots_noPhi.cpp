// computeRplots_noPhi.cpp

// Plot R for A={CxC, Cu, Sn} vs LD2 AND slef R C1 C2 pver z  from 4D tables (no phi_h).
// One PDF per pT^2 bin: R_noPhi_<SCHEME>_ptK.pdf, each page = one (xB,theta) bin.
// Y in [0,1.5], horizontal dashed line at 1. Points only (no connecting lines).
// how 2 compile
// g++ -std=c++17 computeRplots_noPhi.cpp -o computeRplots_noPhi $(root-config --cflags --libs)
//
// how 2 run
// 
// # DAT is default; produces:
// #  R_noPhi_DAT_pt1.pdf ... R_noPhi_DAT_pt5.pdf
// ./computeRplots_noPhi
// 
// # BIS grids (if you produced BIS tables):
// ./computeRplots_noPhi --scheme BIS
// 
// One PDF per pT^2 bin:
//   - Triple vs LD2:  R_noPhi_<SCHEME>_ptK.pdf
//   - Carbon split:   RC12_noPhi_<SCHEME>_ptK.pdf
// Each page = one (xB, theta) bin; Y in [0,1.5]; dotted line at 1; markers only; error bars in marker color.

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

struct Cfg { std::string scheme="DAT"; } cfg;

static void usage(const char* p){ std::cerr<<"Usage: "<<p<<" [--scheme DAT|BIS]\n"; }
static void parseArgs(int argc, char** argv){
  for (int i=1;i<argc;++i){
    std::string a=argv[i];
    if (a=="--scheme" && i+1<argc) cfg.scheme=argv[++i];
    else if (a=="-h"||a=="--help"){ usage(argv[0]); std::exit(0); }
    else { std::cerr<<"[WARN] Unknown option: "<<a<<"\n"; usage(argv[0]); std::exit(1); }
  }
  for (auto& c: cfg.scheme) c = std::toupper(c);
  if (cfg.scheme!="DAT" && cfg.scheme!="BIS"){ std::cerr<<"[FATAL] --scheme must be DAT or BIS\n"; std::exit(2); }
}

static const double Z_EDGES[6]   = {0.30, 0.38, 0.46, 0.54, 0.62, 0.70};
static const double PT2_EDGES[6] = {0.00, 0.24, 0.48, 0.72, 0.96, 1.20};

struct Grid { std::vector<double> x, th; };
static const Grid gridDAT{
  {0.10, 0.11, 0.13, 0.16, 0.20, 0.25, 0.36, 1.00},
  {7.6, 8.4, 10.3, 23.0}
};
static const Grid gridBIS{
  {0.10, 0.11, 0.15, 0.19, 0.29, 1.00},
  {7.6, 8.8, 11.0, 23.0}
};
static inline double mid(const double* e, int i){ return 0.5*(e[i-1]+e[i]); }
static inline double mid(const std::vector<double>& v, int i){ return 0.5*(v[i-1]+v[i]); }

static inline std::string k2(int xb,int th){ return std::to_string(xb)+" "+std::to_string(th); }
static inline std::string k4(int xb,int th,int z,int pt){
  return std::to_string(xb)+" "+std::to_string(th)+" "+std::to_string(z)+" "+std::to_string(pt);
}

struct C2 { double N=0; };
struct C4 { double N=0; };
using Map2 = std::unordered_map<std::string,C2>;
using Map4 = std::unordered_map<std::string,C4>;

static bool readE(const std::string& path, Map2& m){
  std::ifstream f(path); if(!f) return false;
  std::string line; 
  while (std::getline(f,line)){
    if (line.empty() || line[0]=='#') continue;
    std::istringstream iss(line);
    int xb, th; double N;
    if (!(iss>>xb>>th>>N)) continue;
    m[k2(xb,th)].N = N;
  }
  return true;
}
static bool readH4(const std::string& path, Map4& m){
  std::ifstream f(path); if(!f) return false;
  std::string line;
  while (std::getline(f,line)){
    if (line.empty() || line[0]=='#') continue;
    std::istringstream iss(line);
    int xb, th, z, pt; double N;
    if (!(iss>>xb>>th>>z>>pt>>N)) continue;
    m[k4(xb,th,z,pt)].N = N;
  }
  return true;
}

struct Rpt { double z=0, R=0, eR=0; bool ok=false; };

static inline Rpt makeR(double Nah, double Nae, double Ndh, double Nde, double zmid){
  Rpt p; p.z=zmid; p.ok=false;
  if (Nah>0 && Nae>0 && Ndh>0 && Nde>0){
    const double Ra = Nah/Nae, Rd = Ndh/Nde;
    p.R  = Ra/Rd;
    p.eR = p.R * std::sqrt(1.0/Nah + 1.0/Nae + 1.0/Ndh + 1.0/Nde);
    p.ok = true;
  }
  return p;
}
static inline Rpt makeR_C1C2(double N1h,double N1e,double N2h,double N2e,double zmid){
  Rpt p; p.z=zmid; p.ok=false;
  if (N1h>0 && N1e>0 && N2h>0 && N2e>0){
    const double R1 = N1h/N1e, R2 = N2h/N2e;
    p.R  = R1/R2;
    p.eR = p.R * std::sqrt(1.0/N1h + 1.0/N1e + 1.0/N2h + 1.0/N2e);
    p.ok = true;
  }
  return p;
}

static TGraphErrors* makeGraph(const std::vector<Rpt>& V, Color_t color, double dx){
  std::vector<double> X,Y,EY; X.reserve(V.size()); Y.reserve(V.size()); EY.reserve(V.size());
  for (auto& p: V) if (p.ok){ X.push_back(p.z + dx); Y.push_back(p.R); EY.push_back(p.eR); }
  auto* g = new TGraphErrors((int)X.size());
  for (int i=0;i<(int)X.size();++i){ g->SetPoint(i,X[i],Y[i]); g->SetPointError(i,0.0,EY[i]); }
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.0);
  g->SetMarkerColor(color);
  g->SetLineColor(color);     // <-- error bars match marker color
  g->SetLineWidth(1);
  return g;
}

static void drawTriplePage(
  const std::vector<Rpt>& RCxC, const std::vector<Rpt>& RCu, const std::vector<Rpt>& RSn,
  int xb, int th, int ipt, const Grid& G, const std::string& outPdf,
  bool first, bool last)
{
  TCanvas c("c","",900,650);
  gStyle->SetOptStat(0);

  const double zmin = Z_EDGES[0], zmax = Z_EDGES[5] + 0.03;
  TH1F* f = (TH1F*)c.DrawFrame(zmin, 0.0, zmax, 1.5);
  f->GetXaxis()->SetTitle("z");
  f->GetYaxis()->SetTitle("R_{A}^{h}");
  f->SetTitle("");

  const double xL=G.x[xb-1], xH=G.x[xb];
  const double tL=G.th[th-1], tH=G.th[th];
  const double pL=PT2_EDGES[ipt-1], pH=PT2_EDGES[ipt];
  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
  lat.DrawLatex(0.12,0.94,
    Form("Scheme=%s  |  x_{B} [%g,%g]  #theta_{e} [%g,%g] deg  p_{T}^{2} [%g,%g] GeV^{2}/c^{2}",
         cfg.scheme.c_str(), xL,xH, tL,tH, pL,pH));

  TLine ref; ref.SetLineStyle(2); ref.SetLineWidth(2);
  ref.DrawLine(Z_EDGES[0], 1.0, Z_EDGES[5]+0.03, 1.0);

  auto* gCxC = makeGraph(RCxC, kBlack,   0.00);
  auto* gCu  = makeGraph(RCu,  kGreen+2, 0.01);
  auto* gSn  = makeGraph(RSn,  kRed,     0.02);

  if (gCxC->GetN()>0) gCxC->Draw("P E SAME");
  if (gCu->GetN()>0)  gCu->Draw("P E SAME");
  if (gSn->GetN()>0)  gSn->Draw("P E SAME");

  TLegend leg(0.78,0.78,0.93,0.93);
  leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.03);
  leg.AddEntry(gCxC,"CxC","p");
  leg.AddEntry(gCu, "Cu", "p");
  leg.AddEntry(gSn, "Sn", "p");
  leg.Draw();

  if (first) c.Print((outPdf + "(").c_str());
  else if (last) c.Print((outPdf + ")").c_str());
  else c.Print(outPdf.c_str());
}

static void drawC12Page(
  const std::vector<Rpt>& RC12,
  int xb, int th, int ipt, const Grid& G, const std::string& outPdf,
  bool first, bool last)
{
  TCanvas c("c","",900,650);
  gStyle->SetOptStat(0);

  const double zmin = Z_EDGES[0], zmax = Z_EDGES[5] + 0.03;
  TH1F* f = (TH1F*)c.DrawFrame(zmin, 0.5, zmax, 1.5);
  f->GetXaxis()->SetTitle("z");
  f->GetYaxis()->SetTitle("R_{C1/C2}^{h}");
  f->SetTitle("");

  const double xL=G.x[xb-1], xH=G.x[xb];
  const double tL=G.th[th-1], tH=G.th[th];
  const double pL=PT2_EDGES[ipt-1], pH=PT2_EDGES[ipt];
  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
  lat.DrawLatex(0.12,0.94,
    Form("Scheme=%s  |  x_{B} [%g,%g]  #theta_{e} [%g,%g] deg  p_{T}^{2} [%g,%g] GeV^{2}/c^{2}",
         cfg.scheme.c_str(), xL,xH, tL,tH, pL,pH));

  TLine ref; ref.SetLineStyle(2); ref.SetLineWidth(2);
  ref.DrawLine(Z_EDGES[0], 1.0, Z_EDGES[5]+0.03, 1.0);

  auto* gC12 = makeGraph(RC12, kBlack, 0.00);
  if (gC12->GetN()>0) gC12->Draw("P E SAME");

  if (first) c.Print((outPdf + "(").c_str());
  else if (last) c.Print((outPdf + ")").c_str());
  else c.Print(outPdf.c_str());
}

int main(int argc, char** argv){
  parseArgs(argc, argv);
  gROOT->SetBatch(true);

  const Grid& G = (cfg.scheme=="DAT") ? gridDAT : gridBIS;
  const int Nx  = (int)G.x.size()-1;
  const int Nth = (int)G.th.size()-1; // 3
  const int Nz  = 5;
  const int Npt = 5;

  auto E = [&](const std::string& t){ return "Ecounts_"+cfg.scheme+"_"+t+".txt"; };
  auto H = [&](const std::string& t){ return "Hcounts4D_"+cfg.scheme+"_"+t+".txt"; };

  Map2 eLD2, eCxC, eCu, eSn, eC1, eC2;
  Map4 hLD2, hCxC, hCu, hSn, hC1, hC2;

  if (!readE(E("LD2"), eLD2) || !readH4(H("LD2"), hLD2)) { std::cerr<<"[FATAL] Missing LD2 tables\n"; return 10; }
  if (!readE(E("CxC"), eCxC) || !readH4(H("CxC"), hCxC)) { std::cerr<<"[FATAL] Missing CxC tables\n"; return 11; }
  if (!readE(E("Cu"),  eCu ) || !readH4(H("Cu"),  hCu )) { std::cerr<<"[FATAL] Missing Cu tables\n";  return 12; }
  if (!readE(E("Sn"),  eSn ) || !readH4(H("Sn"),  hSn )) { std::cerr<<"[FATAL] Missing Sn tables\n";  return 13; }
  if (!readE(E("C1"),  eC1 ) || !readH4(H("C1"),  hC1 )) { std::cerr<<"[FATAL] Missing C1 tables\n";  return 14; }
  if (!readE(E("C2"),  eC2 ) || !readH4(H("C2"),  hC2 )) { std::cerr<<"[FATAL] Missing C2 tables\n";  return 15; }

  for (int ipt=1; ipt<=Npt; ++ipt){
    const std::string outTriple = "R_noPhi_"   + cfg.scheme + "_pt" + std::to_string(ipt)  + ".pdf";
    const std::string outC12    = "RC12_noPhi_"+ cfg.scheme + "_pt" + std::to_string(ipt)  + ".pdf";
    const int Npages = Nx * Nth;
    int page = 0;

    for (int th=1; th<=Nth; ++th){
      for (int xb=1; xb<=Nx; ++xb){
        std::vector<Rpt> RCxC, RCu, RSn, RC12;
        RCxC.reserve(Nz); RCu.reserve(Nz); RSn.reserve(Nz); RC12.reserve(Nz);

        const std::string kE = k2(xb,th);
        const double Ne_LD2 = eLD2.count(kE) ? eLD2.at(kE).N : 0.0;
        const double Ne_CxC = eCxC.count(kE) ? eCxC.at(kE).N : 0.0;
        const double Ne_Cu  = eCu .count(kE) ? eCu .at(kE).N : 0.0;
        const double Ne_Sn  = eSn .count(kE) ? eSn .at(kE).N : 0.0;
        const double Ne_C1  = eC1 .count(kE) ? eC1 .at(kE).N : 0.0;
        const double Ne_C2  = eC2 .count(kE) ? eC2 .at(kE).N : 0.0;

        for (int iz=1; iz<=Nz; ++iz){
          const std::string k4key = k4(xb,th,iz,ipt);
          const double Nh_CxC = hCxC.count(k4key) ? hCxC.at(k4key).N : 0.0;
          const double Nh_Cu  = hCu .count(k4key) ? hCu .at(k4key).N : 0.0;
          const double Nh_Sn  = hSn .count(k4key) ? hSn .at(k4key).N : 0.0;
          const double Nh_LD2 = hLD2.count(k4key) ? hLD2.at(k4key).N : 0.0;
          const double Nh_C1  = hC1 .count(k4key) ? hC1 .at(k4key).N : 0.0;
          const double Nh_C2  = hC2 .count(k4key) ? hC2 .at(k4key).N : 0.0;

          RCxC.push_back( makeR(Nh_CxC, Ne_CxC, Nh_LD2, Ne_LD2, mid(Z_EDGES,iz)) );
          RCu .push_back( makeR(Nh_Cu , Ne_Cu , Nh_LD2, Ne_LD2, mid(Z_EDGES,iz)) );
          RSn .push_back( makeR(Nh_Sn , Ne_Sn , Nh_LD2, Ne_LD2, mid(Z_EDGES,iz)) );
          RC12.push_back( makeR_C1C2(Nh_C1, Ne_C1, Nh_C2, Ne_C2, mid(Z_EDGES,iz)) );
        }

        const bool first = (page==0);
        const bool last  = (page==Npages-1);
        drawTriplePage(RCxC, RCu, RSn, xb, th, ipt, G, outTriple, first, last);
        drawC12Page   (RC12,         xb, th, ipt, G, outC12,    first, last);
        ++page;
      }
    }
    std::cout << "[OK] Wrote " << outTriple << " and " << outC12
              << " (" << Npages << " pages each)\n";
  }
  return 0;
}
