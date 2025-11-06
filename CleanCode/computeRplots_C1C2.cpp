//computeRplots_C1C2.cpp
// this should do the same as R but instead we replaced the deuterium reference to calculate self ratio on both carbon targets
// the resulting Ratio should be really close to one. 

// g++ -std=c++17 computeRplots_C1C2.cpp -o computeRplots_C1C2 $(root-config --cflags --libs)
//# DAT (default)
//./computeRplots_C1C2
//
//# BIS (if you also have BIS tables for C1/C2)
//./computeRplots_C1C2 --scheme BIS



#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

struct Cfg { std::string scheme = "DAT"; } cfg;

void usage(const char* p){
  std::cerr << "Usage: " << p << " [--scheme DAT|BIS]\n"
            << "Reads Ecounts/Hcounts tables for C1 and C2.\n"
            << "Writes: R_C1overC2_<SCHEME>_phiK.pdf (K=1..5)\n";
}
void parseArgs(int argc, char** argv){
  for (int i=1;i<argc;++i){
    std::string a=argv[i];
    if (a=="--scheme" && i+1<argc) cfg.scheme = argv[++i];
    else if (a=="-h" || a=="--help"){ usage(argv[0]); std::exit(0); }
    else { std::cerr<<"[WARN] Unknown option: "<<a<<"\n"; usage(argv[0]); std::exit(1); }
  }
  for (auto& c: cfg.scheme) c = std::toupper(c);
  if (cfg.scheme!="DAT" && cfg.scheme!="BIS"){
    std::cerr<<"[FATAL] --scheme must be DAT or BIS\n"; std::exit(2);
  }
}

/* ---------- binning (as before) ---------- */
static const double Z_EDGES[6]   = {0.30, 0.38, 0.46, 0.54, 0.62, 0.70};
static const double PT2_EDGES[6] = {0.00, 0.24, 0.48, 0.72, 0.96, 1.20};
static const double PHI_EDGES[6] = {0.0, 72.0, 144.0, 216.0, 288.0, 360.0};

struct Grid { std::vector<double> x, th; };
static const Grid gridDAT{
  {0.10, 0.11, 0.13, 0.16, 0.20, 0.25, 0.36, 1.00},
  {7.6, 8.4, 10.3, 23.0}
};
static const Grid gridBIS{
  {0.10, 0.11, 0.15, 0.19, 0.29, 1.00},
  {7.6, 8.8, 11.0, 23.0}
};

static inline double mid(const double* e, int i){ return 0.5*(e[i-1]+e[i]); } // i=1..5
static inline std::string k2(int xb,int th){ return std::to_string(xb)+" "+std::to_string(th); }
static inline std::string k5(int xb,int th,int z,int pt,int ph){
  return std::to_string(xb)+" "+std::to_string(th)+" "+std::to_string(z)+" "+std::to_string(pt)+" "+std::to_string(ph);
}

/* ---------- table readers ---------- */
struct C2 { double N=0; };
struct C5 { double N=0; };
using Map2 = std::unordered_map<std::string,C2>;
using Map5 = std::unordered_map<std::string,C5>;

bool readE(const std::string& path, Map2& m){
  std::ifstream f(path); if(!f) return false;
  std::string line;
  while (std::getline(f,line)){
    if (line.empty()||line[0]=='#') continue;
    std::istringstream iss(line);
    int xb,th; double N; if(!(iss>>xb>>th>>N)) continue;
    m[k2(xb,th)].N = N;
  } return true;
}
bool readH(const std::string& path, Map5& m){
  std::ifstream f(path); if(!f) return false;
  std::string line;
  while (std::getline(f,line)){
    if (line.empty()||line[0]=='#') continue;
    std::istringstream iss(line);
    int xb,th,z,pt,ph; double N; if(!(iss>>xb>>th>>z>>pt>>ph>>N)) continue;
    m[k5(xb,th,z,pt,ph)].N = N;
  } return true;
}

/* ---------- R and error (C1 over C2) ---------- */
struct Rpt { double z=0, R=0, eR=0; bool ok=false; };
static inline Rpt makeR(double N1h,double N1e,double N2h,double N2e,double zmid){
  Rpt p; p.z=zmid;
  if (N1h>0 && N1e>0 && N2h>0 && N2e>0){
    const double num=N1h/N1e, den=N2h/N2e;
    p.R = num/den;
    p.eR = p.R * std::sqrt(1.0/N1h + 1.0/N2h + 1.0/N1e + 1.0/N2e);
    p.ok=true;
  } else { p.R=0; p.eR=0; p.ok=false; }
  return p;
}

/* ---------- draw one page (single curve) ---------- */
void drawPage(const std::vector<Rpt>& RC12, int xb,int th,int ipt,int iphi,
              const Grid& G, const std::string& outPdf, bool first,bool last)
{
  TCanvas c("c","",900,650);
  gStyle->SetOptStat(0);

  // fixed axes: Y [0,1.5]; allow a tiny x margin for safety
  const double zmin = Z_EDGES[0];
  const double zmax = Z_EDGES[5] + 0.02;
  TH1F* frame = (TH1F*)c.DrawFrame(zmin, 0.0, zmax, 1.5);
  frame->GetXaxis()->SetTitle("z");
  frame->GetYaxis()->SetTitle("R_{C1/C2}^{h}");
  frame->SetTitle("");

  // header with bin intervals
  const double xL=G.x[xb-1], xH=G.x[xb];
  const double tL=G.th[th-1], tH=G.th[th];
  const double pL=PT2_EDGES[ipt-1], pH=PT2_EDGES[ipt];
  const double phL=PHI_EDGES[iphi-1], phH=PHI_EDGES[iphi];
  const bool lastTheta = (th==(int)G.th.size()-1);
  TLatex latex; latex.SetNDC(); latex.SetTextSize(0.035);
  latex.DrawLatex(0.12,0.94,
    Form("R_{C1/C2}^{h}(z)  |  Scheme=%s  |  x_{B} [%g,%g)  #theta_{e} [%g,%g%s deg  p_{T}^{2} [%g,%g)  #phi_{h} [%g,%g)",
         cfg.scheme.c_str(), xL,xH, tL,tH, lastTheta?"]":")", pL,pH, phL,phH));

  // reference line R=1 (dotted)
  TLine ref; ref.SetLineStyle(2); ref.SetLineWidth(2);
  ref.DrawLine(Z_EDGES[0],1.0,Z_EDGES[5],1.0);

  // build graph
  std::vector<double> X,Y,EY;
  for (auto& p: RC12) if (p.ok){ X.push_back(p.z); Y.push_back(p.R); EY.push_back(p.eR); }
  auto* g = new TGraphErrors((int)X.size());
  for (int i=0;i<(int)X.size();++i){ g->SetPoint(i,X[i],Y[i]); g->SetPointError(i,0.0,EY[i]); }
  g->SetMarkerStyle(20); g->SetMarkerSize(1.0);
  g->SetLineColor(kBlack); g->SetMarkerColor(kBlack);
  if (g->GetN()>0) g->Draw("PE SAME");

  if (first) c.Print((outPdf+"(").c_str());
  else if (last) c.Print((outPdf+")").c_str());
  else c.Print(outPdf.c_str());
}

/* ---------- main ---------- */
int main(int argc, char** argv){
  parseArgs(argc, argv);
  gROOT->SetBatch(true);

  const Grid& G = (cfg.scheme=="DAT") ? gridDAT : gridBIS;
  const int Nx = (int)G.x.size()-1;
  const int Nth = (int)G.th.size()-1;
  const int Npt = 5, Nz = 5, Nphi = 5;

  auto needE = [&](const std::string& t){ return "Ecounts_"+cfg.scheme+"_"+t+".txt"; };
  auto needH = [&](const std::string& t){ return "Hcounts_"+cfg.scheme+"_"+t+".txt"; };

  // load C1 and C2
  Map2 eC1,eC2; Map5 hC1,hC2;
  if (!readE(needE("C1"), eC1)) { std::cerr << "[FATAL] Cannot read " << needE("C1") << "\n"; return 10; }
  if (!readE(needE("C2"), eC2)) { std::cerr << "[FATAL] Cannot read " << needE("C2") << "\n"; return 11; }
  if (!readH(needH("C1"), hC1)) { std::cerr << "[FATAL] Cannot read " << needH("C1") << "\n"; return 12; }
  if (!readH(needH("C2"), hC2)) { std::cerr << "[FATAL] Cannot read " << needH("C2") << "\n"; return 13; }

  for (int iphi=1; iphi<=Nphi; ++iphi){
    const std::string outPdf = "R_C1overC2_" + cfg.scheme + "_phi" + std::to_string(iphi) + ".pdf";
    const int Npages = Nx * Nth * Npt;
    int page=0;

    for (int th=1; th<=Nth; ++th){
      for (int xb=1; xb<=Nx; ++xb){
        for (int ipt=1; ipt<=Npt; ++ipt){

          std::vector<Rpt> RC12; RC12.reserve(Nz);
          const std::string ke = k2(xb,th);

          const double N1e = (eC1.count(ke)? eC1.at(ke).N : 0.0);
          const double N2e = (eC2.count(ke)? eC2.at(ke).N : 0.0);

          for (int iz=1; iz<=Nz; ++iz){
            const std::string k = k5(xb,th,iz,ipt,iphi);
            const double N1h = (hC1.count(k)? hC1.at(k).N : 0.0);
            const double N2h = (hC2.count(k)? hC2.at(k).N : 0.0);
            RC12.push_back(makeR(N1h,N1e,N2h,N2e, mid(Z_EDGES,iz)));
          }

          const bool first=(page==0), last=(page==Npages-1);
          drawPage(RC12, xb,th,ipt,iphi, G, outPdf, first,last);
          ++page;
        }
      }
    }
    std::cout << "[OK] Wrote " << outPdf << " (" << Npages << " pages)\n";
  }

  return 0;
}
