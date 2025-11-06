//computeRplots.cpp 
//
////ouput is the plotting of R as a fucnction of z for different kinematic bins
//it reads from the txt files with the coded bin numbers and counts (number of events) of target A and deuterium
//os it will take the target A you want but also the deuterium counts from the LD2 files
//
//g++ computeRplots.cpp -o computeRplots $(root-config --cflags --libs)
//# Default scheme DAT, all targets
//./computeRplots
//
//# BIS scheme, target Sn
//./computeRplots --scheme BIS 
//
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TString.h>
#include <TLine.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

struct Cfg {
  std::string scheme = "DAT"; // DAT or BIS
} cfg;

void usage(const char* p){
  std::cerr << "Usage: " << p << " [--scheme DAT|BIS]\n"
            << "Reads Ecounts/Hcounts tables for CxC, Cu, Sn, and LD2.\n"
            << "Writes: Rcompare_<SCHEME>_phiK.pdf, K=1..5\n";
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

/* ---------- binning edges ---------- */
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

/* ---------- R and error ---------- */
struct Rpt { double z=0, R=0, eR=0; bool ok=false; };
static inline Rpt makeR(double NAh,double NAe,double NDh,double NDe,double zmid){
  Rpt p; p.z=zmid;
  if (NAh>0 && NAe>0 && NDh>0 && NDe>0){
    const double num=NAh/NAe, den=NDh/NDe;
    p.R = num/den;
    p.eR = p.R * std::sqrt(1.0/NAh + 1.0/NDh + 1.0/NAe + 1.0/NDe);
    p.ok = true;
  } else { p.R=0; p.eR=0; p.ok=false; }
  return p;
}

/* ---------- one page draw ---------- */
void drawPage(const std::vector<Rpt>& CxC, const std::vector<Rpt>& Cu,
              const std::vector<Rpt>& Sn, int xb,int th,int ipt,int iphi,
              const Grid& G, const std::string& outPdf,
              bool first,bool last)
{
  TCanvas c("c","",900,650);
  gStyle->SetOptStat(0);

  // Fixed axes: Y in [0,2], X from z_min to z_max (+tiny margin for x-jitter)
  const double zmin = Z_EDGES[0];
  const double zmax = Z_EDGES[5] + 0.03; // allow +0.02 jitter at right edge
  TH1F* frame = (TH1F*)c.DrawFrame(zmin, 0.0, zmax, 1.5);
  frame->GetXaxis()->SetTitle("z");
  frame->GetYaxis()->SetTitle("R_{A}^{h} (A vs LD2)");
  frame->SetTitle("");

  // header
  const double xL=G.x[xb-1], xH=G.x[xb];
  const double tL=G.th[th-1], tH=G.th[th];
  const double pL=PT2_EDGES[ipt-1], pH=PT2_EDGES[ipt];
  const double phL=PHI_EDGES[iphi-1], phH=PHI_EDGES[iphi];
  const bool lastTheta = (th==(int)G.th.size()-1);
  TLatex latex; latex.SetNDC(); latex.SetTextSize(0.035);
  latex.DrawLatex(0.12,0.94,
    Form("R_{A}^{h}(z)  |  Scheme=%s  |  x_{B} [%g,%g]  #theta_{e} [%g,%g%s deg  p_{T}^{2} [%g,%g]  #phi_{h} [%g,%g]",
         cfg.scheme.c_str(), xL,xH, tL,tH, lastTheta?"]":"]", pL,pH, phL,phH));

  // horizontal reference line at R=1
  TLine ref; ref.SetLineStyle(2); ref.SetLineWidth(2);
  ref.DrawLine(zmin, 1.0, Z_EDGES[5]+0.5, 1.0);

  // Helper: build graph with an x-jitter (dx) so points don’t overlap
  auto makeGraph = [](const std::vector<Rpt>& V, int mcolor, double dx){
    std::vector<double> X,Y,EY; X.reserve(5); Y.reserve(5); EY.reserve(5);
    for (auto& p: V) if(p.ok){
      X.push_back(p.z + dx);
      Y.push_back(p.R);
      EY.push_back(p.eR);
    }
    auto* g = new TGraphErrors((int)X.size());
    for (int i=0;i<(int)X.size();++i){ g->SetPoint(i,X[i],Y[i]); g->SetPointError(i,0.0,EY[i]); }
    g->SetMarkerStyle(20); g->SetMarkerSize(1.0);
    g->SetLineColor(mcolor); g->SetMarkerColor(mcolor);
    return g;
  };

  // Jitter: CxC at z, Cu at z+0.01, Sn at z+0.02
  auto* gCxC = makeGraph(CxC, kBlack,     0.00);
  auto* gCu  = makeGraph(Cu,  kGreen+2,   0.01);
  auto* gSn  = makeGraph(Sn,  kRed,       0.02);

  if (gCxC->GetN()>0) gCxC->Draw("PE SAME");
  if (gCu->GetN()>0)  gCu->Draw("PE SAME");
  if (gSn->GetN()>0)  gSn->Draw("PE SAME");

  // Legend
  TLegend leg(0.78,0.78,0.93,0.93);
  leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.03);
  leg.AddEntry(gCxC,"CxC","pe");
  leg.AddEntry(gCu, "Cu", "pe");
  leg.AddEntry(gSn, "Sn", "pe");
  leg.Draw();

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

  // Targets to compare
  const std::vector<std::string> T = {"CxC","Cu","Sn"};

  // Load all maps (CxC, Cu, Sn, LD2)
  std::unordered_map<std::string,Map2> E;
  std::unordered_map<std::string,Map5> H;
  auto needE = [&](const std::string& t){ return "Ecounts_"+cfg.scheme+"_"+t+".txt"; };
  auto needH = [&](const std::string& t){ return "Hcounts_"+cfg.scheme+"_"+t+".txt"; };

  // LD2 reference
  Map2 eLD2; Map5 hLD2;
  if (!readE(needE("LD2"), eLD2)) { std::cerr << "[FATAL] Cannot read " << needE("LD2") << "\n"; return 10; }
  if (!readH(needH("LD2"), hLD2)) { std::cerr << "[FATAL] Cannot read " << needH("LD2") << "\n"; return 11; }

  // nuclear targets
  for (auto& t: T){
    if (!readE(needE(t), E[t])) { std::cerr << "[FATAL] Cannot read " << needE(t) << "\n"; return 12; }
    if (!readH(needH(t), H[t])) { std::cerr << "[FATAL] Cannot read " << needH(t) << "\n"; return 13; }
  }

  // Loop over phi: one PDF per phi bin
  for (int iphi=1; iphi<=Nphi; ++iphi){
    const std::string outPdf = "Rcompare_" + cfg.scheme + "_phi" + std::to_string(iphi) + ".pdf";
    const int Npages = Nx * Nth * Npt;
    int page = 0;

    for (int th=1; th<=Nth; ++th){
      for (int xb=1; xb<=Nx; ++xb){
        for (int ipt=1; ipt<=Npt; ++ipt){

          // Prepare three vectors of 5 points (one per z)
          std::vector<Rpt> R_CxC, R_Cu, R_Sn;
          R_CxC.reserve(Nz); R_Cu.reserve(Nz); R_Sn.reserve(Nz);

          const std::string ke = k2(xb,th);
          const double NDe = (eLD2.count(ke)? eLD2.at(ke).N : 0.0);

          for (int iz=1; iz<=Nz; ++iz){
            const std::string k = k5(xb,th,iz,ipt,iphi);
            const double NDh = (hLD2.count(k)? hLD2.at(k).N : 0.0);
            const double zmid = mid(Z_EDGES, iz);

            // For each target: get Ae (same ke) and Ah (same k)
            {
              const double NAe = (E["CxC"].count(ke)? E["CxC"].at(ke).N : 0.0);
              const double NAh = (H["CxC"].count(k)?  H["CxC"].at(k).N  : 0.0);
              R_CxC.push_back(makeR(NAh,NAe,NDh,NDe,zmid));
            }
            {
              const double NAe = (E["Cu"].count(ke)? E["Cu"].at(ke).N : 0.0);
              const double NAh = (H["Cu"].count(k)?  H["Cu"].at(k).N  : 0.0);
              R_Cu.push_back(makeR(NAh,NAe,NDh,NDe,zmid));
            }
            {
              const double NAe = (E["Sn"].count(ke)? E["Sn"].at(ke).N : 0.0);
              const double NAh = (H["Sn"].count(k)?  H["Sn"].at(k).N  : 0.0);
              R_Sn.push_back(makeR(NAh,NAe,NDh,NDe,zmid));
            }
          }

          const bool first = (page==0);
          const bool last  = (page==Npages-1);
          drawPage(R_CxC,R_Cu,R_Sn, xb,th,ipt,iphi, G, outPdf, first,last);
          ++page;
        }
      }
    }
    std::cout << "[OK] Wrote " << outPdf << " (" << Npages << " pages)\n";
  }

  return 0;
}
