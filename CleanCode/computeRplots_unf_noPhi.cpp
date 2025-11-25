// computeRplots_unf_noPhi.cpp
//
// Plot R_A^h(z) for A = {CxC, Cu, Sn} vs LD2 from 4D tables (no phi_h),
// and overlay unfolding-corrected ratios using unfolded electron counts.
//
// Original nuclear ratios (solid markers):
//   R_A^h = (Hcounts4D_BIS_A / Ecounts_BIS_A) / (Hcounts4D_BIS_LD2 / Ecounts_BIS_LD2)
//
// Unfolded nuclear ratios (hollow markers):
//   R_A^{h,unf} = (Hcounts4D_BIS_A / unfEcounts_BIS_A)
//                 / (Hcounts4D_BIS_LD2 / unfEcounts_BIS_LD2)
//
// Also: C1/C2 self-ratio, original and unfolded:
//
//   R_{C1/C2}^h       = (N_C1^h / N_C1^e) / (N_C2^h / N_C2^e)
//   R_{C1/C2}^{h,unf} = (N_C1^h / N_C1^{e,unf}) / (N_C2^h / N_C2^{e,unf})
//
// - Only BIS scheme is supported.
// - One PDF per pT^2 bin:
//      R_unf_noPhi_BIS_ptK.pdf     (CxC, Cu, Sn vs LD2)
//      RC12_unf_noPhi_BIS_ptK.pdf  (C1/C2 self-ratio)
// - Each page = one (xB, theta) bin; X=z, Y in [0,1.5], dashed line at 1.
//
// Compile:
//   g++ -std=c++17 computeRplots_unf_noPhi.cpp -o computeRplots_unf_noPhi $(root-config --cflags --libs)
//
// Run:
//   ./computeRplots_unf_noPhi
//

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
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

struct Cfg {
  std::string scheme = "BIS"; // fixed to BIS
} cfg;

// Binning (z, pt2)
static const double Z_EDGES[6]   = {0.30, 0.38, 0.46, 0.54, 0.62, 0.70};
static const double PT2_EDGES[6] = {0.00, 0.24, 0.48, 0.72, 0.96, 1.20};

struct Grid {
  std::vector<double> x;
  std::vector<double> th;
};

// BIS grid
static const Grid gridBIS{
  {0.10, 0.11, 0.15, 0.19, 0.29, 1.00},  // 5 xB bins
  {7.6,  8.8,  11.0, 23.0}              // 3 theta bins
};

static inline double mid(const double* e, int i){ return 0.5*(e[i-1]+e[i]); }

// Keys based on bin indices
static inline std::string k2(int xb,int th){
  return std::to_string(xb) + " " + std::to_string(th);
}
static inline std::string k4(int xb,int th,int z,int pt){
  return std::to_string(xb) + " " + std::to_string(th) + " " +
         std::to_string(z)  + " " + std::to_string(pt);
}

// Structures and maps
struct C2 { double N = 0.0; }; // (xb, th)
struct C4 { double N = 0.0; }; // (xb, th, z, pt)

using Map2 = std::unordered_map<std::string,C2>;
using Map4 = std::unordered_map<std::string,C4>;

// Read Ecounts_<SCHEME>_<TARGET>.txt
static bool readE(const std::string& path, Map2& m){
  std::ifstream f(path);
  if (!f) return false;
  std::string line;
  while (std::getline(f,line)){
    if (line.empty() || line[0]=='#') continue;
    std::istringstream iss(line);
    int xb, th; double N;
    if (!(iss >> xb >> th >> N)) continue;
    m[k2(xb,th)].N = N;
  }
  return true;
}

// Read Hcounts4D_<SCHEME>_<TARGET>.txt
static bool readH4(const std::string& path, Map4& m){
  std::ifstream f(path);
  if (!f) return false;
  std::string line;
  while (std::getline(f,line)){
    if (line.empty() || line[0]=='#') continue;
    std::istringstream iss(line);
    int xb, th, z, pt; double N;
    if (!(iss >> xb >> th >> z >> pt >> N)) continue;
    m[k4(xb,th,z,pt)].N = N;
  }
  return true;
}

// Read unfEcounts_BIS_<TARGET>.txt
// Format: "# xb_bin theta_bin N" then "ix iy N"
static bool readUnfE(const std::string& path, Map2& m){
  std::ifstream f(path);
  if (!f) return false;
  std::string line;
  while (std::getline(f,line)){
    if (line.empty() || line[0]=='#') continue;
    std::istringstream iss(line);
    int xb, th; double N;
    if (!(iss >> xb >> th >> N)) continue;
    m[k2(xb,th)].N = N;
  }
  return true;
}

// Ratio point
struct Rpt {
  double z  = 0.0;
  double R  = 0.0;
  double eR = 0.0;
  bool   ok = false;
};

// Nuclear ratio A/LD2
static inline Rpt makeR(double Nah, double Nae,
                        double Ndh, double Nde,
                        double zmid)
{
  Rpt p; p.z = zmid; p.ok = false;
  if (Nah>0 && Nae>0 && Ndh>0 && Nde>0){
    const double Ra = Nah / Nae;
    const double Rd = Ndh / Nde;
    p.R  = Ra / Rd;
    p.eR = p.R * std::sqrt(1.0/Nah + 1.0/Nae + 1.0/Ndh + 1.0/Nde);
    p.ok = true;
  }
  return p;
}

// C1/C2 self-ratio
static inline Rpt makeR_C1C2(double N1h, double N1e,
                             double N2h, double N2e,
                             double zmid)
{
  Rpt p; p.z = zmid; p.ok = false;
  if (N1h>0 && N1e>0 && N2h>0 && N2e>0){
    const double R1 = N1h / N1e;
    const double R2 = N2h / N2e;
    p.R  = R1 / R2;
    p.eR = p.R * std::sqrt(1.0/N1h + 1.0/N1e + 1.0/N2h + 1.0/N2e);
    p.ok = true;
  }
  return p;
}

// Make a TGraphErrors from a vector of Rpt
static TGraphErrors* makeGraph(const std::vector<Rpt>& V,
                               Color_t color, double dx,
                               int markerStyle)
{
  std::vector<double> X, Y, EY;
  X.reserve(V.size());
  Y.reserve(V.size());
  EY.reserve(V.size());

  for (auto& p : V){
    if (!p.ok) continue;
    X.push_back(p.z + dx);
    Y.push_back(p.R);
    EY.push_back(p.eR);
  }

  auto* g = new TGraphErrors((int)X.size());
  for (int i=0;i<(int)X.size();++i){
    g->SetPoint(i,      X[i], Y[i]);
    g->SetPointError(i, 0.0,  EY[i]);
  }

  g->SetMarkerStyle(markerStyle);
  g->SetMarkerSize(1.0);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetLineWidth(1);

  return g;
}

// Draw one page for CxC/Cu/Sn vs LD2 (data + unfolded)
static void drawTriplePage(
  const std::vector<Rpt>& RCxC,      const std::vector<Rpt>& RCu,
  const std::vector<Rpt>& RSn,
  const std::vector<Rpt>& RCxC_unf,  const std::vector<Rpt>& RCu_unf,
  const std::vector<Rpt>& RSn_unf,
  int xb, int th, int ipt, const Grid& G,
  const std::string& outPdf, bool first, bool last)
{
  TCanvas c("c","",900,650);
  gStyle->SetOptStat(0);

  const double zmin = Z_EDGES[0];
  const double zmax = Z_EDGES[5] + 0.03;
  TH1F* frame = (TH1F*)c.DrawFrame(zmin, 0.0, zmax, 1.5);
  frame->GetXaxis()->SetTitle("z");
  frame->GetYaxis()->SetTitle("R_{A}^{h}");
  frame->SetTitle("");

  const double xL = G.x[xb-1];
  const double xH = G.x[xb];
  const double tL = G.th[th-1];
  const double tH = G.th[th];
  const double pL = PT2_EDGES[ipt-1];
  const double pH = PT2_EDGES[ipt];

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.12,0.94,
    Form("Scheme=BIS  |  x_{B} [%g,%g]  #theta_{e} [%g,%g] deg  p_{T}^{2} [%g,%g] GeV^{2}/c^{2}",
         xL,xH, tL,tH, pL,pH));

  TLine ref;
  ref.SetLineStyle(2);
  ref.SetLineWidth(2);
  ref.DrawLine(Z_EDGES[0], 1.0, Z_EDGES[5]+0.03, 1.0);

  // Jitter between targets, but no jitter between data/unf for same target
  auto* gCxC     = makeGraph(RCxC,      kBlack,     0.00, 20); // solid
  auto* gCxC_unf = makeGraph(RCxC_unf,  kBlack,     0.00, 24); // hollow
  auto* gCu      = makeGraph(RCu,       kGreen+2,   0.01, 20); // solid
  auto* gCu_unf  = makeGraph(RCu_unf,   kGreen+2,   0.01, 24); // hollow
  auto* gSn      = makeGraph(RSn,       kRed,       0.02, 20); // solid
  auto* gSn_unf  = makeGraph(RSn_unf,   kRed,       0.02, 24); // hollow

  if (gCxC->GetN()>0)     gCxC    ->Draw("P E SAME");
  if (gCxC_unf->GetN()>0) gCxC_unf->Draw("P E SAME");
  if (gCu->GetN()>0)      gCu     ->Draw("P E SAME");
  if (gCu_unf->GetN()>0)  gCu_unf ->Draw("P E SAME");
  if (gSn->GetN()>0)      gSn     ->Draw("P E SAME");
  if (gSn_unf->GetN()>0)  gSn_unf ->Draw("P E SAME");

  TLegend leg(0.60,0.68,0.93,0.93);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.03);
  leg.AddEntry(gCxC,    "CxC (data)", "p");
  leg.AddEntry(gCxC_unf,"CxC (unf.)", "p");
  leg.AddEntry(gCu,     "Cu (data)",  "p");
  leg.AddEntry(gCu_unf, "Cu (unf.)",  "p");
  leg.AddEntry(gSn,     "Sn (data)",  "p");
  leg.AddEntry(gSn_unf, "Sn (unf.)",  "p");
  leg.Draw();

  if (first)      c.Print((outPdf + "(").c_str());
  else if (last)  c.Print((outPdf + ")").c_str());
  else            c.Print(outPdf.c_str());
}

// Draw one page for C1/C2 (data + unfolded)
static void drawC12Page(
  const std::vector<Rpt>& RC12,
  const std::vector<Rpt>& RC12_unf,
  int xb, int th, int ipt, const Grid& G,
  const std::string& outPdf, bool first, bool last)
{
  TCanvas c("c","",900,650);
  gStyle->SetOptStat(0);

  const double zmin = Z_EDGES[0];
  const double zmax = Z_EDGES[5] + 0.03;
  TH1F* frame = (TH1F*)c.DrawFrame(zmin, 0.5, zmax, 1.5);
  frame->GetXaxis()->SetTitle("z");
  frame->GetYaxis()->SetTitle("R_{C1/C2}^{h}");
  frame->SetTitle("");

  const double xL = G.x[xb-1];
  const double xH = G.x[xb];
  const double tL = G.th[th-1];
  const double tH = G.th[th];
  const double pL = PT2_EDGES[ipt-1];
  const double pH = PT2_EDGES[ipt];

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.12,0.94,
    Form("Scheme=BIS  |  x_{B} [%g,%g]  #theta_{e} [%g,%g] deg  p_{T}^{2} [%g,%g] GeV^{2}/c^{2}",
         xL,xH, tL,tH, pL,pH));

  TLine ref;
  ref.SetLineStyle(2);
  ref.SetLineWidth(2);
  ref.DrawLine(Z_EDGES[0], 1.0, Z_EDGES[5]+0.03, 1.0);

  auto* gC12     = makeGraph(RC12,     kBlack, 0.00, 20); // solid
  auto* gC12_unf = makeGraph(RC12_unf, kBlack, 0.00, 24); // hollow

  if (gC12->GetN()>0)     gC12    ->Draw("P E SAME");
  if (gC12_unf->GetN()>0) gC12_unf->Draw("P E SAME");

  TLegend leg(0.60,0.78,0.93,0.93);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.03);
  leg.AddEntry(gC12,    "C1/C2 (data)", "p");
  leg.AddEntry(gC12_unf,"C1/C2 (unf.)", "p");
  leg.Draw();

  if (first)      c.Print((outPdf + "(").c_str());
  else if (last)  c.Print((outPdf + ")").c_str());
  else            c.Print(outPdf.c_str());
}

// ---------------------- main -------------------------

int main(int argc, char** argv){
  gROOT->SetBatch(true);

  const Grid& G = gridBIS;
  const int Nx  = (int)G.x.size()-1;  // 5 xB bins
  const int Nth = (int)G.th.size()-1; // 3 theta bins
  const int Nz  = 5;                  // 5 z bins
  const int Npt = 5;                  // 5 pt2 bins

  auto E      = [&](const std::string& t){ return "Ecounts_"     + cfg.scheme + "_" + t + ".txt"; };
  auto H4     = [&](const std::string& t){ return "Hcounts4D_"   + cfg.scheme + "_" + t + ".txt"; };
  auto E_unf  = [&](const std::string& t){ return "unfEcounts_"  + cfg.scheme + "_" + t + ".txt"; };

  // Maps for electrons (data and unfolded) and hadrons
  Map2 eLD2, eCxC, eCu, eSn, eC1, eC2;
  Map2 eLD2_unf, eCxC_unf, eCu_unf, eSn_unf, eC1_unf, eC2_unf;
  Map4 hLD2, hCxC, hCu, hSn, hC1, hC2;

  // Original electron counts
  if (!readE(E("LD2"), eLD2)) { std::cerr << "[FATAL] Missing " << E("LD2") << "\n"; return 10; }
  if (!readE(E("CxC"), eCxC)) { std::cerr << "[FATAL] Missing " << E("CxC") << "\n"; return 11; }
  if (!readE(E("Cu"),  eCu )) { std::cerr << "[FATAL] Missing " << E("Cu")  << "\n"; return 12; }
  if (!readE(E("Sn"),  eSn )) { std::cerr << "[FATAL] Missing " << E("Sn")  << "\n"; return 13; }
  if (!readE(E("C1"),  eC1 )) { std::cerr << "[FATAL] Missing " << E("C1")  << "\n"; return 14; }
  if (!readE(E("C2"),  eC2 )) { std::cerr << "[FATAL] Missing " << E("C2")  << "\n"; return 15; }

  // Unfolded electron counts
  if (!readUnfE(E_unf("LD2"), eLD2_unf)) { std::cerr << "[FATAL] Missing " << E_unf("LD2") << "\n"; return 16; }
  if (!readUnfE(E_unf("CxC"), eCxC_unf)) { std::cerr << "[FATAL] Missing " << E_unf("CxC") << "\n"; return 17; }
  if (!readUnfE(E_unf("Cu"),  eCu_unf )) { std::cerr << "[FATAL] Missing " << E_unf("Cu")  << "\n"; return 18; }
  if (!readUnfE(E_unf("Sn"),  eSn_unf )) { std::cerr << "[FATAL] Missing " << E_unf("Sn")  << "\n"; return 19; }
  if (!readUnfE(E_unf("C1"),  eC1_unf )) { std::cerr << "[FATAL] Missing " << E_unf("C1")  << "\n"; return 20; }
  if (!readUnfE(E_unf("C2"),  eC2_unf )) { std::cerr << "[FATAL] Missing " << E_unf("C2")  << "\n"; return 21; }

  // Hadron counts
  if (!readH4(H4("LD2"), hLD2)) { std::cerr << "[FATAL] Missing " << H4("LD2") << "\n"; return 22; }
  if (!readH4(H4("CxC"), hCxC)) { std::cerr << "[FATAL] Missing " << H4("CxC") << "\n"; return 23; }
  if (!readH4(H4("Cu"),  hCu )) { std::cerr << "[FATAL] Missing " << H4("Cu")  << "\n"; return 24; }
  if (!readH4(H4("Sn"),  hSn )) { std::cerr << "[FATAL] Missing " << H4("Sn")  << "\n"; return 25; }
  if (!readH4(H4("C1"),  hC1 )) { std::cerr << "[FATAL] Missing " << H4("C1")  << "\n"; return 26; }
  if (!readH4(H4("C2"),  hC2 )) { std::cerr << "[FATAL] Missing " << H4("C2")  << "\n"; return 27; }

  // Loop over pT^2 bins
  for (int ipt=1; ipt<=Npt; ++ipt){
    const std::string outTriple = "R_unf_noPhi_BIS_pt"    + std::to_string(ipt) + ".pdf";
    const std::string outC12    = "RC12_unf_noPhi_BIS_pt" + std::to_string(ipt) + ".pdf";
    const int Npages = Nx * Nth;
    int page = 0;

    for (int th=1; th<=Nth; ++th){
      for (int xb=1; xb<=Nx; ++xb){
        std::vector<Rpt> RCxC, RCu, RSn;
        std::vector<Rpt> RCxC_unf, RCu_unf, RSn_unf;
        std::vector<Rpt> RC12, RC12_unf;

        RCxC.reserve(Nz);     RCu.reserve(Nz);     RSn.reserve(Nz);
        RCxC_unf.reserve(Nz); RCu_unf.reserve(Nz); RSn_unf.reserve(Nz);
        RC12.reserve(Nz);     RC12_unf.reserve(Nz);

        const std::string keyE = k2(xb,th);

        const double Ne_LD2     = eLD2.count(keyE)     ? eLD2.at(keyE).N     : 0.0;
        const double Ne_CxC     = eCxC.count(keyE)     ? eCxC.at(keyE).N     : 0.0;
        const double Ne_Cu      = eCu .count(keyE)     ? eCu .at(keyE).N     : 0.0;
        const double Ne_Sn      = eSn .count(keyE)     ? eSn .at(keyE).N     : 0.0;
        const double Ne_C1      = eC1 .count(keyE)     ? eC1 .at(keyE).N     : 0.0;
        const double Ne_C2      = eC2 .count(keyE)     ? eC2 .at(keyE).N     : 0.0;

        const double Ne_LD2_unf = eLD2_unf.count(keyE) ? eLD2_unf.at(keyE).N : 0.0;
        const double Ne_CxC_unf = eCxC_unf.count(keyE) ? eCxC_unf.at(keyE).N : 0.0;
        const double Ne_Cu_unf  = eCu_unf .count(keyE) ? eCu_unf .at(keyE).N : 0.0;
        const double Ne_Sn_unf  = eSn_unf .count(keyE) ? eSn_unf .at(keyE).N : 0.0;
        const double Ne_C1_unf  = eC1_unf .count(keyE) ? eC1_unf .at(keyE).N : 0.0;
        const double Ne_C2_unf  = eC2_unf .count(keyE) ? eC2_unf .at(keyE).N : 0.0;

        for (int iz=1; iz<=Nz; ++iz){
          const std::string k4key = k4(xb,th,iz,ipt);

          const double Nh_CxC = hCxC.count(k4key) ? hCxC.at(k4key).N : 0.0;
          const double Nh_Cu  = hCu .count(k4key) ? hCu .at(k4key).N : 0.0;
          const double Nh_Sn  = hSn .count(k4key) ? hSn .at(k4key).N : 0.0;
          const double Nh_LD2 = hLD2.count(k4key) ? hLD2.at(k4key).N : 0.0;
          const double Nh_C1  = hC1 .count(k4key) ? hC1 .at(k4key).N : 0.0;
          const double Nh_C2  = hC2 .count(k4key) ? hC2 .at(k4key).N : 0.0;

          const double zmid = mid(Z_EDGES, iz);

          // Nuclear ratios (data)
          RCxC.push_back( makeR(Nh_CxC, Ne_CxC, Nh_LD2, Ne_LD2, zmid) );
          RCu .push_back( makeR(Nh_Cu , Ne_Cu , Nh_LD2, Ne_LD2, zmid) );
          RSn .push_back( makeR(Nh_Sn , Ne_Sn , Nh_LD2, Ne_LD2, zmid) );

          // Nuclear ratios (unfolded electrons)
          RCxC_unf.push_back( makeR(Nh_CxC, Ne_CxC_unf, Nh_LD2, Ne_LD2_unf, zmid) );
          RCu_unf .push_back( makeR(Nh_Cu , Ne_Cu_unf , Nh_LD2, Ne_LD2_unf, zmid) );
          RSn_unf .push_back( makeR(Nh_Sn , Ne_Sn_unf , Nh_LD2, Ne_LD2_unf, zmid) );

          // C1/C2 self-ratios (data and unfolded)
          RC12    .push_back( makeR_C1C2(Nh_C1, Ne_C1,     Nh_C2, Ne_C2,     zmid) );
          RC12_unf.push_back( makeR_C1C2(Nh_C1, Ne_C1_unf, Nh_C2, Ne_C2_unf, zmid) );
        }

        const bool first = (page == 0);
        const bool last  = (page == Npages-1);

        drawTriplePage(RCxC, RCu, RSn,
                       RCxC_unf, RCu_unf, RSn_unf,
                       xb, th, ipt, G, outTriple, first, last);
        drawC12Page   (RC12, RC12_unf,
                       xb, th, ipt, G, outC12,   first, last);

        ++page;
      }
    }

    std::cout << "[OK] Wrote " << outTriple << " and " << outC12
              << " (" << Npages << " pages each)\n";
  }

  return 0;
}
