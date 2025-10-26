// make2Dkinfromhad.cpp
//
//   Produces xB–Q^2 and xB–theta 2D histograms from hadron TTrees, binned in (z, pT2, phih);
//   outputs are 5 PDFs (one per phi bin), each with 25 xB–Q^2 + 25 xB–theta pages.
//   Optionally exports one (xB–Q^2) + (xB–theta) PNG for a given 1..125 (z, pT2, phih) bin.
//
// How to compile (with ROOT):
//   g++ -O2 -std=c++17 make2Dkinfromhad.cpp `root-config --cflags --libs` -o make2Dkinfromhad
//
// How to run (example):
//   ./make2Dkinfromhad
//
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <iostream>
#include <string>
#include <algorithm>

// ===============================
// ========== INPUTS =============
// ===============================
namespace INPUTS {
  const std::string kFile = "build/Trees_Cu_RGD.root";
  const std::string kTree = "tEv_had_Cu_RGD";
  const int SINGLE_BIN_INDEX = 25;  // e.g. 73 for bin #73, -1 disables
  const std::string kOutDir = "plots_full_bins";
}

// ===============================
// ======= BIN DEFINITIONS =======
// ===============================
static const double Z_EDGES[6]   = {0.30, 0.38, 0.46, 0.54, 0.62, 0.70};
static const double PT2_EDGES[6] = {0.00, 0.24, 0.48, 0.72, 0.96, 1.20};
static const double PHI_EDGES[6] = {0.0, 72.0, 144.0, 216.0, 288.0, 360.0};

static const int XB_NBINS = 20;
static const double XB_MIN = 0.0, XB_MAX = 1;

static const int Q2_NBINS = 20;
static const double Q2_MIN = 0.0, Q2_MAX = 5.0;

static const int TH_NBINS = 20;
static const double TH_MIN = 7.0, TH_MAX = 23.0;

// ===============================
// ======= HELPER FUNCTIONS ======
// ===============================
inline void IndexToTriplet125(int idx, int &iz, int &ip, int &iph) {
  const int idx0 = idx - 1;
  iph = idx0 / 25;
  const int r = idx0 % 25;
  ip = r / 5;
  iz = r % 5;
}

TString BuildSelCell(int iz, int ip, int iph) {
  const bool lastPhi = (iph == 4);
  TString sel;
  sel.Form("z>=%g && z<%g && pt2>=%g && pt2<%g && %s",
           Z_EDGES[iz], Z_EDGES[iz + 1],
           PT2_EDGES[ip], PT2_EDGES[ip + 1],
           lastPhi ? Form("phih>=%g && phih<=%g", PHI_EDGES[iph], PHI_EDGES[iph + 1])
                   : Form("phih>=%g && phih<%g",  PHI_EDGES[iph], PHI_EDGES[iph + 1]));
  return sel;
}

TString CellTitleSuffix(int iz, int ip, int iph) {
  const bool lastPhi = (iph == 4);
  return Form("z [%.2f, %.2f)   p_{T}^{2} [%.2f, %.2f)   #phi_{h} [%.0f^{#circ}, %.0f^{#circ}%s",
              Z_EDGES[iz], Z_EDGES[iz + 1],
              PT2_EDGES[ip], PT2_EDGES[ip + 1],
              PHI_EDGES[iph], PHI_EDGES[iph + 1],
              lastPhi ? "]" : ")");
}

// --- DEBUG: print selection + a few (xb, Q2, th) tuples to stdout
void DebugPrintSample(TTree* T, const TString& sel, int maxPrint = 5) {
  Long64_t n = T->Draw("xb:Q2:th", sel, "goff");
  std::cout << "  -> selected rows: " << n << std::endl;
  const int nShow = (int)std::min<Long64_t>(n, maxPrint);
  const double* xbArr = T->GetVal(0);
  const double* q2Arr = T->GetVal(1);
  const double* thArr = T->GetVal(2);
  for (int i = 0; i < nShow; ++i) {
    std::cout << "     xb=" << xbArr[i]
              << "  Q2=" << q2Arr[i]
              << "  th=" << thArr[i] << " deg"
              << std::endl;
  }
}

#include <TLatex.h>

// Draw a 2D hist (xb vs Q2) for a selection; returns pointer so caller can delete
TH2F* DrawXBQ2(TTree *T, const TString &sel, const TString &title) {
  // Book the target histogram (so we control binning)
  TH2F *h = new TH2F("hxbq2", "", XB_NBINS, XB_MIN, XB_MAX, Q2_NBINS, Q2_MIN, Q2_MAX);
  h->GetXaxis()->SetTitle("x_{B}");
  h->GetYaxis()->SetTitle("Q^{2}  [GeV^{2}]");

  // Fill via TTree::Draw and check how many rows were selected
  Long64_t nsel = T->Draw("Q2:xb>>hxbq2", sel, "COLZ");
  gPad->SetRightMargin(0.12);

  if (nsel <= 0) {
    // No entries: draw a clean empty frame with axes and an info label
    gPad->Clear();
    TH2F *frame = new TH2F("frame_xbq2","",XB_NBINS, XB_MIN, XB_MAX, Q2_NBINS, Q2_MIN, Q2_MAX);
    frame->GetXaxis()->SetTitle("x_{B}");
    frame->GetYaxis()->SetTitle("Q^{2}  [GeV^{2}]");
    frame->SetTitle(title);
    frame->Draw();                // draw axes
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.18, 0.88, title);
    lat.DrawLatex(0.18, 0.82, "#bf{No entries in this (z, p_{T}^{2}, #phi_{h}) bin}");
    gPad->Update();
    delete h;                     // nothing filled
    return frame;                 // caller can delete this instead
  }

  // Entries exist: now set the nice title (after Draw so ROOT won’t overwrite)
  h->SetTitle(title);
  gPad->Update();
  return h;
}

// Draw a 2D hist (xb vs theta) for a selection; returns pointer so caller can delete
TH2F* DrawXBTheta(TTree *T, const TString &sel, const TString &title) {
  TH2F *h = new TH2F("hxbth", "", XB_NBINS, XB_MIN, XB_MAX, TH_NBINS, TH_MIN, TH_MAX);
  h->GetXaxis()->SetTitle("x_{B}");
  h->GetYaxis()->SetTitle("#theta_{e}  [deg]");

  Long64_t nsel = T->Draw("th:xb>>hxbth", sel, "COLZ");
  gPad->SetRightMargin(0.12);

  if (nsel <= 0) {
    gPad->Clear();
    TH2F *frame = new TH2F("frame_xbth","",XB_NBINS, XB_MIN, XB_MAX, TH_NBINS, TH_MIN, TH_MAX);
    frame->GetXaxis()->SetTitle("x_{B}");
    frame->GetYaxis()->SetTitle("#theta_{e}  [deg]");
    frame->SetTitle(title);
    frame->Draw();
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.18, 0.88, title);
    lat.DrawLatex(0.18, 0.82, "#bf{No entries in this (z, p_{T}^{2}, #phi_{h}) bin}");
    gPad->Update();
    delete h;
    return frame;
  }

  h->SetTitle(title);
  gPad->Update();
  return h;
}


// ===============================
// ========== MAIN ===============
// ===============================
int main() {
  gStyle->SetOptStat(0);
  gSystem->mkdir(INPUTS::kOutDir.c_str(), kTRUE);

  TFile fin(INPUTS::kFile.c_str(), "READ");
  if (fin.IsZombie()) {
    std::cerr << "Error: cannot open file " << INPUTS::kFile << std::endl;
    return 1;
  }
  TTree *T = (TTree *)fin.Get(INPUTS::kTree.c_str());
  if (!T) {
    std::cerr << "Error: tree " << INPUTS::kTree << " not found." << std::endl;
    return 1;
  }

  // ---- Loop over phi bins ----
  for (int iph = 0; iph < 5; ++iph) {
    TString pdf = Form("%s/phiBin%d.pdf", INPUTS::kOutDir.c_str(), iph + 1);
    TCanvas c("c", "", 900, 800);
    c.Print(pdf + "[");

    // 25 (z,pT²) bins per φh
    for (int ip = 0; ip < 5; ++ip) {
      for (int iz = 0; iz < 5; ++iz) {
        TString sel = BuildSelCell(iz, ip, iph);
        TString suffix = CellTitleSuffix(iz, ip, iph);

        // ---- DEBUG: print the bin ranges + some sample values to stdout
        //std::cout << "[phi-bin " << (iph+1) << "] "
        //          << "z in [" << Z_EDGES[iz] << ", " << Z_EDGES[iz+1] << ")  "
        //          << "pt2 in [" << PT2_EDGES[ip] << ", " << PT2_EDGES[ip+1] << ")  "
        //          << "phih in [" << PHI_EDGES[iph] << ", " << PHI_EDGES[iph+1] << (iph==4?"]":")")
        //          << std::endl;
        //DebugPrintSample(T, sel, /*maxPrint=*/5);

        // xB vs Q2
        TString title1 = "x_{B} vs Q^{2}  |  " + suffix;
        c.Clear();
        TH2F* h1 = DrawXBQ2(T, sel, title1);
        c.Print(pdf);
        delete h1;

        // xB vs theta
        TString title2 = "x_{B} vs #theta_{e}  |  " + suffix;
        c.Clear();
        TH2F* h2 = DrawXBTheta(T, sel, title2);
        c.Print(pdf);
        delete h2;
      }
    }
    c.Print(pdf + "]");
  }

  // ---- Optional single-bin PNG export ----
  if (INPUTS::SINGLE_BIN_INDEX >= 1 && INPUTS::SINGLE_BIN_INDEX <= 125) {
    int iz, ip, iph;
    IndexToTriplet125(INPUTS::SINGLE_BIN_INDEX, iz, ip, iph);
    TString selExact = BuildSelCell(iz, ip, iph);
    TString suffix = CellTitleSuffix(iz, ip, iph);

    //std::cout << "[single-bin DEBUG] "
    //          << "z in [" << Z_EDGES[iz] << ", " << Z_EDGES[iz+1] << ")  "
    //          << "pt2 in [" << PT2_EDGES[ip] << ", " << PT2_EDGES[ip+1] << ")  "
    //          << "phih in [" << PHI_EDGES[iph] << ", " << PHI_EDGES[iph+1] << (iph==4?"]":")")
    //          << std::endl;
    //DebugPrintSample(T, selExact, /*maxPrint=*/10);

    TCanvas c("c2", "", 900, 800);
    TString title = "x_{B} vs Q^{2}  |  " + suffix;
    TH2F* h1 = DrawXBQ2(T, selExact, title);
    TString png1 = Form("%s/bin%03d_xbQ2.png", INPUTS::kOutDir.c_str(), INPUTS::SINGLE_BIN_INDEX);
    c.SaveAs(png1);
    delete h1;

    c.Clear();
    TString title2 = "x_{B} vs #theta_{e}  |  " + suffix;
    TH2F* h2 = DrawXBTheta(T, selExact, title2);
    TString png2 = Form("%s/bin%03d_xbTheta.png", INPUTS::kOutDir.c_str(), INPUTS::SINGLE_BIN_INDEX);
    c.SaveAs(png2);
    delete h2;
  }

  std::cout << "Done. PDFs and (optional) PNGs saved in " << INPUTS::kOutDir << std::endl;
  return 0;
}
