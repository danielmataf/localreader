// make4unfold.cpp
//----------------------------------
// Create 125 xb–theta 2D histograms per binning (DAT, BIS) from hadron tree,
// organized in 5 directories (phi_1..phi_5). Also add ONE top-level 2D histo
// for electrons-only per binning. Target hardcoded below.
// Compile: g++ -O2 -std=c++17 make4unfold.cpp `root-config --cflags --libs` -o make4unfold
// Run:     ./make4unfold
//--------------------------------------
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <string>

//pre sets for inputs and ouputs; define  target HERE  
static const std::string kInFile   = "build/Trees_CxC_RGD.root";
static const std::string kTarget   = "CxC";
static const std::string kTreeHad  = "tEv_had_CxC_RGD";
static const std::string kTreeOnly = "tEv_onlye_CxC_RGD";

static const double Z_EDGES[6]   = {0.30, 0.38, 0.46, 0.54, 0.62, 0.70};
static const double PT2_EDGES[6] = {0.00, 0.24, 0.48, 0.72, 0.96, 1.20};
static const double PHI_EDGES[6] = {0.0, 72.0, 144.0, 216.0, 288.0, 360.0};

//DAT binning    //custom binning defined by UNFOLDING PROCEDURE
static const double xEdgesDAT[]  = {0.10, 0.11, 0.13, 0.16, 0.20, 0.25, 0.36, 1.00};
static const int    NX_DAT       = sizeof(xEdgesDAT)/sizeof(double) - 1;
static const double thEdgesDAT[] = {7.6,  8.4,  10.3, 23.0};
static const int    NY_DAT       = sizeof(thEdgesDAT)/sizeof(double) - 1;

//BIS binning
static const double xEdgesBIS[]  = {0.10, 0.11, 0.15, 0.19, 0.29, 1.00};
static const int    NX_BIS       = sizeof(xEdgesBIS)/sizeof(double) - 1;
static const double thEdgesBIS[] = {7.6,  8.8,  11.0, 23.0};
static const int    NY_BIS       = sizeof(thEdgesBIS)/sizeof(double) - 1;


static TString cellCut(int iz, int ip, int iph) {
    //triplet selection range value
  const bool lastPhi = (iph == 4);
  return Form("z>=%.6g && z<%.6g && pt2>=%.6g && pt2<%.6g && %s",
              Z_EDGES[iz], Z_EDGES[iz+1],
              PT2_EDGES[ip], PT2_EDGES[ip+1],
              lastPhi ? Form("phih>=%.6g && phih<=%.6g", PHI_EDGES[iph], PHI_EDGES[iph+1])
                      : Form("phih>=%.6g && phih<%.6g",  PHI_EDGES[iph], PHI_EDGES[iph+1]));
}

static void fill25PerPhi(TTree* Thad,TFile* fout, const char* subdirName, const double* xEdges, int nx, const double* thEdges, int ny, const char* prefix) {
  // Make /phi_k directory  and cd into it so Draw creates hists there
  TDirectory* dphi = fout->mkdir(subdirName);
  dphi->cd();
    //will give 5 directories in the root file
  // Build 5x5 (z,pt2) cells
  for (int ip = 0; ip < 5; ++ip) {
    for (int iz = 0; iz < 5; ++iz) {

      //give a unique name per cell
      TString hname = Form("%s_z%d_pt%d", prefix, iz+1, ip+1);
      TH2F* h = new TH2F(hname, hname, nx, xEdges, ny, thEdges);
      h->GetXaxis()->SetTitle("x_{B}");
      h->GetYaxis()->SetTitle("#theta_{e}  [deg]");

      //assign histo to this directory
      h->SetDirectory(dphi);

      //cell selection
      TString sel = cellCut(iz, ip, atoi(subdirName+4) - 1); // "phi_%d" -> take index
      //fill with draw via Draw 
      Thad->Draw(TString::Format("th:xb>>%s", hname.Data()), sel, "goff");
      dphi->WriteTObject(h);
      delete h; // keep memory clean; key remains in file
    }
  }

  //exits directory and enters the next phih directory
  fout->cd();
}

//additional function to make only one histo for only_e counts 
static void writeOnlyE(TTree* Te, TFile* fout, const double* xEdges, int nx, const double* thEdges, int ny, const char* name) {
  //make sure we are not in any directory aka top lvl
  fout->cd();
  TH2F* h = new TH2F(name, name, nx, xEdges, ny, thEdges);
  h->GetXaxis()->SetTitle("x_{B}");
  h->GetYaxis()->SetTitle("#theta_{e}  [deg]");
  h->SetDirectory(fout);

  //normal fill w/o selection
  Te->Draw(TString::Format("th:xb>>%s", name), "", "goff");
  fout->WriteTObject(h);
  delete h;
}

// ###################################

int main() {
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE); // added this to prevent auto-attachment slip ups

  TFile* fin = TFile::Open(kInFile.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Cannot open input file: " << kInFile << std::endl;
    return 1;
  }
    //sanity checks in trees 
  TTree* Thad = (TTree*)fin->Get(kTreeHad.c_str());
  TTree* Te   = (TTree*)fin->Get(kTreeOnly.c_str());
  if (!Thad) { std::cerr << "Missing hadron tree: " << kTreeHad << std::endl; return 1; }
  if (!Te)   { std::cerr << "Missing only-e tree: " << kTreeOnly << std::endl; return 1; }

  //DAT file
  {
    TString out = Form("unfold_DAT_%s.root", kTarget.c_str());
    TFile* fout = TFile::Open(out, "RECREATE");
    if (!fout || fout->IsZombie()) { std::cerr << "Cannot create " << out << std::endl; return 1; }

    //5 phi directories
    for (int iph = 0; iph < 5; ++iph) {
      TString sub = Form("phi_%d", iph+1);
      //setting the prefix with binning tag, target, and phi index (optional)
      TString prefix = Form("hDAT_%s_phi%d", kTarget.c_str(), iph+1);
      //Important: pass iph via directory name to selection builder
      //We set current directory inside fill25PerPhi
      fill25PerPhi(Thad, fout, sub.Data(), xEdgesDAT, NX_DAT, thEdgesDAT, NY_DAT, prefix.Data());
    }

    //then go back to top-lvl to make  onlye hist
    writeOnlyE(Te, fout, xEdgesDAT, NX_DAT, thEdgesDAT, NY_DAT,
               Form("hDAT_onlye_%s", kTarget.c_str()));
    fout->Close();
    delete fout;
    std::cout << "Wrote: " << out << std::endl;
  }

  //BIS file
  {
    //same principels as above
    TString out = Form("unfold_BIS_%s.root", kTarget.c_str());
    TFile* fout = TFile::Open(out, "RECREATE");
    if (!fout || fout->IsZombie()) { std::cerr << "unabl 2 create " << out << std::endl; return 1; }

    for (int iph = 0; iph < 5; ++iph) {
      TString sub = Form("phi_%d", iph+1);
      TString prefix = Form("hBIS_%s_phi%d", kTarget.c_str(), iph+1);
      fill25PerPhi(Thad, fout, sub.Data(), xEdgesBIS, NX_BIS, thEdgesBIS, NY_BIS, prefix.Data());
    }

    writeOnlyE(Te, fout, xEdgesBIS, NX_BIS, thEdgesBIS, NY_BIS,
               Form("hBIS_onlye_%s", kTarget.c_str()));

    fout->Close();
    delete fout;
    std::cout << "Wrote: " << out << std::endl;
  }

  fin->Close();
  delete fin;
  return 0;
}
