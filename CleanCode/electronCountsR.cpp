// electronCountsR.cpp
// ------------------------------*
// Plots electron counts vs xB for all 4 targets, per theta bin.
// Reads from R_unfold_TARGET.root files.
// Outputs multi-page PDF (one page per theta bin).
//how to compile 
// g++ electronCountsR.cpp -o electronCountsR $(root-config --cflags --libs)
//# DAT binning (default)
//./electronCountsR --scheme DAT

//# BIS binning
//./electronCountsR --scheme BIS

//# if needed :save PNGs for each θ page
//./electronCountsR --scheme DAT --png

// electronCountsR.cpp
// Plot only-electron counts vs xB for all four targets (CxC, Cu, LD2, Sn)
// One page per theta bin, points only with Poisson errors.
//
// not normalized!!!! maybe consider normalizing with intgral 1 . Errors are small 
// Build: g++ electronCountsR.cpp -o electronCountsR $(root-config --cflags --libs)
// Run  : ./electronCountsR --scheme DAT
//        ./electronCountsR --scheme BIS
//        ./electronCountsR --scheme DAT --png

#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TStyle.h>

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <cctype>

struct Config {
  std::string scheme = "DAT"; // DAT or BIS
  bool save_png = false;      // PDF by default, PNG if true
} cfg;
// whatever options we could have either BIS or DATA 
void usage(const char* prog){
  std::cerr
    << "Usage: " << prog << " --scheme DAT|BIS [--png]\n"
    << "Reads unfold_<SCHEME>_<TARGET>.root (TARGET = CxC, Cu, LD2, Sn)\n"
    << "Plots N(xB,theta) with Poisson errors, one page per theta bin.\n"
    << "Outputs electroncount_<SCHEME>.pdf or PNGs.\n";
}

void parseArgs(int argc, char** argv){
  for (int i=1;i<argc;++i){
    std::string a = argv[i];
    if (a=="--scheme" && i+1<argc) cfg.scheme = argv[++i];
    else if (a=="--png") cfg.save_png = true;
    else if (a=="-h" || a=="--help"){ usage(argv[0]); std::exit(0); }
    else { std::cerr << "[WARN] Unknown option: " << a << "\n"; usage(argv[0]); std::exit(1); }
  }
  for (auto& c : cfg.scheme) c = std::toupper(c);
  if (cfg.scheme!="DAT" && cfg.scheme!="BIS"){
    std::cerr << "[FATAL] --scheme must be DAT or BIS\n"; std::exit(2);
  }
}

struct TargetData {
  std::string name;
  std::unique_ptr<TFile> file;
  TH2* hist = nullptr;
  int color;
};

int main(int argc, char** argv){
  parseArgs(argc, argv);

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  const std::vector<std::string> targets = {"CxC","Cu","LD2","Sn"};
  const std::vector<int> colors = {kBlack, kGreen+2, kBlue, kRed};

  std::vector<TargetData> data;
  data.reserve(targets.size());

  //get hists
  for (size_t i=0;i<targets.size();++i){
    TargetData t;
    t.name = targets[i];
    t.color = colors[i];

    std::string fname = "R_unfold_" + cfg.scheme + "_" + t.name + ".root";
    t.file.reset(TFile::Open(fname.c_str(), "READ"));
    if (!t.file || t.file->IsZombie()){
      std::cerr << "[Err] error openeing file " << fname << "\n";
      return 10;
    }

    std::string hname = "h" + cfg.scheme + "_onlye_" + t.name;
    t.hist = dynamic_cast<TH2*>(t.file->Get(hname.c_str()));
    if (!t.hist){
      std::cerr << "[Err] Missing histogram " << hname << " in " << fname << "\n";
      return 11;
    }

    data.push_back(std::move(t));
  }

  // use binning of first hist
  TH2* href = data.front().hist;
  TAxis* ax = href->GetXaxis();
  TAxis* ay = href->GetYaxis();
  int ny = ay->GetNbins();

  std::string outName = "electroncount_" + cfg.scheme + (cfg.save_png ? "" : ".pdf");

  bool firstPage = true;
  for (int iy=1; iy<=ny; ++iy){
    auto* c = new TCanvas(Form("c_theta_%d", iy), "", 900, 650);
    //c->SetGrid(); //no GRID

    double ymin = 1e99, ymax = -1e99;

    std::vector<TGraphErrors*> graphs;
  //double ymin = +1e99, ymax = -1e99;

for (auto& t : data) {
  TGraphErrors* g = new TGraphErrors();
  const int nx = ax->GetNbins();

  // 
  double total = 0.0;
  for (int ix=1; ix<=nx; ++ix)   total += t.hist->GetBinContent(ix, iy); // this should sum all the counts each xb in this theata bin

    g->Set(nx);

  for (int ix=1; ix<=nx; ++ix) {  //looping a second time over the xb bins to fill graphs now that we have the total count per bin info 
    const double x    = ax->GetBinCenter(ix);   //standard xB bin center
    const double Nraw = t.hist->GetBinContent(ix, iy);  // "raw" counts in this (xB, theta) bin 
    const double eRaw = std::sqrt(std::max(0.0, Nraw));  // Poisson raw error, so just the sqrt

    const double N    = (total > 0.0) ? (Nraw / total) : 0.0;   // total is known cause we looped twice // acts as an "average" count
    const double eN   = (total > 0.0) ? (eRaw / total) : 0.0;   // total is known cause we looped twice// acts as an "average" count

    g->SetPoint(ix-1, x, N);
    g->SetPointError(ix-1, 0.0, eN);

    ymin = std::min(ymin, N - eN);
    ymax = std::max(ymax, N + eN);
  }

  g->SetMarkerColor(t.color);
  g->SetLineColor(t.color);   //same color
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.0);
  g->SetLineWidth(1);

  graphs.push_back(g);
}

    double x1 = ax->GetBinLowEdge(1);
    double x2 = ax->GetBinUpEdge(ax->GetNbins());
    double pad = 0.15 * ymax;
    double y1 = std::max(0.0, ymin - pad);
    double y2 = ymax + pad;

    auto* frame = c->DrawFrame(x1, y1, x2, y2);
    frame->GetXaxis()->SetTitle("x_{B}");
    frame->GetYaxis()->SetTitle("normalized electron counts");

    double thL = ay->GetBinLowEdge(iy);
    double thH = ay->GetBinUpEdge(iy);
    frame->SetTitle(Form("#theta_{e} #in [%.3g, %.3g] deg", thL, thH));

    for (auto* g : graphs) g->Draw("PE SAME");

    auto* leg = new TLegend(0.65, 0.72, 0.88, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
    leg->AddEntry(graphs[0], "CxC", "p");
    leg->AddEntry(graphs[1], "Cu", "p");
    leg->AddEntry(graphs[2], "LD2", "p");
    leg->AddEntry(graphs[3], "Sn", "p");
    leg->Draw();

    // output
    if (cfg.save_png) {
      c->Print(Form("electroncount_%s_theta%d.png", cfg.scheme.c_str(), iy));
    } else {
      if (firstPage) c->Print((outName + "(").c_str()), firstPage=false;
      else if (iy==ny) c->Print((outName + ")").c_str());
      else c->Print(outName.c_str());
    }

    for (auto* g : graphs) delete g;
    delete c;
  }

  if (cfg.save_png)
    std::cout << "[bump] saved PNGs electroncount_" << cfg.scheme << "_theta*.png\n";
  else
    std::cout << "[bump] wrote " << outName << " with " << ny << " pages.\n";

  return 0;
}
