// electronCountsR.cpp
// Plot only-electron counts vs xB for all four targets (CxC, Cu, LD2, Sn)
// One page per theta bin, points only with Poisson errors.
//
// Normalized per-theta-bin distributions (sum_x N = 1).
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

void usage(const char* prog){
  std::cerr
    << "Usage: " << prog << " --scheme DAT|BIS [--png]\n"
    << "Reads R_unfold_<SCHEME>_<TARGET>.root (TARGET = CxC, Cu, LD2, Sn)\n"
    << "and unf2D_<TARGET>outw4.root with h_xB_thetael_UNF.\n"
    << "Plots N(xB,theta) (normalized) with Poisson errors, one page per theta bin.\n"
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
  const std::vector<int> colors  = {kBlack, kGreen+2, kBlue, kRed};

  // --- Original data (R_unfold_...) ---
  std::vector<TargetData> data;
  data.reserve(targets.size());

  for (size_t i=0;i<targets.size();++i){
    TargetData t;
    t.name  = targets[i];
    t.color = colors[i];

    std::string fname = "R_unfold_" + cfg.scheme + "_" + t.name + ".root";
    t.file.reset(TFile::Open(fname.c_str(), "READ"));
    if (!t.file || t.file->IsZombie()){
      std::cerr << "[Err] error opening file " << fname << "\n";
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

  // --- Unfolded data (unf2D_...outw4.root, h_xB_thetael_UNF) ---
  std::vector<TargetData> dataUnf;
  dataUnf.reserve(targets.size());

for (size_t i=0;i<targets.size();++i){
  TargetData t;
  t.name  = targets[i];
  t.color = colors[i]; // same color as original

  // Directory where the UNF files live
  const std::string unfDir = "/home/matamoros/Downloads/RooUnfold/";

  std::string fname = unfDir + "unf2D_" + t.name + "outw4.root";
  t.file.reset(TFile::Open(fname.c_str(), "READ"));
  if (!t.file || t.file->IsZombie()){
    std::cerr << "[Err] error opening UNF file " << fname << "\n";
    return 20;
  }

  std::string hname = "h_xB_thetael_UNF";
  t.hist = dynamic_cast<TH2*>(t.file->Get(hname.c_str()));
  if (!t.hist){
    std::cerr << "[Err] Missing UNF histogram " << hname << " in " << fname << "\n";
    return 21;
  }

  dataUnf.push_back(std::move(t));
}

  // --- Use binning of first original hist for page structure ---
  TH2* href = data.front().hist;
  TAxis* ax = href->GetXaxis();
  TAxis* ay = href->GetYaxis();
  int ny = ay->GetNbins();   // number of theta bins (pages)

  std::string outName = "electroncount_" + cfg.scheme + (cfg.save_png ? "" : ".pdf");

  bool firstPage = true;
  for (int iy=1; iy<=ny; ++iy){
    auto* c = new TCanvas(Form("c_theta_%d", iy), "", 900, 650);

    double ymin = 1e99, ymax = -1e99;

    std::vector<TGraphErrors*> graphs;     // original
    std::vector<TGraphErrors*> graphsUnf;  // unfolded

    // ---------- ORIGINAL DATA GRAPHS ----------
    for (auto& t : data) {
      TGraphErrors* g = new TGraphErrors();
      const int nx = ax->GetNbins();

      double total = 0.0;
      for (int ix=1; ix<=nx; ++ix)
        total += t.hist->GetBinContent(ix, iy);

      g->Set(nx);

      for (int ix=1; ix<=nx; ++ix) {
        const double x    = ax->GetBinCenter(ix);
        const double Nraw = t.hist->GetBinContent(ix, iy);
        const double eRaw = std::sqrt(std::max(0.0, Nraw));

        const double N  = (total > 0.0) ? (Nraw / total) : 0.0;
        const double eN = (total > 0.0) ? (eRaw / total) : 0.0;

        g->SetPoint(ix-1, x, N);
        g->SetPointError(ix-1, 0.0, eN);

        ymin = std::min(ymin, N - eN);
        ymax = std::max(ymax, N + eN);
      }

      g->SetMarkerColor(t.color);
      g->SetLineColor(t.color);
      g->SetMarkerStyle(20);   // solid marker (data)
      g->SetMarkerSize(1.0);
      g->SetLineWidth(1);

      graphs.push_back(g);
    }

    // ---------- UNFOLDED DATA GRAPHS ----------
    // Same theta-bin index 'iy', but use their own xB binning.
    for (auto& t : dataUnf) {
      TGraphErrors* g = new TGraphErrors();
      TAxis* axUnf = t.hist->GetXaxis();
      const int nxUnf = axUnf->GetNbins();

      double total = 0.0;
      for (int ix=1; ix<=nxUnf; ++ix)
        total += t.hist->GetBinContent(ix, iy);

      g->Set(nxUnf);

      for (int ix=1; ix<=nxUnf; ++ix) {
        const double x    = axUnf->GetBinCenter(ix);      // UNF xB bin center (5 bins)
        const double Nraw = t.hist->GetBinContent(ix, iy);
        const double eRaw = std::sqrt(std::max(0.0, Nraw));

        const double N  = (total > 0.0) ? (Nraw / total) : 0.0;
        const double eN = (total > 0.0) ? (eRaw / total) : 0.0;

        g->SetPoint(ix-1, x, N);
        g->SetPointError(ix-1, 0.0, eN);

        ymin = std::min(ymin, N - eN);
        ymax = std::max(ymax, N + eN);
      }

      g->SetMarkerColor(t.color);
      g->SetLineColor(t.color);
      g->SetMarkerStyle(24);   // hollow marker (unfolded)
      g->SetMarkerSize(1.0);
      g->SetLineWidth(1);

      graphsUnf.push_back(g);
    }

    // ---------- FRAME & AXES ----------
    double x1 = 0.0; // force X-axis to start at 0, as you requested
    double x2 = ax->GetBinUpEdge(ax->GetNbins());
    double pad = (ymax > 0.0) ? 0.15 * ymax : 0.1;
    double y1 = std::max(0.0, ymin - pad);
    double y2 = ymax + pad;

    auto* frame = c->DrawFrame(x1, y1, x2, y2);
    frame->GetXaxis()->SetTitle("x_{B}");
    frame->GetYaxis()->SetTitle("normalized electron counts");

    double thL = ay->GetBinLowEdge(iy);
    double thH = ay->GetBinUpEdge(iy);
    frame->SetTitle(Form("#theta_{e} #in [%.3g, %.3g] deg", thL, thH));

    // Draw original first, then UNF over them
    for (auto* g : graphs)    g->Draw("PE SAME");
    for (auto* g : graphsUnf) g->Draw("PE SAME");

    // ---------- LEGEND ----------
    auto* leg = new TLegend(0.60, 0.62, 0.88, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.032);

    // assumes same ordering of targets in graphs and graphsUnf
    leg->AddEntry(graphs[0],    "CxC (data)", "p");
    leg->AddEntry(graphsUnf[0], "CxC (unf.)", "p");
    leg->AddEntry(graphs[1],    "Cu (data)",  "p");
    leg->AddEntry(graphsUnf[1], "Cu (unf.)",  "p");
    leg->AddEntry(graphs[2],    "LD2 (data)", "p");
    leg->AddEntry(graphsUnf[2], "LD2 (unf.)", "p");
    leg->AddEntry(graphs[3],    "Sn (data)",  "p");
    leg->AddEntry(graphsUnf[3], "Sn (unf.)",  "p");
    leg->Draw();

    // ---------- OUTPUT ----------
    if (cfg.save_png) {
      c->Print(Form("electroncount_%s_theta%d.png", cfg.scheme.c_str(), iy));
    } else {
      if (firstPage) {
        c->Print((outName + "(").c_str());
        firstPage = false;
      }
      else if (iy==ny) {
        c->Print((outName + ")").c_str());
      }
      else {
        c->Print(outName.c_str());
      }
    }

    // cleanup
    for (auto* g : graphs)    delete g;
    for (auto* g : graphsUnf) delete g;
    delete c;
  }

  if (cfg.save_png)
    std::cout << "[bump] saved PNGs electroncount_" << cfg.scheme << "_theta*.png\n";
  else
    std::cout << "[bump] wrote " << outName << " with " << ny << " pages.\n";

  return 0;
}
