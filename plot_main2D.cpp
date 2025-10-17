// plot_main2D.cpp
// Usage examples:
//   root -l -b -q 'plot_main2D.cpp("mon_LD2.root","LD2")'
// or compile:
//   g++ -O2 -std=c++17 plot_main2D.cpp `root-config --cflags --libs` -o plot_main2D
//   ./plot_main2D mon_LD2.root LD2
//   ./plot_main2D CleanCode/build/pass1LD2.root LD2_RGD

#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TClass.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TString.h"
#include "TError.h"

namespace {
bool DrawAndSaveHistogram(TObject* obj, const std::string& outPng) {
    if (!obj) { 
        std::cerr << "[ERR] Null object, cannot draw: " << outPng << "\n";
        return false;
    }

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    TH2* h2 = dynamic_cast<TH2*>(obj);
    TH1* h1 = dynamic_cast<TH1*>(obj);
    if (!h1) {
        std::cerr << "[ERR] Object is not a TH1/TH2: " << obj->GetName() << "\n";
        return false;
    }

    // Common: force X to start at 0
    double xmax = h1->GetXaxis()->GetXmax();
    if (xmax < 0) xmax = 0;                    // safety
    h1->GetXaxis()->SetRangeUser(0.0, xmax);

    // Canvas (wider for 2D due to palette)
    TCanvas c("c","c", h2 ? 900 : 800, h2 ? 750 : 600);

    if (h2) {
        // For 2D: also force Y to start at 0
        double ymax = h2->GetYaxis()->GetXmax();
        if (ymax < 0) ymax = 0;
        h2->GetYaxis()->SetRangeUser(0.0, ymax);

        c.SetRightMargin(0.14);
        h2->Draw("COLZ");
    } else {
        // For 1D: start Y (counts) at 0
        h1->SetMinimum(0.0);
        h1->Draw("HIST");
    }

    c.SaveAs(outPng.c_str());
    std::cout << "[OK] Wrote " << outPng << "\n";
    return true;
}

//getters
template <typename HType>
HType* GetHistChecked(TFile* f, const std::string& name) {
    if (!f || f->IsZombie()) return nullptr;
    TObject* obj = f->Get(name.c_str());
    if (!obj) {
        std::cerr << "[ERR] Not found in file: " << name << "\n";
        return nullptr;
    }
    HType* h = dynamic_cast<HType*>(obj);
    if (!h) {
        std::cerr << "[ERR] Object exists but is not of expected type for: " << name << "\n";
        return nullptr;
    }
    return h;
}

} // namespace

/** Main plotting function you asked for. */
void PlotMain2D(const std::string& rootFilePath, const std::string& targetName) {
    // The saved keys are "xQ2_<targetName>" and "pt2z_<targetName>"
    const std::string key_xQ2  = "xQ2_"  + targetName;
    const std::string key_pt2z = "pt2z_" + targetName;

    std::unique_ptr<TFile> fin(TFile::Open(rootFilePath.c_str(), "READ"));
    if (!fin || fin->IsZombie()) {
        std::cerr << "[FATAL] Cannot open file: " << rootFilePath << "\n";
        return;
    }

    TH2* hxQ2  = GetHistChecked<TH2>(fin.get(), key_xQ2);
    TH2* hpt2z = GetHistChecked<TH2>(fin.get(), key_pt2z);

    if (!hxQ2 || !hpt2z) {
        std::cerr << "[ABORT] Missing required histograms. Check names & targetName.\n";
        return;
    }

    // Save with the exact filenames you specified
    DrawAndSaveHistogram(hxQ2,  "manuscript_main2DxQ.png");
    DrawAndSaveHistogram(hpt2z, "manuscript_main2Dzpt.png");
}

/** Optional: a tiny CLI for convenience. */
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.root> <targetName>\n"
                  << "Example: " << argv[0] << " mon_LD2.root LD2\n";
        return 1;
    }
    PlotMain2D(argv[1], argv[2]);
    return 0;
}
