// plot_xQ2.cpp
#include <iostream>
#include <string>
#include <memory>

//g++ -O2 plot_xQ2.cpp $(root-config --cflags --libs) -o plot_xQ2
//./plot_xQ2 pass1LD2.root


#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TError.h"

using std::cerr;
using std::cout;
using std::endl;

void setRootStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleSize(0.045, "XYZ");
    gStyle->SetLabelSize(0.04,  "XYZ");
    gStyle->SetPadRightMargin(0.12); // room for Z color bar
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.12);
}

void forceAxesStartAtZero(TH2* h) {
    if (!h) return;
    auto *x = h->GetXaxis();
    auto *y = h->GetYaxis();
    // If current minima are > 0 already, keep them; else set to 0.
    const double xmax = x->GetXmax();
    const double ymax = y->GetXmax();
    double xmin = x->GetXmin();
    double ymin = y->GetXmin();
    if (xmin < 0) xmin = 0;
    if (ymin < 0) ymin = 0;
    // Only apply if 0 is within the original axis limits
    if (xmax > 0) x->SetRangeUser(0., xmax);
    if (ymax > 0) y->SetRangeUser(0., ymax);
}

bool drawAndSave(TH2* h, const std::string& outBase, const std::string& prettyTitle) {
    if (!h) return false;
    forceAxesStartAtZero(h);

    std::unique_ptr<TCanvas> c(new TCanvas(("c_"+outBase).c_str(), prettyTitle.c_str(), 950, 800));
    c->SetTicks(1,1);
    c->SetGrid(0,0);

    h->SetTitle(prettyTitle.c_str());
    h->GetZaxis()->SetTitleOffset(1.1);

    // Prefer linear z by default; switch to logz if very wide range
    double zmin = h->GetMinimum(1); // ignore empty bins
    double zmax = h->GetMaximum();
    bool useLogZ = (zmax > 0 && zmin > 0 && zmax/zmin > 1e3);

    if (useLogZ) c->SetLogz();

    h->Draw("COLZ");

    c->SaveAs((outBase + ".png").c_str());
    c->SaveAs((outBase + ".pdf").c_str());
    return true;
}

int main(int argc, char** argv) {
    setRootStyle();
    gErrorIgnoreLevel = kWarning; // keep output tidy

    std::string inFile = "build/pass1LD2.root";
    if (argc > 1) inFile = argv[1];

    std::unique_ptr<TFile> f(TFile::Open(inFile.c_str(), "READ"));
    if (!f || f->IsZombie()) {
        cerr << "Error: cannot open ROOT file: " << inFile << endl;
        return 1;
    }

    // ROOT ignores the ";1" cycle when using Get()
    const char* H1 = "xQ2_LD2_RGD";
    const char* H2 = "xQ2pos_LD2_RGD";

    TH2* hxq2   = dynamic_cast<TH2*>(f->Get(H1));
    TH2* hxq2ps = dynamic_cast<TH2*>(f->Get(H2));

    if (!hxq2 && !hxq2ps) {
        cerr << "Error: neither '" << H1 << "' nor '" << H2 << "' found in " << inFile << endl;
        f->ls(); // list contents to help debug
        return 2;
    }

    bool ok1 = drawAndSave(hxq2,   "xQ2_LD2_RGD",   "x vs Q^{2} (LD2, RGD)");
    bool ok2 = drawAndSave(hxq2ps, "xQ2pos_LD2_RGD","x vs Q^{2} (positive region, LD2, RGD)");

    // Optional: side-by-side comparison if both exist
    if (ok1 && ok2) {
        std::unique_ptr<TCanvas> c(new TCanvas("c_compare", "xQ2 comparison", 1300, 700));
        c->Divide(2,1);
        c->cd(1);
        if (hxq2) {
            forceAxesStartAtZero(hxq2);
            gPad->SetRightMargin(0.12);
            hxq2->SetTitle("x vs Q^{2} — LD2 (RGD)");
            hxq2->Draw("COLZ");
        }
        c->cd(2);
        if (hxq2ps) {
            forceAxesStartAtZero(hxq2ps);
            gPad->SetRightMargin(0.12);
            hxq2ps->SetTitle("x vs Q^{2} — positive region, LD2 (RGD)");
            hxq2ps->Draw("COLZ");
        }
        c->SaveAs("xQ2_compare_LD2_RGD.png");
        c->SaveAs("xQ2_compare_LD2_RGD.pdf");
    }

    cout << "Done. Wrote: xQ2_LD2_RGD.(png|pdf), xQ2pos_LD2_RGD.(png|pdf)"
            << (ok1 && ok2 ? ", and xQ2_compare_LD2_RGD.(png|pdf)." : ".")
            << endl;

    return 0;
}
