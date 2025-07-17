#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLine.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <cmath>
#include <iostream>
#include <vector>

// Compile & run:
// g++ plot_xQ2_overlay.cpp -o plot_xQ2_overlay `root-config --cflags --glibs`
// ./plot_xQ2_overlay

constexpr double M = 0.9396;
constexpr double E = 10.5;

constexpr double THETA_HI  = 27.0;
constexpr double THETA_MID = 14.7;
constexpr double THETA_LOW = 10.2;

inline double deg2rad(double d) { return d * M_PI / 180.0; }

double Q2_of_x_theta(double x, double thetaDeg) {
    double c = 1.0 - std::cos(deg2rad(thetaDeg));
    double A = 2.0 * M * E * E * c;
    double B = E * c;
    return A * x / (B + M * x);
}

double x_of_Q2_theta(double Q2, double thetaDeg) {
    double c = 1.0 - std::cos(deg2rad(thetaDeg));
    return (E * c * Q2) / (M * (2.0 * E * E * c - Q2));
}

double Q2_of_x_W2(double x, double W2) {
    return (W2 - M * M) * x / (1.0 - x);
}

int main() {
    TFile* fin = TFile::Open("DptinputFiles/REcutC2_test.root", "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open input ROOT file.\n";
        return 1;
    }

    auto* h = dynamic_cast<TH2F*>(fin->Get("xQ2_C2_RGD"));
    if (!h) {
        std::cerr << "ERROR: histogram xQ2_C2_RGD not found.\n";
        return 2;
    }

    // Final xB lines (user defined + two chosen): total 4
    //std::vector<double> xBins = {0.18, 0.32, 0.44, 0.50};  // → 5 xB regions
    std::vector<double> xBins = {0.10, 0.15, 0.20, 0.25, 0.31, 0.44};  // 6 bin edges → 5 xB regions
        std::vector<double> Q2LinesToDraw = {1.0, 5.0};        // keep only Q² = 1 and Q² = 5

    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c", "xQ2 with overlays", 1400, 800);
    c->SetRightMargin(0.15);

    // Create θ curves
    auto make_theta_curve = [](double thetaDeg, const char* name, int color) {
        double c = 1.0 - std::cos(deg2rad(thetaDeg));
        double A = 2.0 * M * E * E * c;
        double B = E * c;
        TF1* f = new TF1(name, "[0]*x/([1]+0.9396*x)", 0.0, 0.7);
        f->SetParameters(A, B);
        f->SetLineWidth(2);
        f->SetLineColor(color);
        f->SetLineStyle(1);
        f->SetRange(0.0, 0.7);
        f->SetMinimum(1.0);
        return f;
    };

    TF1* fThetaHi  = make_theta_curve(THETA_HI , "fThetaHi" , kRed);
    TF1* fThetaMid = make_theta_curve(THETA_MID, "fThetaMid", kRed + 1);
    TF1* fThetaLow = make_theta_curve(THETA_LOW, "fThetaLow", kRed + 2);

    // W² = 4 GeV² line
    TF1* fW2 = new TF1("fW2", "(4 - 0.9396*0.9396)*x/(1 - x)", 0.0, 0.7);
    fW2->SetLineColor(kOrange + 1);
    fW2->SetLineWidth(2);
    fW2->SetLineStyle(1);
    fW2->SetMinimum(1.0);

    // Vertical lines: full height
    std::vector<TLine*> xLines;
    for (double x : xBins) {
        TLine* lx = new TLine(x, 1.0, x, 5.0);
        lx->SetLineColor(kBlue + 2);
        lx->SetLineWidth(2);
        lx->SetLineStyle(1);
        xLines.push_back(lx);
    }

    // Horizontal lines: only Q² = 1 and Q² = 5
    std::vector<TLine*> Q2Lines;
    double x_min = 0.05;
    double x_max = 0.6;
    for (double Q2val : Q2LinesToDraw) {
        TLine* lQ = new TLine(x_min, Q2val, x_max, Q2val);
        lQ->SetLineColor(kGreen + 2);
        lQ->SetLineWidth(2);
        lQ->SetLineStyle(1);
        Q2Lines.push_back(lQ);
    }

    // Draw everything
    h->GetXaxis()->SetTitle("x_{B}");
    h->GetYaxis()->SetTitle("Q^{2}  [GeV^{2}]");
    h->Draw("colz");

    fThetaHi ->Draw("same");
    fThetaMid->Draw("same");
    fThetaLow->Draw("same");
    fW2      ->Draw("same");

    for (auto* l : xLines)  l->Draw("same");
    for (auto* l : Q2Lines) l->Draw("same");

    TLegend* leg = new TLegend(0.70, 0.60, 0.93, 0.87);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(fThetaHi , Form("#theta_{e}=%.1f^{#circ}", THETA_HI ), "l");
    leg->AddEntry(fThetaMid, Form("#theta_{e}=%.1f^{#circ}", THETA_MID), "l");
    leg->AddEntry(fThetaLow, Form("#theta_{e}=%.1f^{#circ}", THETA_LOW), "l");
    leg->AddEntry(fW2, "W^{2} = 4 GeV^{2}", "l");
    leg->Draw();

    c->SaveAs("xQ2_overlay_final.png");

    // Print bin regions
    std::cout << "\n===== Region Definitions (xB between lines) =====\n";
    for (size_t i = 0; i < xBins.size() - 1; ++i) {
        double xb1 = xBins[i];
        double xb2 = xBins[i + 1];
        std::cout << "x_B ∈ [" << xb1 << ", " << xb2 << "]\n";
    }

    fin->Close();
    return 0;
}
