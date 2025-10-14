#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <algorithm>     // for std::max
#include <iostream>
#include <vector>
#include <tuple>

// g++ compHistoJun34.cpp $(root-config --cflags --libs) -o compHistoJun34
// ./compHistoJun34

void CompareHistograms(const char* target) {
    // fix: add trailing slash in paths
    std::string file1 = std::string("/home/matamoros/pass1") + target + ".root";
    std::string file2 = std::string("/home/matamoros/sep") + target + "_sim.root";

    TFile* rootFile1 = new TFile(file1.c_str(), "READ");
    if (!rootFile1->IsOpen()) { std::cerr << "Error: Cannot open file1!\n"; return; }
    TFile* rootFile2 = new TFile(file2.c_str(), "READ");
    if (!rootFile2->IsOpen()) { std::cerr << "Error: Cannot open file2!\n"; return; }

    // histogram name pairs + X titles
    const std::vector<std::tuple<std::string, std::string, std::string>> histogramPairs = {
        {"Q2_"      + std::string(target) + "_RGD", "Q2_"      + std::string(target) + "_sim", "Q^{2} [GeV^{2}]"},
        {"W2_"      + std::string(target) + "_RGD", "W2_"      + std::string(target) + "_sim", "W^{2} [GeV^{2}]"},
        {"nu_"      + std::string(target) + "_RGD", "nu_"      + std::string(target) + "_sim", "#nu [GeV]"},
        {"phih_"    + std::string(target) + "_RGD", "phih_"    + std::string(target) + "_sim", "#phi_{h} [deg]"},
        {"xb_"      + std::string(target) + "_RGD", "xb_"      + std::string(target) + "_sim", "x_{B}"},
        {"y_"       + std::string(target) + "_RGD", "y_"       + std::string(target) + "_sim", "y"},
        {"z_"       + std::string(target) + "_RGD", "z_"       + std::string(target) + "_sim", "z"},
        {"targetVz_"+ std::string(target) + "_RGD", "targetVz_"+ std::string(target) + "_sim", "Target V_{z} [cm]"},
        {"pt2_"     + std::string(target) + "_RGD", "pt2_"     + std::string(target) + "_sim", "p_{T}^{2} [GeV^{2}]"}
    };

    // global style: no titles/stats/fit box
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TCanvas* pdfCanvas = new TCanvas("pdfCanvas", "Combined Histogram Comparison", 1000, 800);
    pdfCanvas->Divide(3, 3);
    int canvasIndex = 1, histCount = 0;

    pdfCanvas->Print("Comparison.pdf[");

    for (const auto& [hist1Name, hist2Name, xAxisTitle] : histogramPairs) {
        TH1F* h1 = dynamic_cast<TH1F*>(rootFile1->Get(hist1Name.c_str())); // data (blue)
        TH1F* h2 = dynamic_cast<TH1F*>(rootFile2->Get(hist2Name.c_str())); // sim  (red)
        if (!h1 || !h2) { std::cerr << "Error: Cannot retrieve " << hist1Name << " / " << hist2Name << "\n"; continue; }

        // per-hist style
        h1->SetStats(0); h2->SetStats(0);
        h1->SetTitle(""); h2->SetTitle("");
        h1->SetLineColor(kBlue); h2->SetLineColor(kRed);
        h1->SetLineWidth(2);     h2->SetLineWidth(2);

        // normalize to unit area (shape comparison)
        double i1 = h1->Integral();
        double i2 = h2->Integral();
        if (i1) h1->Scale(1.0 / i1);
        if (i2) h2->Scale(1.0 / i2);

        // axis cosmetics (match your other script)
        h1->GetXaxis()->SetTitle(xAxisTitle.c_str());
        h1->GetXaxis()->SetTitleOffset(0.8);  // try 0.9 → 0.8 depending on taste
        h1->GetXaxis()->SetTitleSize(0.055);
        h1->GetXaxis()->SetLabelSize(0.045);
        h1->GetYaxis()->SetLabelSize(0.045);
        h1->GetYaxis()->SetTitleSize(0.050);  // optional: uncomment the next line if you want a Y title
        // h1->GetYaxis()->SetTitle("Normalized counts");

        // y-range
        float ymax = std::max(h1->GetMaximum(), h2->GetMaximum());
        h1->SetMaximum(ymax * 1.1);
        h1->SetMinimum(0);

        // canvas for PNG (blank window title)
        TCanvas* canvas = new TCanvas(Form("ComparisonCanvas_%s", hist1Name.c_str()), "", 800, 600);
        h1->Draw("hist");
        h2->Draw("hist same");

        // legend (same format/pos/size as before)
        TLegend* leg = new TLegend(0.67, 0.70, 0.92, 0.88);
        leg->SetTextSize(0.05);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h1, "Data", "l");
        leg->AddEntry(h2, "Simulation", "l");
        leg->Draw();

        canvas->SaveAs(Form("Comparison_%s_vs_%s.png", hist1Name.c_str(), hist2Name.c_str()));

        // draw on PDF pad
        pdfCanvas->cd(canvasIndex);
        h1->Draw("hist");
        h2->Draw("hist same");
        leg->DrawClone();  // clone for the PDF pad

        canvasIndex++;
        histCount++;
        if (canvasIndex > 9) {
            pdfCanvas->Print("Comparison.pdf");
            pdfCanvas->Clear();
            pdfCanvas->Divide(3, 3);
            canvasIndex = 1;
        }

        delete leg;
        delete canvas;
    }

    if (histCount % 9 != 0) pdfCanvas->Print("Comparison.pdf");
    pdfCanvas->Print("Comparison.pdf]");
    delete pdfCanvas;

    rootFile1->Close();
    rootFile2->Close();
}

int main() {
    const char* target = "LD2";
    CompareHistograms(target);
    return 0;
}

