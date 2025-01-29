#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <map>
// How to compile and run
// g++ -o threecomp threecomp.cpp $(root-config --cflags --libs)
// ./threecomp

void Compare3Histograms(const char* target) {
    std::string file1 = std::string("/home/matamoros/RGDv7") + target + "_test.root";
    std::string file2 = std::string("/home/matamoros/jan") + target + "_sim.root";
    std::string file3 = std::string("/home/matamoros/symm") + target + "_sim.root";

    TFile* rootFile1 = new TFile(file1.c_str(), "READ");
    if (!rootFile1->IsOpen()) {
        std::cerr << "Error: Cannot open file1!" << std::endl;
        return;
    }

    TFile* rootFile2 = new TFile(file2.c_str(), "READ");
    if (!rootFile2->IsOpen()) {
        std::cerr << "Error: Cannot open file2!" << std::endl;
        return;
    }

    TFile* rootFile3 = new TFile(file3.c_str(), "READ");
    if (!rootFile3->IsOpen()) {
        std::cerr << "Error: Cannot open file3!" << std::endl;
        return;
    }

    const std::vector<std::string> histogramNames = {
        "Q2_", "W2_", "nu_", "phih_", "xb_", "y_", "z_", "targetVz_", "pt2_",
        "ptot_ele_", "px_ele_", "py_ele_", "pz_ele_", "E_el", "E_pi", "theta_el", "phi_el",
        "ptot_pro_", "px_pi_", "py_pi_", "pz_pi_", "theta_pi", "phi_pi", "lu_el", "lv_el",
        "lw_el", "epcal_el", "Nphe15_", "Nphe16_", "chi2_el_", "chi2_pi_", "helicity_",
        "helicity_raw_"
    };

    const std::vector<std::string> xTitles = {
        "Q^{2} [GeV^{2}]", "W^{2} [GeV^{2}]", "#nu [GeV]", "#phi_{h} [deg]", "x_{B}", "y", "z", 
        "Target V_{z} [cm]", "p_{T}^{2} [GeV^{2}]", "p_{tot} Electron [GeV]", "p_{x} Electron [GeV]",
        "p_{y} Electron [GeV]", "p_{z} Electron [GeV]", "E Electron [GeV]", "E Pion [GeV]", 
        "#theta Electron [deg]", "#phi Electron [deg]", "p_{tot} Pion [GeV]", "p_{x} Pion [GeV]",
        "p_{y} Pion [GeV]", "p_{z} Pion [GeV]", "#theta Pion [deg]", "#phi Pion [deg]", 
        "lu Electron", "lv Electron", "lw Electron", "E_{PCAL} Electron [GeV]", "Nphe15", "Nphe16", 
        "#chi^{2} Electron", "#chi^{2} Pion", "Helicity", "Helicity (Raw)"
    };

    TCanvas* pdfCanvas = new TCanvas("pdfCanvas", "Combined Histogram Comparison", 1000, 800);
    pdfCanvas->Divide(3, 3);
    int canvasIndex = 1;
    int histCount = 0;

    pdfCanvas->Print("Comparison.pdf[");

    for (size_t i = 0; i < histogramNames.size(); ++i) {
        std::string histName1 = histogramNames[i] + std::string(target) + "_RGD";
        std::string histName2 = histogramNames[i] + std::string(target) + "_sim";
        std::string histName3 = histogramNames[i] + std::string(target) + "_sim";

        TH1F* h1 = dynamic_cast<TH1F*>(rootFile1->Get(histName1.c_str()));
        TH1F* h2 = dynamic_cast<TH1F*>(rootFile2->Get(histName2.c_str()));
        TH1F* h3 = dynamic_cast<TH1F*>(rootFile3->Get(histName3.c_str()));

        if (!h1 || !h2 || !h3) {
            std::cerr << "Error retrieving histograms from files" << std::endl;
            continue;
        }

        TCanvas* canvas = new TCanvas(Form("ComparisonCanvas_%s", histName1.c_str()),
                                      Form("Histogram Comparison %s", histName1.c_str()), 800, 600);

        h1->SetLineColor(kBlue);
        h2->SetLineColor(kRed);
        h3->SetLineColor(kGreen);

        double integral1 = h1->Integral();
        double integral2 = h2->Integral();
        double integral3 = h3->Integral();

        if (integral1 != 0) h1->Scale(1.0 / integral1);
        if (integral2 != 0) h2->Scale(1.0 / integral2);
        if (integral3 != 0) h3->Scale(1.0 / integral3);

        float max = std::max({h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum()});
        h1->SetMaximum(max * 1.1);
        h1->SetMinimum(0);

        h1->GetXaxis()->SetTitle(xTitles[i].c_str());
        h2->GetXaxis()->SetTitle(xTitles[i].c_str());
        h3->GetXaxis()->SetTitle(xTitles[i].c_str());

        h1->Draw("hist");
        h2->Draw("hist same");
        h3->Draw("hist same");

        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(h1, "RGD Data 0v7", "l");
        legend->AddEntry(h2, "FullTorus", "l");
        legend->AddEntry(h3, "SymmTorus", "l");
        legend->Draw();

        canvas->SaveAs(Form("Comparison_%s.png", histName1.c_str()));

        pdfCanvas->cd(canvasIndex);
        h1->Draw("hist");
        h2->Draw("hist same");
        h3->Draw("hist same");
        legend->Draw();

        canvasIndex++;
        histCount++;

        if (canvasIndex > 9) {
            pdfCanvas->Print("Comparison.pdf");
            pdfCanvas->Clear();
            pdfCanvas->Divide(3, 3);
            canvasIndex = 1;
        }

        delete canvas;
    }

    if (histCount % 9 != 0) {
        pdfCanvas->Print("Comparison.pdf");
    }

    pdfCanvas->Print("Comparison.pdf]");
    delete pdfCanvas;

    rootFile1->Close();
    rootFile2->Close();
    rootFile3->Close();
}

int main() {
    const char* target = "C2";
    Compare3Histograms(target);
    return 0;
}
