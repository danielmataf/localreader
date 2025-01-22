#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <map>
//how to compile and run 
//g++ -o threecomp threecomp.cpp $(root-config --cflags --libs)
//./threecomp

void Compare3Histograms(const char* target, const std::vector<std::string>& plotTitles, const std::vector<std::string>& xTitles) {
    std::string file1 = std::string("/home/matamoros/ful") + target + "_test.root";
    std::string file2 = std::string("/home/matamoros/ful") + target + "_sim.root";
    std::string file3 = std::string("/home/matamoros/other") + target + "_sim.root";

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

    TCanvas* pdfCanvas = new TCanvas("pdfCanvas", "Combined Histogram Comparison", 1000, 800);
    pdfCanvas->Divide(3, 3);
    int canvasIndex = 1;
    int histCount = 0;

    pdfCanvas->Print("Comparison.pdf[");

    for (const auto& name : histogramNames) {
        std::string histName1 = name + std::string(target) + "_RGD";
        std::string histName2 = name + std::string(target) + "_sim";
        std::string histName3 = name + std::string(target) + "_sim";

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

        float max1 = h1->GetMaximum();
        float max2 = h2->GetMaximum();
        float max3 = h3->GetMaximum();
        float max = std::max({max1, max2, max3});

        h1->SetMaximum(max * 1.1);
        h1->SetMinimum(0);

        h1->Draw("hist");
        h2->Draw("hist same");
        h3->Draw("hist same");

        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(h1, "RGD Data", "l");
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
    Compare3Histograms(target, {}, {});
    return 0;
}
