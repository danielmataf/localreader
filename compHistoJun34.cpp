#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <tuple>

void CompareHistograms(const char* target) {
    std::string file1 = std::string("/home/matamoros/jan") + target + "_test.root";
    std::string file2 = std::string("/home/matamoros/jan") + target + "_sim.root";

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

    // Mapping histogram names to x-axis labels
    const std::vector<std::tuple<std::string, std::string, std::string>> histogramPairs = {
        {"Q2_" + std::string(target) + "_RGD", "Q2_" + std::string(target) + "_sim", "Q^{2} [GeV^{2}]"},
        {"W2_" + std::string(target) + "_RGD", "W2_" + std::string(target) + "_sim", "W^{2} [GeV^{2}]"},
        {"nu_" + std::string(target) + "_RGD", "nu_" + std::string(target) + "_sim", "#nu [GeV]"},
        {"phih_" + std::string(target) + "_RGD", "phih_" + std::string(target) + "_sim", "#phi_{h} [deg]"},
        {"xb_" + std::string(target) + "_RGD", "xb_" + std::string(target) + "_sim", "x_{B}"},
        {"y_" + std::string(target) + "_RGD", "y_" + std::string(target) + "_sim", "y"},
        {"z_" + std::string(target) + "_RGD", "z_" + std::string(target) + "_sim", "z"},
        {"targetVz_" + std::string(target) + "_RGD", "targetVz_" + std::string(target) + "_sim", "Target V_{z} [cm]"},
        {"pt2_" + std::string(target) + "_RGD", "pt2_" + std::string(target) + "_sim", "p_{T}^{2} [GeV^{2}]"},
        {"ptot_ele_" + std::string(target) + "_RGD", "ptot_ele_" + std::string(target) + "_sim", "p_{tot} Electron [GeV]"},
        {"px_ele_" + std::string(target) + "_RGD", "px_ele_" + std::string(target) + "_sim", "p_{x} Electron [GeV]"},
        {"py_ele_" + std::string(target) + "_RGD", "py_ele_" + std::string(target) + "_sim", "p_{y} Electron [GeV]"},
        {"pz_ele_" + std::string(target) + "_RGD", "pz_ele_" + std::string(target) + "_sim", "p_{z} Electron [GeV]"},
        {"E_el" + std::string(target) + "_RGD", "E_el" + std::string(target) + "_sim", "E Electron [GeV]"},
        {"E_pi" + std::string(target) + "_RGD", "E_pi" + std::string(target) + "_sim", "E Pion [GeV]"},
        {"theta_el" + std::string(target) + "_RGD", "theta_el" + std::string(target) + "_sim", "#theta Electron [deg]"},
        {"phi_el" + std::string(target) + "_RGD", "phi_el" + std::string(target) + "_sim", "#phi Electron [deg]"},
        {"ptot_pro_" + std::string(target) + "_RGD", "ptot_pro_" + std::string(target) + "_sim", "p_{tot} Pion [GeV]"},
        {"px_pi_" + std::string(target) + "_RGD", "px_pi_" + std::string(target) + "_sim", "p_{x} Pion [GeV]"},
        {"py_pi_" + std::string(target) + "_RGD", "py_pi_" + std::string(target) + "_sim", "p_{y} Pion [GeV]"},
        {"pz_pi_" + std::string(target) + "_RGD", "pz_pi_" + std::string(target) + "_sim", "p_{z} Pion [GeV]"},
        {"theta_pi" + std::string(target) + "_RGD", "theta_pi" + std::string(target) + "_sim", "#theta Pion [deg]"},
        {"phi_pi" + std::string(target) + "_RGD", "phi_pi_" + std::string(target) + "_sim", "#phi Pion [deg]"},
        {"chi2_el_" + std::string(target) + "_RGD", "chi2_el_" + std::string(target) + "_sim", "#chi^{2} Electron"},
        {"chi2_pi_" + std::string(target) + "_RGD", "chi2_pi_" + std::string(target) + "_sim", "#chi^{2} Pion"},
        {"helicity_" + std::string(target) + "_RGD", "helicity_" + std::string(target) + "_sim", "Helicity"},
        {"helicity_raw_" + std::string(target) + "_RGD", "helicity_raw_" + std::string(target) + "_sim", "Helicity (Raw)"}
    };

    TCanvas* pdfCanvas = new TCanvas("pdfCanvas", "Combined Histogram Comparison", 1000, 800);
    pdfCanvas->Divide(3, 3);
    int canvasIndex = 1;
    int histCount = 0;

    pdfCanvas->Print("Comparison.pdf[");

    for (const auto& [hist1Name, hist2Name, xAxisTitle] : histogramPairs) {
        TH1F* h1 = dynamic_cast<TH1F*>(rootFile1->Get(hist1Name.c_str()));
        TH1F* h2 = dynamic_cast<TH1F*>(rootFile2->Get(hist2Name.c_str()));

        if (!h1 || !h2) {
            std::cerr << "Error: Cannot retrieve histograms from files!" << std::endl;
            continue;
        }

        TCanvas* canvas = new TCanvas(Form("ComparisonCanvas_%s", hist1Name.c_str()),
                                      Form("Histogram Comparison %s", hist1Name.c_str()), 800, 600);

        h1->SetLineColor(kBlue);
        h2->SetLineColor(kRed);

        double integral1 = h1->Integral();
        double integral2 = h2->Integral();
        if (integral1 != 0) h1->Scale(1.0 / integral1);
        if (integral2 != 0) h2->Scale(1.0 / integral2);

        float max = std::max(h1->GetMaximum(), h2->GetMaximum());
        h1->SetMaximum(max * 1.1);
        h1->SetMinimum(0);

        h1->GetXaxis()->SetTitle(xAxisTitle.c_str());
        h2->GetXaxis()->SetTitle(xAxisTitle.c_str());

        h1->Draw("hist");
        h2->Draw("hist same");

        TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.69);
        legend->AddEntry(h1, "Data", "l");
        legend->AddEntry(h2, "Simulation", "l");
        legend->Draw();

        canvas->SaveAs(Form("Comparison_%s_vs_%s.png", hist1Name.c_str(), hist2Name.c_str()));

        pdfCanvas->cd(canvasIndex);
        h1->Draw("hist");
        h2->Draw("hist same");
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
}

int main() {
    const char* target = "C2"; 
    CompareHistograms(target);
    return 0;
}
