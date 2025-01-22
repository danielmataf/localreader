#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <map>

//compile and run 
//g++ -o comphist compHistoJun34.cpp root-config --cflags --libs
//g++ -o comphist compHistoJun34.cpp $(root-config --cflags --libs)
//./comphist
#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <map>

void CompareHistograms(const char* target, const std::vector<std::string>& plotTitles, const std::vector<std::string>& xTitles) {
    std::string file1 = std::string("/home/matamoros/fullT") + target + "_test.root";
    std::string file2 = std::string("/home/matamoros/fullT") + target + "_sim.root";

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

    //const std::vector<std::pair<std::string, std::string>> histogramPairs = {
    const std::vector<std::tuple<std::string, std::string, std::string>> histogramPairs = {
        //{"histogramNameFromFile1", "histogramNameFromFile2", "X-Axis Title "}
        {"Q2_" + std::string(target) + "_RGD", "Q2_" + std::string(target) + "_sim","temp"},
        {"W2_" + std::string(target) + "_RGD", "W2_" + std::string(target) + "_sim","temp"},
        {"nu_" + std::string(target) + "_RGD", "nu_" + std::string(target) + "_sim","temp"},
        {"phih_" + std::string(target) + "_RGD", "phih_" + std::string(target) + "_sim","temp"},
        {"xb_" + std::string(target) + "_RGD", "xb_" + std::string(target) + "_sim","temp"},
        {"y_" + std::string(target) + "_RGD", "y_" + std::string(target) + "_sim","temp"},
        {"z_" + std::string(target) + "_RGD", "z_" + std::string(target) + "_sim","temp"},
        {"targetVz_" + std::string(target) + "_RGD", "targetVz_" + std::string(target) + "_sim","cm"},
        {"pt2_" + std::string(target) + "_RGD", "pt2_" + std::string(target) + "_sim","temp"},
        {"ptot_ele_" + std::string(target) + "_RGD", "ptot_ele_" + std::string(target) + "_sim","temp"},
        {"px_ele_" + std::string(target) + "_RGD", "px_ele_" + std::string(target) + "_sim","temp"},
        {"py_ele_" + std::string(target) + "_RGD", "py_ele_" + std::string(target) + "_sim","temp"},
        {"pz_ele_" + std::string(target) + "_RGD", "pz_ele_" + std::string(target) + "_sim","temp"},
        {"E_el" + std::string(target) + "_RGD", "E_el" + std::string(target) + "_sim","temp"},
        {"E_pi" + std::string(target) + "_RGD", "E_pi" + std::string(target) + "_sim","temp"},
        {"theta_el" + std::string(target) + "_RGD", "theta_el" + std::string(target) + "_sim","temp"},
        {"phi_el" + std::string(target) + "_RGD", "phi_el" + std::string(target) + "_sim","temp"},
        {"ptot_pro_" + std::string(target) + "_RGD", "ptot_pro_" + std::string(target) + "_sim","temp"},
        {"px_pi_" + std::string(target) + "_RGD", "px_pi_" + std::string(target) + "_sim","temp"},
        {"py_pi_" + std::string(target) + "_RGD", "py_pi_" + std::string(target) + "_sim","temp"},
        {"pz_pi_" + std::string(target) + "_RGD", "pz_pi_" + std::string(target) + "_sim","temp"},
        {"theta_pi" + std::string(target) + "_RGD", "theta_pi" + std::string(target) + "_sim","temp"},
        {"phi_pi" + std::string(target) + "_RGD", "phi_pi" + std::string(target) + "_sim","temp"},
        {"lu_el"+ std::string(target) + "_RGD", "lu_el" + std::string(target) + "_sim","temp"},
        {"lv_el"+ std::string(target) + "_RGD", "lv_el" + std::string(target) + "_sim","temp"},
        {"lw_el"+ std::string(target) + "_RGD", "lw_el" + std::string(target) + "_sim","temp"},
        {"epcal_el"+ std::string(target) + "_RGD", "epcal_el" + std::string(target) + "_sim","temp"},
        {"Nphe15_"+ std::string(target) + "_RGD", "Nphe15_" + std::string(target) + "_sim","temp"},
        {"Nphe16_"+ std::string(target) + "_RGD", "Nphe16_" + std::string(target) + "_sim","temp"},
        {"chi2_el_"+ std::string(target) + "_RGD", "chi2_el_" + std::string(target) + "_sim","temp"},
        {"chi2_pi_"+ std::string(target) + "_RGD", "chi2_pi_" + std::string(target) + "_sim","temp"},
        {"helicity_"+ std::string(target) + "_RGD", "helicity_" + std::string(target) + "_sim","temp"},
        {"helicity_raw_"+ std::string(target) + "_RGD", "helicity_raw_" + std::string(target) + "_sim","temp"},
        {"targetVz_pi_"+std::string(target) + "_RGD","targetVz_pi_" + std::string(target) + "_sim","cm"},


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

    // Canvas for individual PNG output
    TCanvas* canvas = new TCanvas(Form("ComparisonCanvas_%s_%s", hist1Name.c_str(), hist2Name.c_str()),
                                  Form("Histogram Comparison %s vs %s", hist1Name.c_str(), hist2Name.c_str()), 800, 600);

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

    //TLegend* legend = new TLegend(0.7, 0.75, 0.9, 0.85);  //too high
    TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.69);   //good size and lvl

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
    const char* target = "C2"; // can be changed to targets Cu, Sn, LD2, C1, C2
    CompareHistograms(target, {}, {});
    return 0;
}
