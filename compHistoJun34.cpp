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
    std::string file1 = std::string("/home/matamoros/midnovv6") + target + "_test.root";
    std::string file2 = std::string("/home/matamoros/midnovv6") + target + "_sim.root";

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

    const std::vector<std::pair<std::string, std::string>> histogramPairs = {
        {"Q2_" + std::string(target) + "_RGD", "Q2_" + std::string(target) + "_sim"},
        {"W2_" + std::string(target) + "_RGD", "W2_" + std::string(target) + "_sim"},
        {"nu_" + std::string(target) + "_RGD", "nu_" + std::string(target) + "_sim"},
        {"phih_" + std::string(target) + "_RGD", "phih_" + std::string(target) + "_sim"},
        {"xb_" + std::string(target) + "_RGD", "xb_" + std::string(target) + "_sim"},
        {"y_" + std::string(target) + "_RGD", "y_" + std::string(target) + "_sim"},
        {"z_" + std::string(target) + "_RGD", "z_" + std::string(target) + "_sim"},
        {"targetVz_" + std::string(target) + "_RGD", "targetVz_" + std::string(target) + "_sim"},
        {"pt2_" + std::string(target) + "_RGD", "pt2_" + std::string(target) + "_sim"},
        {"ptot_ele_" + std::string(target) + "_RGD", "ptot_ele_" + std::string(target) + "_sim"},
        {"px_ele_" + std::string(target) + "_RGD", "px_ele_" + std::string(target) + "_sim"},
        {"py_ele_" + std::string(target) + "_RGD", "py_ele_" + std::string(target) + "_sim"},
        {"pz_ele_" + std::string(target) + "_RGD", "pz_ele_" + std::string(target) + "_sim"},
        {"E_el" + std::string(target) + "_RGD", "E_el" + std::string(target) + "_sim"},
        {"E_pi" + std::string(target) + "_RGD", "E_pi" + std::string(target) + "_sim"},
        {"theta_el" + std::string(target) + "_RGD", "theta_el" + std::string(target) + "_sim"},
        {"phi_el" + std::string(target) + "_RGD", "phi_el" + std::string(target) + "_sim"},
        {"ptot_pro_" + std::string(target) + "_RGD", "ptot_pro_" + std::string(target) + "_sim"},
        {"px_pro_" + std::string(target) + "_RGD", "px_pro_" + std::string(target) + "_sim"},
        {"py_pro_" + std::string(target) + "_RGD", "py_pro_" + std::string(target) + "_sim"},
        {"pz_pro_" + std::string(target) + "_RGD", "pz_pro_" + std::string(target) + "_sim"},
        {"theta_pi" + std::string(target) + "_RGD", "theta_pi" + std::string(target) + "_sim"},
        {"phi_pi" + std::string(target) + "_RGD", "phi_pi" + std::string(target) + "_sim"},
        {"lu_el"+ std::string(target) + "_RGD", "lu_el" + std::string(target) + "_sim"},
        {"lv_el"+ std::string(target) + "_RGD", "lv_el" + std::string(target) + "_sim"},
        {"lw_el"+ std::string(target) + "_RGD", "lw_el" + std::string(target) + "_sim"},
        {"epcal_el"+ std::string(target) + "_RGD", "epcal_el" + std::string(target) + "_sim"},
        {"Nphe15_"+ std::string(target) + "_RGD", "Nphe15_" + std::string(target) + "_sim"},
        {"Nphe16_"+ std::string(target) + "_RGD", "Nphe16_" + std::string(target) + "_sim"},
        {"chi2_el_"+ std::string(target) + "_RGD", "chi2_el_" + std::string(target) + "_sim"},
        {"chi2_pi_"+ std::string(target) + "_RGD", "chi2_pi_" + std::string(target) + "_sim"},
        {"helicity_"+ std::string(target) + "_RGD", "helicity_" + std::string(target) + "_sim"},
        {"helicity_raw_"+ std::string(target) + "_RGD", "helicity_raw_" + std::string(target) + "_sim"},


    };

    TCanvas* pdfCanvas = new TCanvas("pdfCanvas", "Combined Histogram Comparison", 1000, 800);
    pdfCanvas->Divide(3, 3); 
    int canvasIndex = 1;
    int histCount = 0;

    pdfCanvas->Print("Comparison.pdf["); 

    for (size_t i = 0; i < histogramPairs.size(); ++i) {
        const auto& pair = histogramPairs[i];

        TH1F* h1 = dynamic_cast<TH1F*>(rootFile1->Get(pair.first.c_str()));     // data
        TH1F* h2 = dynamic_cast<TH1F*>(rootFile2->Get(pair.second.c_str()));    // sim

        if (!h1 || !h2) {
            std::cerr << "Error: Cannot retrieve histograms from files!" << std::endl;
            continue;
        }

        // possibility of creating individual canvas for PNG output
        TCanvas* canvas = new TCanvas(Form("ComparisonCanvas_%s_%s", pair.first.c_str(), pair.second.c_str()),
                                      Form("Histogram Comparison %s vs %s", pair.first.c_str(), pair.second.c_str()), 800, 600);

        h1->SetLineColor(kBlue);    // data
        h2->SetLineColor(kRed);     // sim

        // scaling
        double integral1 = h1->Integral();
        double integral2 = h2->Integral();
        if (integral1 != 0) h1->Scale(1.0 / integral1);
        if (integral2 != 0) h2->Scale(1.0 / integral2);

        // getting max value for y-axis scaling
        float max1 = h1->GetMaximum();
        float max2 = h2->GetMaximum();
        float max = (max1 > max2) ? max1 : max2;
        h1->SetMaximum(max * 1.1);
        h1->SetMinimum(0);

        h1->Draw("hist");
        h2->Draw("hist same");

        TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.69);
        legend->AddEntry(h1, "Data", "l");
        legend->AddEntry(h2, "Simulation", "l");
        legend->Draw();

        // save PNG inloop
        canvas->SaveAs(Form("Comparison_%s_vs_%s.png", pair.first.c_str(), pair.second.c_str()));

        // add histogram in PDF 
        pdfCanvas->cd(canvasIndex);
        h1->Draw("hist");
        h2->Draw("hist same");
        legend->Draw();

        canvasIndex++;
        histCount++;

        if (canvasIndex > 9) {  //index to advance in cd(i)
            pdfCanvas->Print("Comparison.pdf"); // page
            pdfCanvas->Clear(); // clear for next page
            pdfCanvas->Divide(3, 3); // redefine 3x3
            canvasIndex = 1; // reset  index
        }

        delete canvas;  //cleanup
    }

    if (histCount % 9 != 0) {
        pdfCanvas->Print("Comparison.pdf"); // print  last page anyway (if less than 9 )
    }

    pdfCanvas->Print("Comparison.pdf]"); 
    delete pdfCanvas; //clean up 

    rootFile1->Close();
    rootFile2->Close();
}

int main() {
    const char* target = "C2"; // can be changed to targets Cu, Sn, LD2, C1, C2
    CompareHistograms(target, {}, {});
    return 0;
}
