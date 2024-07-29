#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <map>

//compile and run 
//g++ -o comphist compHistoJun34.cpp root-config --cflags --libs
//./comphist

void CompareHistograms(const char* target, const std::vector<std::string>& plotTitles, const std::vector<std::string>& xTitles) {
    std::string file1 = std::string("/home/matamoros/pos") + target + "_data.root";
    std::string file2 = std::string("/home/matamoros/pos") + target + "_sim.root";

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
        {"Q2_" + std::string(target) + "_true", "Q2_" + std::string(target) + "_sim"},
        {"W2_" + std::string(target) + "_true", "W2_" + std::string(target) + "_sim"},
        {"nu_" + std::string(target) + "_true", "nu_" + std::string(target) + "_sim"},
        {"phih_" + std::string(target) + "_true", "phih_" + std::string(target) + "_sim"},
        {"xb_" + std::string(target) + "_true", "xb_" + std::string(target) + "_sim"},
        {"y_" + std::string(target) + "_true", "y_" + std::string(target) + "_sim"},
        {"z_" + std::string(target) + "_true", "z_" + std::string(target) + "_sim"},
        {"targetVz_" + std::string(target) + "_true", "targetVz_" + std::string(target) + "_sim"},
        {"pt2_" + std::string(target) + "_true", "pt2_" + std::string(target) + "_sim"},
    };

    std::map<std::string, int> targetMap = {
        {"LD2", 1},
        {"Sn", 2},
        {"C1", 3},
        {"Cu", 4},
        {"C2", 5}
    };

    double lowvz, highvz;

    switch (targetMap[target]) {
        case 1: // LD2
            lowvz = -7.5;
            highvz = -2.5;
            break;
        case 2: // Sn
        case 3: // C1
            lowvz = -3.5;
            highvz = -1.5;
            break;
        case 4: // Cu
        case 5: // C2
            lowvz = -8.5;
            highvz = -6.5;
            break;
        default:
            std::cerr << "Error: Unknown target!" << std::endl;
            return;
    }

    const std::map<std::string, std::vector<double>> xStartValues = {
        {"Q2_" + std::string(target) + "_true", {1.5, 100.0}},
        {"W2_" + std::string(target) + "_true", {0.0, 100.0}},
        {"nu_" + std::string(target) + "_true", {100.0, 100.0}}, // No lines for "nu"
        {"phih_" + std::string(target) + "_true", {100.0, 100.0}}, // No lines for "phih"
        {"xb_" + std::string(target) + "_true", {100.0, 100.0}}, // No lines for "xb"
        {"y_" + std::string(target) + "_true", {100.0, 100.0}}, // No lines for "y2"
        {"z_" + std::string(target) + "_true", {0.3, 0.7}},
        {"targetVz_" + std::string(target) + "_true", {lowvz, highvz}},
        {"pt2_" + std::string(target) + "_true", {-10.0, 1.5}}
    };

    for (size_t i = 0; i < histogramPairs.size(); ++i) {
        const auto& pair = histogramPairs[i];
        TCanvas* canvas = new TCanvas(Form("ComparisonCanvas_%s_%s", pair.first.c_str(), pair.second.c_str()), Form("Histogram Comparison %s vs %s", pair.first.c_str(), pair.second.c_str()), 800, 600);
        
        TH1F* h1 = dynamic_cast<TH1F*>(rootFile1->Get(pair.first.c_str()));     //data
        TH1F* h2 = dynamic_cast<TH1F*>(rootFile2->Get(pair.second.c_str()));    //sim

        if (!h1 || !h2) {
            std::cerr << "Error: Cannot retrieve histograms from files!" << std::endl;
            return;
        }
        if (pair.first == "pt2_" + std::string(target) + "_data" ){ // || pair.first == "U_z_" + std::string(target) + "_data") {
            canvas->SetLogy();
        }
        h1->SetLineColor(kBlue);    //data
        h2->SetLineColor(kRed);     //sim

        std::string plotTitle = plotTitles[i];
        std::string xTitle = xTitles[i];
        //std::string yTitle = yTitles[i];
        
        if (pair.first == "pt2_" + std::string(target) + "_true") {
            plotTitle = "Comparison of p_{t}^{2}";
            xTitle = "p_{t}^{2} (GeV^{2})";
        }
        else if (pair.first == "Q2_" + std::string(target) + "_true") {
            xTitle = "Q^{2} (GeV^{2})";
        }
        else if (pair.first == "W2_" + std::string(target) + "_true") {
            xTitle = "W^{2} (GeV^{2})";
        }
        else if (pair.first == "xb_" + std::string(target) + "_true") {
            xTitle = "x_{B}";
        }
        else if (pair.first == "y_" + std::string(target) + "_true") {
            xTitle = "y";
        }
        else if (pair.first == "z_" + std::string(target) + "_true") {
            xTitle = "z";
        }
        else if (pair.first == "targetVz_" + std::string(target) + "_true") {
            xTitle = "V_{z} (cm)";
        }
        else if (pair.first == "phih_" + std::string(target) + "_true") {
            xTitle = "\\phi_{h} (deg)";
        }

        h1->SetTitle(plotTitle.c_str());
        h1->GetXaxis()->SetTitle(xTitle.c_str());
        //h1->GetYaxis()->SetTitle(yTitle.c_str());

        //4 scaling 
        double integral1 = h1->Integral();
        double integral2 = h2->Integral();
        if (integral1 != 0) h1->Scale(1.0 / integral1);
        if (integral2 != 0) h2->Scale(1.0 / integral2);

        //get max post scaling 2 properly plot with y-axis
        float max1 = h1->GetMaximum();
        float max2 = h2->GetMaximum();
        float max = (max1 > max2) ? max1 : max2;

        //y-axis range
        h1->SetMaximum(max * 1.1);  // multiplying max to leave some space in the axis
        h1->SetMinimum(0);          // start from 0

        h1->Draw("hist");
        h2->Draw("hist same");

        TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.69);  
        legend->SetFillColorAlpha(0, 0.5);  // see through
        legend->SetTextSize(0.03);  
        legend->AddEntry(h1, "Data", "l");
        legend->AddEntry(h2, "Simulation", "l");
        legend->Draw();

        if (xStartValues.find(pair.first) != xStartValues.end()) {
            const std::vector<double>& xStarts = xStartValues.at(pair.first);
            for (const double xStart : xStarts) {
                if (xStart != 100.0) { 
                    TLine* line = new TLine(xStart, 0, xStart, max);
                    line->SetLineColor(kBlack);
                    line->SetLineWidth(2);
                    line->SetLineStyle(2); 
                    line->Draw("");
                }
            }
        }

        TLatex* prelimText = new TLatex();
        prelimText->SetTextSize(0.08); 
        prelimText->SetTextAngle(45);
        prelimText->SetTextColorAlpha(kGray + 1, 0.3);  // alpha is for transparency
        prelimText->SetNDC();
        prelimText->SetTextAlign(22);  // centered
        prelimText->DrawLatex(0.5, 0.5, "preliminary");

        canvas->SaveAs(Form("%s_vs_%s_comparison.pdf", pair.first.c_str(), pair.second.c_str()));

        delete legend;
        delete canvas;
    }

    delete rootFile1;
    delete rootFile2;
}

int main() {
    const char* target = "Sn";  //  "Cu", "LD2", "C1", or "C2" 
    
    // Titles for plots
    std::vector<std::string> plotTitles = {
        "Comparison of Q^{2}", "Comparison of W^{2}", "Comparison of \\nu", 
        "Comparison of \\phi_{h}", "Comparison of x_{B}", "Comparison of y", 
        "Comparison of z", "Comparison of V_{z}", "Comparison of p_{t}^{2}"
    };

    // X-axis titles
    std::vector<std::string> xTitles = {
        "Q^{2} (GeV^{2})", "W^{2} (GeV^{2})", "\\nu (GeV)", 
        "\\phi_{h} (deg)", "x_{B}", "y", 
        "z", "V_{z} (cm)", "p_{t}^{2} (GeV^{2})"
    };

    // Y-axis titles

    CompareHistograms(target, plotTitles, xTitles);

    return 0;
}