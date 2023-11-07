#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <iostream>
#include <map>
#include <vector>

void CompareHistograms(const char* file1, const char* file2) {
    //open root files + error handling
    TFile* rootFile1 = new TFile(file1, "READ");
    if (!rootFile1->IsOpen()) {
        std::cerr << "Error: Cannot open file1!" << std::endl;
        return;
    }
    TFile* rootFile2 = new TFile(file2, "READ");
    if (!rootFile2->IsOpen()) {
        std::cerr << "Error: Cannot open file2!" << std::endl;
        return;
    }

    // histogram names to recover 
    const std::vector<std::pair<const char*, const char*>> histogramPairs = {
        {"Q2pre", "Q2pre"},
        {"W2pre", "W2pre"},
        {"nupre", "nupre"},
        {"phihpre", "phihpre"},
        {"xbpre", "xbpre"},
        {"ypre", "ypre"},
        {"zpre", "zpre"},
        {"Vzpre", "Vzpre"},
        {"pt2pre", "pt2pre"}
    };

    //some info for plots
    //( name , { cutlow, cuthigh, title, titleXaxis, titleYaxis, rangeXmin,rangeXmax, bool4logscale})
    const std::map<const char*, std::tuple<double, double, std::string, std::string, std::string, double, double, bool>> xStartValues = {
        {"Q2pre", {1.5, 100.0, "Q^{2} Distribution in LD2", "Q^{2}(Gev^{2})", "", 0.0, 6.0, false}},
        {"W2pre", {0.0, 100.0, "W^{2} Distribution in LD2", "W^{2}", "", 0.0, 20.0, false}},
        {"nupre", {100.0, 100.0, "#nu Distribution in LD2", "#nu", "", 0.0, 12.0, false}},
        {"phihpre", {100.0, 100.0, "#phi_{h} Distribution in LD2", "#phi_{h}", "", 0.0, 360.0, false}},
        {"xbpre", {100.0, 100.0, "x_{b} Distribution in LD2", " x_{b}", "", 0.0, 1.0, false}},
        {"ypre", {0.25, 0.85, "y Distribution in LD2", "y (GeV)", "", 0.2, 1.0, false}},
        {"zpre", {0.3, 0.7, "z Distribution in LD2", " z", "", 0.0, 1.0, false}},
        {"Vzpre", {-9.0, -1.5, "V_{z} Distribution in LD2", "V_{z} (cm) ", "", -20.0, 10.0, false}},
        {"pt2pre", {100.0, 100.0, "p_{T}^{2} Distribution in LD2", "p_{T}^{2} (Gev^{2})", "", 0.0, 3.0, true}}
    };

    //create canvas 4 each variable comp
    for (const auto& pair : histogramPairs) {
        TCanvas* canvas = new TCanvas(Form("ComparisonCanvas_%s_%s", pair.first, pair.second), Form("Histogram Comparison %s vs %s", pair.first, pair.second), 800, 600);

        //recover hists
        TH1F* h1 = dynamic_cast<TH1F*>(rootFile1->Get(pair.first));
        TH1F* h2 = dynamic_cast<TH1F*>(rootFile2->Get(pair.second));

        //error handler
        if (!h1 || !h2) {
            std::cerr << "Error: Cannot retrieve histograms from files!" << std::endl;
            return;
        }

        //set histogram specifications
        h1->SetLineColor(kRed);
        h2->SetLineColor(kBlue);
        h1->SetTitle(Form("%s ", std::get<2>(xStartValues.at(pair.first)).c_str()));
        h2->SetTitle(Form("%s ", std::get<2>(xStartValues.at(pair.second)).c_str()));

        //axis titles
        h1->GetXaxis()->SetTitle(std::get<3>(xStartValues.at(pair.first)).c_str());
        h1->GetYaxis()->SetTitle(std::get<4>(xStartValues.at(pair.first)).c_str());
        h2->GetXaxis()->SetTitle(std::get<3>(xStartValues.at(pair.second)).c_str());
        h2->GetYaxis()->SetTitle(std::get<4>(xStartValues.at(pair.second)).c_str());

        //x-axis range
        const double xMin = std::get<5>(xStartValues.at(pair.first));
        const double xMax = std::get<6>(xStartValues.at(pair.first));
        h1->GetXaxis()->SetRangeUser(xMin, xMax);
        h2->GetXaxis()->SetRangeUser(xMin, xMax);

        //log scale if needed
        if (std::get<7>(xStartValues.at(pair.first))) {
            canvas->SetLogy();
        }

        //scaling manually histogram to the same maximum value
        float max1 = h1->GetMaximum();
        float max2 = h2->GetMaximum();
        float max = (max1 > max2) ? max1 : max2;
        h1->Scale(max / max1);
        h2->Scale(max / max2);

        //drawing
        h1->Draw("hist");
        h2->Draw("hist, same");

        //creating and drawing Cuts
        const double xStart1 = std::get<0>(xStartValues.at(pair.first));
        const double xStart2 = std::get<1>(xStartValues.at(pair.second));
        if (xStart1 != 100.0 && xStart2 != 100.0) { // Skip the exaggerated values
            TLine* line1 = new TLine(xStart1, 0, xStart1, max);
            line1->SetLineColor(kBlack);
            line1->SetLineWidth(2);
            line1->SetLineStyle(2); //dashed
            line1->Draw("");
            TLine* line2 = new TLine(xStart2, 0, xStart2, max);
            line2->SetLineColor(kBlack);
            line2->SetLineWidth(2);
            line2->SetLineStyle(2); //dashed
            line2->Draw("");
        }

        // Save the canvas as a PDF file
        canvas->SaveAs(Form("%s_vs_%s_comparison.pdf", pair.first, pair.second));
        delete canvas;
    }

    delete rootFile1;
    delete rootFile2;
}

int main() {
    // add two files to compare; suggestion file1 = MC, file2 = REC
    const char* file1 = "/home/matamoros/Desktop/reboot/CleanCode/build/chop_WINsimLD2.root";
    const char* file2 = "/home/matamoros/Desktop/reboot/CleanCode/build/chop_WINrecLD2.root";

    CompareHistograms(file1, file2);

    return 0;
}
