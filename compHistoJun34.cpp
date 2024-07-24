#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <iostream>
#include <vector>
#include <map>

//compile and run 
//g++ -o comphist compHistoJun34.cpp `root-config --cflags --libs`
//./comphist

void CompareHistograms(const char* target) {
    std::string file1 = std::string("/home/matamoros/") + target + "_true.root";
    std::string file2 = std::string("/home/matamoros/") + target + "_sim.root";

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
        {"pt2_" + std::string(target) + "_true", "pt2_" + std::string(target) + "_sim"}
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

    for (const auto& pair : histogramPairs) {
        TCanvas* canvas = new TCanvas(Form("ComparisonCanvas_%s_%s", pair.first.c_str(), pair.second.c_str()), Form("Histogram Comparison %s vs %s", pair.first.c_str(), pair.second.c_str()), 800, 600);
        
        TH1F* h1 = dynamic_cast<TH1F*>(rootFile1->Get(pair.first.c_str()));
        TH1F* h2 = dynamic_cast<TH1F*>(rootFile2->Get(pair.second.c_str()));

        if (!h1 || !h2) {
            std::cerr << "Error: Cannot retrieve histograms from files!" << std::endl;
            return;
        }
        //if (pair.first == "U_pt2_" + std::string(target) + "_true" || pair.first == "U_z_" + std::string(target) + "_true") {
            canvas->SetLogy();
        //}
        h1->SetLineColor(kBlue);
        h2->SetLineColor(kRed);
        h1->SetTitle(Form("Histogram from File 1 (%s)", pair.first.c_str()));
        h2->SetTitle(Form("Histogram from File 2 (%s)", pair.second.c_str()));

        float max1 = h1->GetMaximum();
        float max2 = h2->GetMaximum();
        float max = (max1 > max2) ? max1 : max2;
        h1->Scale(max / max1);
        h2->Scale(max / max2);

        h1->Draw("hist");
        h2->Draw("hist same");

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

        delete canvas;
    }

    delete rootFile1;
    delete rootFile2;
}

void Create2DHistogram(const char* inputFile, const char* outputFile) {
    TFile file(inputFile, "READ");

    TH1F* Q2pre = dynamic_cast<TH1F*>(file.Get("U_Q2_Sn_simf"));
    TH1F* xbpre = dynamic_cast<TH1F*>(file.Get("U_xb_Sn_simf"));

    if (!Q2pre || !xbpre) {
        std::cerr << "Error: Could not read TH1 histograms from the file." << std::endl;
        file.Close();
        return;
    }

    TH2D* th2d = new TH2D("2D_Histogram", "2D Histogram", xbpre->GetNbinsX(), xbpre->GetXaxis()->GetXmin(), xbpre->GetXaxis()->GetXmax(), Q2pre->GetNbinsX(), Q2pre->GetXaxis()->GetXmin(), Q2pre->GetXaxis()->GetXmax());

    for (int xbin = 1; xbin <= xbpre->GetNbinsX(); xbin++) {
        for (int ybin = 1; ybin <= Q2pre->GetNbinsX(); ybin++) {
            double xValue = xbpre->GetBinContent(xbin);
            double yValue = Q2pre->GetBinContent(ybin);
            th2d->Fill(xValue, yValue);
        }
    }

    TCanvas canvas("canvas", "2D Histogram", 800, 600);
    canvas.cd();
    
    th2d->Draw("colz");

    canvas.SaveAs(outputFile);

    file.Close();
}

int main() {
    const char* target = "C1";  //  "Cu", "LD2", "C1", or "C2" 
    CompareHistograms(target);
    //Create2DHistogram(file1, "Qpre_vs_xbpre_File1.pdf");
    //Create2DHistogram(file2, "test2D.pdf");

    return 0;
}
