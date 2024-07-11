#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>  

#include <iostream>


//compile and run 
//g++ -o comphist compHistoJun34.cpp `root-config --cflags --libs`
//./comphist

void CompareHistograms(const char* file1, const char* file2) {
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
    const std::vector<std::pair<const char*, const char*>> histogramPairs = {
        {"U_Q2_LD2_true", "U_Q2_LD2_simf"},
        {"U_W2_LD2_true", "U_W2_LD2_simf"},
        {"U_nu_LD2_true", "U_nu_LD2_simf"},
        {"U_phih_LD2_true", "U_phih_LD2_simf"},
        {"U_xb_LD2_true", "U_xb_LD2_simf"},
        {"U_y_LD2_true", "U_y_LD2_simf"},
        {"U_z_LD2_true", "U_z_LD2_simf"},
        {"U_targetVz_LD2_true", "U_targetVz_LD2_simf"},
        {"U_pt2_LD2_true", "U_pt2_LD2_simf"}
    };
    const std::map<const char*, std::vector<double>> xStartValues = {
        {"U_Q2_LD2_true", {1.5, 100.0}},
        {"U_W2_LD2_true", {0.0, 100.0}},
        {"U_nu_LD2_true", {100.0, 100.0}}, // No lines for "nu"
        {"U_phih_LD2_true", {100.0, 100.0}}, // No lines for "phih"
        {"U_xb_LD2_true", {100.0, 100.0}}, // No lines for "xb"
        {"U_y_LD2_true", {100.0, 100.0}}, // No lines for "y2"
        {"U_z_LD2_true", {0.3, 0.7}},
        {"U_targetVz_LD2_true", {-7.5, -2.5}},
        {"U_pt2_LD2_true", {-10.0, 1.5}}
    };
    for (const auto& pair : histogramPairs) {
        TCanvas* canvas = new TCanvas(Form("ComparisonCanvas_%s_%s", pair.first, pair.second), Form("Histogram Comparison %s vs %s", pair.first, pair.second), 800, 600);
        
        TH1F* h1 = dynamic_cast<TH1F*>(rootFile1->Get(pair.first));
        TH1F* h2 = dynamic_cast<TH1F*>(rootFile2->Get(pair.second));

        if (!h1 || !h2) {
            std::cerr << "Error: Cannot retrieve histograms from files!" << std::endl;
            return;
        }
        //if (std::string(pair.first) == "U_pt2_C1_true" || std::string(pair.first) == "U_z_C1_true") {
            canvas->SetLogy();
        //}
        h1->SetLineColor(kBlue);
        h2->SetLineColor(kRed);
        h1->SetTitle(Form("Histogram from File 1 (%s)", pair.first));
        h2->SetTitle(Form("Histogram from File 2 (%s)", pair.second));

        float max1 = h1->GetMaximum();
        float max2 = h2->GetMaximum();
        float max = (max1 > max2) ? max1 : max2;
        h1->Scale(max / max1);
        h2->Scale(max / max2);

        h1->Draw("hist");
        h2->Draw("hist, same");

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

        canvas->SaveAs(Form("%s_vs_%s_comparison.pdf", pair.first, pair.second));

        delete canvas;
    }

    delete rootFile1;
    delete rootFile2;
}


void Create2DHistogram(const char* inputFile, const char* outputFile) {
    TFile file(inputFile, "READ");

    TH1F* Q2pre = dynamic_cast<TH1F*>(file.Get("U_Q2_C1_simf"));
    TH1F* xbpre = dynamic_cast<TH1F*>(file.Get("U_xb_C1_simf"));

    if (!Q2pre || !xbpre) {
        // Check if the dynamic_cast was successful
        // Handle the case where the histograms are not TH1D
        std::cerr << "Error: Could not read TH1 histograms from the file." << std::endl;
        file.Close();
        return;
    }

    TH2D* th2d = new TH2D("2D_Histogram", "2D Histogram", xbpre->GetNbinsX(), xbpre->GetXaxis()->GetXmin(), xbpre->GetXaxis()->GetXmax(), Q2pre->GetNbinsX(), Q2pre->GetXaxis()->GetXmin(), Q2pre->GetXaxis()->GetXmax());

    //for (int xbin = 1; xbin <= xbpre->GetNbinsX(); xbin++) {
    //    for (int ybin = 1; ybin <= Q2pre->GetNbinsX(); ybin++) {
    //        double content_xb = xbpre->GetBinContent(xbin);
    //        double content_Q2 = Q2pre->GetBinContent(ybin);
    //        th2d->SetBinContent(xbin, ybin, content_Q2);
    //    }
    //}

    for (int xbin = 1; xbin <=xbpre->GetNbinsX(); xbin++) {
        for (int ybin = 1; ybin <= Q2pre->GetNbinsX(); ybin++) {
            double xValue = xbpre->GetBinContent(xbin);
            double yValue = Q2pre->GetBinContent(ybin);
            th2d->Fill(xValue, yValue);
        }
    }

    TCanvas canvas("canvas", "2D Histogram", 800, 600);
    canvas.cd();
    
    th2d->Draw("colz");
    //TLatex* prelimText = new TLatex();
    //prelimText->SetTextSize(0.08);
    //prelimText->SetTextAngle(45);
    //prelimText->SetTextColorAlpha(kGray + 1, 0.3);  
    //prelimText->SetNDC();
    //prelimText->SetTextAlign(22);  
    //prelimText->DrawLatex(0.5, 0.5, "preliminary");


    canvas.SaveAs(outputFile);

    file.Close();
}



int main() {
    const char* file1 = "/home/matamoros/LD2true.root";
    const char* file2 = "/home/matamoros/LD2sim.root";
    CompareHistograms(file1, file2);
    //Create2DHistogram(file1, "Qpre_vs_xbpre_File1.pdf");
    //Create2DHistogram(file2, "test2D.pdf");

    return 0;
}