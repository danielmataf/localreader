#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <iostream>

void CompareHistograms(const char* file1, const char* file2) {
    // Open the first ROOT file
    TFile* rootFile1 = new TFile(file1, "READ");
    if (!rootFile1->IsOpen()) {
        std::cerr << "Error: Cannot open file1!" << std::endl;
        return;
    }

    // Open the second ROOT file
    TFile* rootFile2 = new TFile(file2, "READ");
    if (!rootFile2->IsOpen()) {
        std::cerr << "Error: Cannot open file2!" << std::endl;
        return;
    }

    // Define the pairs of histogram names you want to compare
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

    // Define the x-values for TLines for each histogram
    const std::map<const char*, std::vector<double>> xStartValues = {
        {"Q2pre", {1.5, 100.0}},
        {"W2pre", {0.0, 100.0}},
        {"nupre", {100.0, 100.0}}, // No lines for "nu"
        {"phihpre", {100.0, 100.0}}, // No lines for "phih"
        {"xbpre", {100.0, 100.0}}, // No lines for "xb"
        {"ypre", {100.0, 100.0}}, // No lines for "y2"
        {"zpre", {0.3, 0.7}},
        {"Vzpre", {-7.5, -2.5}},
        {"pt2pre", {100.0, 100.0}}
    };

    // Create a canvas for each pair of histograms
    for (const auto& pair : histogramPairs) {
        TCanvas* canvas = new TCanvas(Form("ComparisonCanvas_%s_%s", pair.first, pair.second), Form("Histogram Comparison %s vs %s", pair.first, pair.second), 800, 600);

        // Retrieve the histograms from both files using dynamic_cast
        TH1F* h1 = dynamic_cast<TH1F*>(rootFile1->Get(pair.first));
        TH1F* h2 = dynamic_cast<TH1F*>(rootFile2->Get(pair.second));

        if (!h1 || !h2) {
            std::cerr << "Error: Cannot retrieve histograms from files!" << std::endl;
            return;
        }

        // Set the histogram colors and titles
        h1->SetLineColor(kRed);
        h2->SetLineColor(kBlue);
        h1->SetTitle(Form("Histogram from File 1 (%s)", pair.first));
        h2->SetTitle(Form("Histogram from File 2 (%s)", pair.second));

        //scaling manually histo to same maximum value
        float max1 = h1->GetMaximum();
        float max2 = h2->GetMaximum();
        float max = (max1 > max2) ? max1 : max2;
        h1->Scale(max / max1);
        h2->Scale(max / max2);

        // Draw the histograms on the canvas
        h1->Draw("hist");
        h2->Draw("hist, same");

        // Draw vertical lines (if required)
        if (xStartValues.find(pair.first) != xStartValues.end()) {
            const std::vector<double>& xStarts = xStartValues.at(pair.first);
            for (const double xStart : xStarts) {
                if (xStart != 100.0) { // Skip the exaggerated values
                    TLine* line = new TLine(xStart, 0, xStart, max);
                    line->SetLineColor(kBlack);
                    line->SetLineWidth(2);
                    line->SetLineStyle(2); // Set the line style to dashed
                    line->Draw("");
                }
            }
        }

        // Save the canvas as a PDF file
        canvas->SaveAs(Form("%s_vs_%s_comparison.pdf", pair.first, pair.second));

        delete canvas;
    }

    delete rootFile1;
    delete rootFile2;
}
/*
void Create2DHistogram(const char* file, const char* outputName) {
    TFile* rootFile = new TFile(file, "READ");
    if (!rootFile->IsOpen()) {
        std::cerr << "Error: Can't open file" << std::endl;
        return;
    }

    // Retrieve the histograms using dynamic_cast
    TH1F* h_Qpre = dynamic_cast<TH1F*>(rootFile->Get("xbpre"));
    TH1F* h_xbpre = dynamic_cast<TH1F*>(rootFile->Get("Q2pre"));

    if (!h_Qpre || !h_xbpre) {
        std::cerr << "Error: Can't retrieve histograms from file" << std::endl;
        return;
    }

    // Create a 2D histogram for Qpre vs xbpre
    const int nBinsX = h_xbpre->GetNbinsX();
    const int nBinsY = h_Qpre->GetNbinsX();
    TH2F* h_Qxb = new TH2F("Qpre_vs_xbpre", "Qpre vs xbpre", nBinsX, h_xbpre->GetBinLowEdge(1), h_xbpre->GetBinLowEdge(nBinsX + 1), nBinsY, h_Qpre->GetBinLowEdge(1), h_Qpre->GetBinLowEdge(nBinsY + 1));
    //TH2F* h_Qxb = new TH2F("Qpre_vs_xbpre", "Qpre vs xbpre", 100, 0, 1,100, 0, 6);

    for (int i = 1; i <= nBinsX; i++) {
        for (int j = 1; j <= nBinsY; j++) {
            h_Qxb->SetBinContent(i, j, h_Qpre->GetBinContent(j) * h_xbpre->GetBinContent(i));
        }
    }

    // Create a canvas for the 2D histogram
    TCanvas* canvas = new TCanvas(outputName, outputName, 800, 600);
    canvas->cd();
    h_Qxb->Draw("colz");

    // Save the 2D histogram as a PDF file
    canvas->SaveAs(outputName);

    // Clean up memory
    delete h_Qxb;
    delete rootFile;

}
*/

void Create2DHistogram(const char* inputFile, const char* outputFile) {
    // Open the ROOT file
    TFile file(inputFile, "READ");

    // Read the TH1 histograms using dynamic_cast
    TH1F* Q2pre = dynamic_cast<TH1F*>(file.Get("Q2pre"));
    TH1F* xbpre = dynamic_cast<TH1F*>(file.Get("xbpre"));

    if (!Q2pre || !xbpre) {
        // Check if the dynamic_cast was successful
        // Handle the case where the histograms are not TH1D
        std::cerr << "Error: Could not read TH1 histograms from the file." << std::endl;
        file.Close();
        return;
    }

    // Create a new TH2D using the values from Q2pre and xbpre
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

    // Create a canvas to plot the TH2D
    TCanvas canvas("canvas", "2D Histogram", 800, 600);
    canvas.cd();
    
    // Plot the TH2D with the "colz" representation
    th2d->Draw("colz");

    // Save the plot as a PDF file
    canvas.SaveAs(outputFile);

    // Close the ROOT file
    file.Close();
}



int main() {
    const char* file1 = "/home/matamoros/Desktop/reboot/CleanCode/build/chop_ULT_MCLD2.root";
    const char* file2 = "/home/matamoros/Desktop/reboot/CleanCode/build/chop_ULT_LD2.root";
    CompareHistograms(file1, file2);
    //Create2DHistogram(file1, "Qpre_vs_xbpre_File1.pdf");
    Create2DHistogram(file2, "test2D.pdf");

    return 0;
}
