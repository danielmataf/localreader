#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <iostream>
#include <string>
#include "constants.h" // Include your actual constants header



// g++ thetabins.cpp `root-config --cflags --libs` -o thetabins
// ./thetabins REcutC2_test.root C2 RGD
void draw_xQ2_with_grid(const std::string& filename, const std::string& targetName, const std::string& tag) {
    TFile* f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    //get histo 
    std::string histName = "xQ2_" + targetName + "_" + tag;
    auto* hist = dynamic_cast<TH2F*>(f->Get(histName.c_str()));
    if (!hist) {
        std::cerr << "Error: Histogram " << histName << " not found in file." << std::endl;
        f->Close();
        return;
    }

    TCanvas* c = new TCanvas("c", "xQ2 with Grid", 800, 700);
    gStyle->SetOptStat(0);
    hist->Draw("COLZ");
    hist->GetXaxis()->SetTitle("x_{B}");
    hist->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
    //hist->GetXaxis()->CenterTitle();
    //hist->GetYaxis()->CenterTitle();
    //hist->GetXaxis()->SetTitleSize(0.05);
    //hist->GetYaxis()->SetTitleSize(0.05);
    //hist->GetXaxis()->SetLabelSize(0.04);
    //hist->GetYaxis()->SetLabelSize(0.04);

    double xmin = Constants::Rcutminx;
    double xmax = Constants::Rcutmaxx;
    double Qmin = Constants::RcutminQ;
    double Qmax = Constants::RcutmaxQ;

    int nLines = 6;
    double dx = (xmax - xmin) / nLines;
    double dQ = (Qmax - Qmin) / nLines;

    //vertical
    for (int i = 1; i < nLines; ++i) {
        double x = xmin + i * dx;
        TLine* vline = new TLine(x, Qmin, x, Qmax);
        vline->SetLineColor(kBlack);
        vline->SetLineStyle(2);
        vline->Draw("same");
    }

    //horizontal 
    for (int j = 1; j < nLines; ++j) {
        double q = Qmin + j * dQ;
        TLine* hline = new TLine(xmin, q, xmax, q);
        hline->SetLineColor(kBlack);
        hline->SetLineStyle(2);
        hline->Draw("same");
    }

    std::string outName = "xQ2_withGrid_" + targetName + "_" + tag + ".pdf";
    c->SaveAs(outName.c_str());

    f->Close();
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: ./thetabins <ROOT file> <TargetName (e.g., C2)> <Tag (e.g., RGD)>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    std::string targetName = argv[2];
    std::string tag = argv[3];

    draw_xQ2_with_grid(filename, targetName, tag);
    return 0;
}