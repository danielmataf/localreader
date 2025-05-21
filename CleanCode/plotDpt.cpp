#include <TFile.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TNamed.h>
#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "constants.h"
#include <THnSparse.h>
//How to copile and run 
//g++ plotDpt.cpp -o plotDpt `root-config --cflags --libs`
//./plotDpt   -- (default) uses root input
//./plotDpt 1 -- force txt input (not good, herited from R)

const int Rbin = Constants::Rbin_nu; // xB-based binning for Dpt
const int FIXED_SIZE = Constants::Rbin_nu; //acquired from R 

// == REMEMBER ==
// 1. variables in 3D are xb, Q2, z (no hadron or electron specific restricitions)
// 2. variables in 5D are xb, Q2, nu, z, phih? (no need for electron restriction counts)

struct DptMatrixStruct {
    double val[50][50][50];
    double err[50][50][50];
};

struct DptMatrix {
    std::vector<std::vector<std::vector<double>>> valbis;  //logic remains as result and calced err
    std::vector<std::vector<std::vector<double>>> errbis;  //logic remains as result and calced err
    DptMatrix() :
        valbis(FIXED_SIZE, std::vector<std::vector<double>>(FIXED_SIZE, std::vector<double>(FIXED_SIZE, 0.0))),
        errbis(FIXED_SIZE, std::vector<std::vector<double>>(FIXED_SIZE, std::vector<double>(FIXED_SIZE, 0.0))) {}
};


void computeDpt(const std::string& fileName, const std::string& tag, DptMatrixStruct& matrix, TH3F*& refHisto) {
    std::cout << "[computeDpt] Opening file: " << fileName << std::endl;
    TFile* f = new TFile(fileName.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "[computeDpt] Failed to open file: " << fileName << std::endl;
        return;
    }
    //sanity checks
    std::string name_ptD   = "wpt2_D_" + tag;
    std::string name_ptA   = "wpt2_A_" + tag;
    std::string name_sqD   = "sqpt2_D_" + tag;
    std::string name_sqA   = "sqpt2_A_" + tag;
    std::string name_3dD   = "pt2_3D_D_" + tag;
    std::string name_3dA   = "pt2_3D_A_" + tag;

    auto* h_ptD    = dynamic_cast<TH3F*>(f->Get(name_ptD.c_str()));
    auto* h_ptA    = dynamic_cast<TH3F*>(f->Get(name_ptA.c_str()));
    auto* h_sqptD  = dynamic_cast<TH3F*>(f->Get(name_sqD.c_str()));
    auto* h_sqptA  = dynamic_cast<TH3F*>(f->Get(name_sqA.c_str()));
    auto* h_countD = dynamic_cast<TH3F*>(f->Get(name_3dD.c_str()));
    auto* h_countA = dynamic_cast<TH3F*>(f->Get(name_3dA.c_str()));

    std::cout << "[computeDpt] Looking for histograms:\n  "
              << name_ptD << "\n  " << name_ptA << "\n  "
              << name_sqD << "\n  " << name_sqA << "\n  "
              << name_3dD << "\n  " << name_3dA << std::endl;

    if (!h_ptD || !h_ptA || !h_sqptD || !h_sqptA || !h_countD || !h_countA) {
        std::cerr << "[computeDpt] Missing one or more required histograms in file: " << fileName << std::endl;
        if (!h_ptD)    std::cerr << "  --> " << name_ptD << " missing" << std::endl;
        if (!h_ptA)    std::cerr << "  --> " << name_ptA << " missing" << std::endl;
        if (!h_sqptD)  std::cerr << "  --> " << name_sqD << " missing" << std::endl;
        if (!h_sqptA)  std::cerr << "  --> " << name_sqA << " missing" << std::endl;
        if (!h_countD) std::cerr << "  --> " << name_3dD << " missing" << std::endl;
        if (!h_countA) std::cerr << "  --> " << name_3dA << " missing" << std::endl;
        f->Close();
        return;
    }

    if (!refHisto) {
        refHisto = h_ptD;
        refHisto->SetDirectory(nullptr); //acquired from R 
    }

    int NX = h_ptD->GetNbinsX();
    int NY = h_ptD->GetNbinsY();
    int NZ = h_ptD->GetNbinsZ();    //ok

    for (int x = 1; x <= NX; ++x) {
        for (int y = 1; y <= NY; ++y) {
            for (int z = 1; z <= NZ; ++z) {
                //proceeding straight to third loop since we dont have strict electron counts here as in R 
                double valD = h_ptD->GetBinContent(x, y, z);
                double valA = h_ptA->GetBinContent(x, y, z);
                double countD = h_countD->GetBinContent(x, y, z);
                double countA = h_countA->GetBinContent(x, y, z);
                double sqD = h_sqptD->GetBinContent(x, y, z);
                double sqA = h_sqptA->GetBinContent(x, y, z);
                //lgtm
                double avgD = (countD > 0) ? valD / countD : 0.0;
                double avgA = (countA > 0) ? valA / countA : 0.0;
                double dpt = avgA - avgD;

                double varD = (countD > 0) ? sqD / countD - avgD * avgD : 0.0;
                double varA = (countA > 0) ? sqA / countA - avgA * avgA : 0.0;
                double errD = (countD > 0) ? std::sqrt(varD / countD) : 0.0;
                double errA = (countA > 0) ? std::sqrt(varA / countA) : 0.0;
                double errDpt = std::sqrt(errD * errD + errA * errA);

                if (!std::isfinite(dpt) || !std::isfinite(errDpt)) {
                    std::cerr << "[computeDpt::NaN] x=" << x << " y=" << y << " z=" << z
                              << " | avgA=" << avgA << ", avgD=" << avgD
                              << " | countA=" << countA << ", countD=" << countD
                              << " | dpt=" << dpt << ", err=" << errDpt << std::endl;
                }

                matrix.val[x - 1][y - 1][z - 1] = dpt;
                matrix.err[x - 1][y - 1][z - 1] = errDpt;
            }
        }
    }

    f->Close();
    std::cout << "[computeDpt] Successfully loaded and computed from " << fileName << std::endl;
}


void drawDptPdf(const DptMatrixStruct& A, const DptMatrixStruct& B, const DptMatrixStruct& C, TH3F* binref) {
    if (!binref) {
        std::cerr << "[drawDptPdf] No reference histogram provided. Cannot plot.\n";
        return;
    }

    std::cout << "[drawDptPdf] Plotting using reference: " << binref->GetName() << std::endl;
    TCanvas canvas("c", "canvas", 1200, 800);
    std::string outPdf = "Dptvalues/DptTripleTarget_fromRoot.pdf";

    int NX = binref->GetNbinsX();
    int NY = binref->GetNbinsY();
    int NZ = binref->GetNbinsZ();

    for (int x = 0; x < NX; ++x) {
        canvas.Clear();
        canvas.Divide(3, 2);

        for (int y = 0; y < NY; ++y) {
            canvas.cd(y + 1);
            auto* gC = new TGraphErrors();
            auto* gCu = new TGraphErrors();
            auto* gSn = new TGraphErrors();

            for (int z = 0; z < NZ; ++z) {
                double zVal = binref->GetZaxis()->GetBinCenter(z + 1);
                double valC = A.val[x][y][z];
                double valCu = B.val[x][y][z];
                double valSn = C.val[x][y][z];
                double errC = A.err[x][y][z];
                double errCu = B.err[x][y][z];
                double errSn = C.err[x][y][z];
                gSn->SetPoint(z, zVal , valSn); gSn->SetPointError(z, 0, errSn);
                gCu->SetPoint(z, zVal + 0.01, valCu); gCu->SetPointError(z, 0, errCu);
                gC->SetPoint(z, zVal + 0.02, valC); gC->SetPointError(z, 0, errC);
            }

            TMultiGraph* mg = new TMultiGraph();
            gC->SetMarkerStyle(20); gC->SetMarkerColor(kBlack);
            gCu->SetMarkerStyle(20); gCu->SetMarkerColor(kGreen);
            gSn->SetMarkerStyle(20); gSn->SetMarkerColor(kOrange);
            mg->Add(gC); mg->Add(gCu); mg->Add(gSn);

            mg->Draw("APE1");
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("#Delta<p_{T}^{2}>");
            mg->GetYaxis()->SetRangeUser(-0.05, 0.15  );

            TLegend* legend = new TLegend(0.15, 0.15, 0.35, 0.30);
            legend->SetTextSize(0.035);
            legend->AddEntry(gC, "C", "lp");
            legend->AddEntry(gCu, "Cu", "lp");
            legend->AddEntry(gSn, "Sn", "lp");
            legend->Draw("same");

            TLine* line = new TLine(0.0, 0.0, 1.0, 0.0);
            line->SetLineStyle(2);
            line->Draw("same");

            TLatex text;
            text.SetTextSize(0.045);
            text.DrawLatexNDC(0.45, 0.22, Form("%.2f < x_{B} < %.2f", binref->GetXaxis()->GetBinLowEdge(x + 1), binref->GetXaxis()->GetBinUpEdge(x + 1)));
            text.DrawLatexNDC(0.45, 0.17, Form("%.2f < Q^{2} < %.2f", binref->GetYaxis()->GetBinLowEdge(y + 1), binref->GetYaxis()->GetBinUpEdge(y + 1)));

            TLatex watermark;
            watermark.SetTextSize(0.08);
            watermark.SetTextAngle(45);
            watermark.SetTextColorAlpha(kGray + 1, 0.3);
            watermark.SetNDC();
            watermark.SetTextAlign(22);
            watermark.DrawLatex(0.5, 0.5, "very preliminary");
        }
        if (x == 0)
            canvas.Print((outPdf + "(").c_str());
        else if (x == NX - 1)
            canvas.Print((outPdf + ")").c_str());
        else
            canvas.Print(outPdf.c_str());
    }
}

int main() {
    DptMatrixStruct dptC, dptCu, dptSn;
    TH3F* refHist = nullptr;

    computeDpt("DptinputFiles/Dhist_C2_RGD.root", "C2_RGD", dptC, refHist);
    computeDpt("DptinputFiles/Dhist_Cu_RGD.root", "Cu_RGD", dptCu, refHist);
    computeDpt("DptinputFiles/Dhist_Sn_RGD.root", "Sn_RGD", dptSn, refHist);

    if (refHist) {
        drawDptPdf(dptC, dptCu, dptSn, refHist);
        std::cout << "[plotDpt] Dpt PDF generated in ../Dptvalues/" << std::endl;
    } else {
        std::cerr << "[plotDpt] No valid reference histogram.\n";
    }
    return 0;
}