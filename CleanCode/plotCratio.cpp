#include <TFile.h>
#include <TH3F.h>
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

const int Rbin = Constants::Rbin_nu; // xB-based binning for Cratio

struct CratioMatrixStruct {
    double val[50][50][50];
    double err[50][50][50];
};
void computeCratio(const std::string& fileName, const std::string& tag, CratioMatrixStruct& matrix, TH3F*& refHisto) {
    std::cout << "[computeCratio] Opening file: " << fileName << std::endl;
    TFile* f = new TFile(fileName.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "[computeCratio] Failed to open file: " << fileName << std::endl;
        return;
    }

    std::string nameD = "wCphi_D_" + tag;
    std::string nameA = "wCphi_A_" + tag;
    std::string nameSqD = "sqCphi_D_" + tag;
    std::string nameSqA = "sqCphi_A_" + tag;
    std::string name3dD = "cphi_3D_D_" + tag;
    std::string name3dA = "cphi_3D_A_" + tag;

    std::cout << "[computeCratio] Looking for histograms:\n  "
              << nameD << "\n  " << nameA << "\n  "
              << nameSqD << "\n  " << nameSqA << "\n  "
              << name3dD << "\n  " << name3dA << std::endl;

    auto* hD = dynamic_cast<TH3F*>(f->Get(nameD.c_str()));
    auto* hA = dynamic_cast<TH3F*>(f->Get(nameA.c_str()));
    auto* hSqD = dynamic_cast<TH3F*>(f->Get(nameSqD.c_str()));
    auto* hSqA = dynamic_cast<TH3F*>(f->Get(nameSqA.c_str()));
    auto* h3dD = dynamic_cast<TH3F*>(f->Get(name3dD.c_str()));
    auto* h3dA = dynamic_cast<TH3F*>(f->Get(name3dA.c_str()));

    if (!hD || !hA || !hSqD || !hSqA || !h3dD || !h3dA) {
        std::cerr << "[computeCratio] One or more histograms missing in file: " << fileName << std::endl;
        if (!hD)    std::cerr << "  --> MISSING: " << nameD << std::endl;
        if (!hA)    std::cerr << "  --> MISSING: " << nameA << std::endl;
        if (!hSqD)  std::cerr << "  --> MISSING: " << nameSqD << std::endl;
        if (!hSqA)  std::cerr << "  --> MISSING: " << nameSqA << std::endl;
        if (!h3dD)  std::cerr << "  --> MISSING: " << name3dD << std::endl;
        if (!h3dA)  std::cerr << "  --> MISSING: " << name3dA << std::endl;
        f->Close();
        return;
    }

    if (!refHisto) {
        refHisto = hD;
        refHisto->SetDirectory(nullptr);
    }

    int NX = hD->GetNbinsX();
    int NY = hD->GetNbinsY();
    int NZ = hD->GetNbinsZ();

    for (int x = 1; x <= NX; ++x) {
        for (int y = 1; y <= NY; ++y) {
            for (int z = 1; z <= NZ; ++z) {
                double valD = hD->GetBinContent(x, y, z);
                double valA = hA->GetBinContent(x, y, z);
                double countD = h3dD->GetBinContent(x, y, z);
                double countA = h3dA->GetBinContent(x, y, z);
                double sqD = hSqD->GetBinContent(x, y, z);
                double sqA = hSqA->GetBinContent(x, y, z);

                double avgD = (countD > 0) ? valD / countD : 0.0;
                double avgA = (countA > 0) ? valA / countA : 0.0;
                double Cratio = (avgD != 0.0) ? avgA / avgD : 0.0;

                double varD = (countD > 0) ? sqD / countD - avgD * avgD : 0.0;
                double varA = (countA > 0) ? sqA / countA - avgA * avgA : 0.0;
                double errD = (countD > 0) ? std::sqrt(varD / countD) : 0.0;
                double errA = (countA > 0) ? std::sqrt(varA / countA) : 0.0;

                double errCratio = (avgD > 0 && avgA > 0) ?
                    Cratio * std::sqrt(std::pow(errA / avgA, 2) + std::pow(errD / avgD, 2)) : 0.0;

                matrix.val[x - 1][y - 1][z - 1] = Cratio;
                matrix.err[x - 1][y - 1][z - 1] = errCratio;
            }
        }
    }

    f->Close();
    std::cout << "[computeCratio] Calculation complete from: " << fileName << std::endl;
}


void drawCratioPdf(const CratioMatrixStruct& A, const CratioMatrixStruct& B, const CratioMatrixStruct& C, TH3F* binref) {
    if (!binref) {
        std::cerr << "[drawCratioPdf] No reference histogram provided. Cannot plot.\n";
        return;
    }

    std::cout << "[drawCratioPdf] Starting plot using reference: " << binref->GetName() << std::endl;
    TCanvas canvas("c", "canvas", 1200, 800);
    std::string outPdf = "Cratiovalues/CratioTripleTarget_fromRoot.pdf";

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
                gC->SetPoint(z, zVal, valC); gC->SetPointError(z, 0, errC);
                gCu->SetPoint(z, zVal + 0.01, valCu); gCu->SetPointError(z, 0, errCu);
                gSn->SetPoint(z, zVal + 0.02, valSn); gSn->SetPointError(z, 0, errSn);
            }

            TMultiGraph* mg = new TMultiGraph();
            gC->SetMarkerStyle(20); gC->SetMarkerColor(kOrange);
            gCu->SetMarkerStyle(20); gCu->SetMarkerColor(kGreen);
            gSn->SetMarkerStyle(20); gSn->SetMarkerColor(kBlack);
            mg->Add(gC); mg->Add(gCu); mg->Add(gSn);

            mg->Draw("APE1");
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("Cratio");
            mg->GetYaxis()->SetRangeUser(0.5, 1.5);

            TLegend* legend = new TLegend(0.15, 0.15, 0.35, 0.30);
            legend->SetTextSize(0.035);
            legend->AddEntry(gC, "C", "lp");
            legend->AddEntry(gCu, "Cu", "lp");
            legend->AddEntry(gSn, "Sn", "lp");
            legend->Draw("same");

            TLine* line = new TLine(0.0, 1.0, 1.0, 1.0);
            line->SetLineStyle(2);
            line->Draw("same");

            TLatex text;
            text.SetTextSize(0.045);
            text.DrawLatexNDC(0.45, 0.22, Form("%.2f < x < %.2f", binref->GetXaxis()->GetBinLowEdge(x + 1), binref->GetXaxis()->GetBinUpEdge(x + 1)));
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
    CratioMatrixStruct cC, cCu, cSn;
    TH3F* refHist = nullptr;

    computeCratio("CratioInputFiles/CratioHist_C2_RGD.root", "C2_RGD", cC, refHist);
    computeCratio("CratioInputFiles/CratioHist_Cu_RGD.root", "Cu_RGD", cCu, refHist);
    computeCratio("CratioInputFiles/CratioHist_Sn_RGD.root", "Sn_RGD", cSn, refHist);

    if (refHist) {
        drawCratioPdf(cC, cCu, cSn, refHist);
        std::cout << "[plotCratio] Cratio PDF generated in ../Cratiovalues/" << std::endl;
    } else {
        std::cerr << "[plotCratio] No valid reference histogram.\n";
    }

    return 0;
}
