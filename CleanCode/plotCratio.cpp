#include <TFile.h>
#include <TH3F.h>
#include <TH2F.h>
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
const int FIXED_SIZE = Constants::Rbin_nu;

/*
struct CratioMatrixStruct {
    double val[50][50][50];
    double err[50][50][50];
};

struct RatioMatrix {
    std::vector<std::vector<std::vector<double>>> valbis;
    std::vector<std::vector<std::vector<double>>> errbis;
    CratioMatrix() :
        valbis(FIXED_SIZE, std::vector<std::vector<double>>(FIXED_SIZE, std::vector<double>(FIXED_SIZE, 0.0))),
        errbis(FIXED_SIZE, std::vector<std::vector<double>>(FIXED_SIZE, std::vector<double>(FIXED_SIZE, 0.0))) {}
};

void computeCratio(const std::string& fileName, const std::string& tag, CratioMatrixStruct& matrix, TH3F*& refHisto) {
    std::cout << "[computeCratio] Opening file: " << fileName << std::endl;
    TFile* f = new TFile(fileName.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "[computeCratio] Failed to open file: " << fileName << std::endl;
        return;
    }
    //starting sanity checks
    std::string nameD = "wCphi_D_" + tag;   //weighted   
    std::string nameA = "wCphi_A_" + tag;   //weighted
    std::string nameSqD = "sqCphi_D_" + tag;    //weighted square
    std::string nameSqA = "sqCphi_A_" + tag;    //weighted square
    std::string name3dD = "cphi_3D_D_" + tag;   //count
    std::string name3dA = "cphi_3D_A_" + tag;   //count

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
        refHisto->SetDirectory(nullptr); //acquired from R
    }
    // ok as usual
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

                double avgD = (countD > 0) ? valD / countD : 0.0;   //weighted average
                double avgA = (countA > 0) ? valA / countA : 0.0;   //weighted average
                double Cratio = (avgD != 0.0) ? avgA / avgD : 0.0;  //ratio of two weighted avgs

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
*/
// REGIONS compute and drawing //
struct CratioMatrixSimple {
    double val[4][4] = {{0.0}};
    double err[4][4] = {{0.0}};
};

void computeCratioRegion4x4(const std::string& fileName, const std::string& tag, char region, CratioMatrixSimple& matrix) {
    std::cout << "[computeCratioRegion4x4] Opening file: " << fileName << std::endl;
    TFile* f = new TFile(fileName.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "[computeCratioRegion4x4] Failed to open file." << std::endl;
        return;
    }

    std::ostringstream nameA, nameD, nameSqA, nameSqD, nameCountA, nameCountD;
    nameA       << "region" << region << "_A_" << tag;
    nameD       << "region" << region << "_D_" << tag;
    nameSqA     << "region" << region << "_w2A_" << tag;
    nameSqD     << "region" << region << "_w2D_" << tag;
    nameCountA  << "region" << region << "_count_A_" << tag;
    nameCountD  << "region" << region << "_count_D_" << tag;

    auto* hA = dynamic_cast<TH2F*>(f->Get(nameA.str().c_str()));
    auto* hD = dynamic_cast<TH2F*>(f->Get(nameD.str().c_str()));
    auto* hSqA = dynamic_cast<TH2F*>(f->Get(nameSqA.str().c_str()));
    auto* hSqD = dynamic_cast<TH2F*>(f->Get(nameSqD.str().c_str()));
    auto* hCountA = dynamic_cast<TH2F*>(f->Get(nameCountA.str().c_str()));
    auto* hCountD = dynamic_cast<TH2F*>(f->Get(nameCountD.str().c_str()));

    if (!hA || !hD || !hSqA || !hSqD || !hCountA || !hCountD) {
    std::cerr << "[computeCratioRegion4x4] Missing one or more histograms.\n";
        if (!hA)     std::cerr << " - Missing: " << nameA.str()     << "\n";
        if (!hD)     std::cerr << " - Missing: " << nameD.str()     << "\n";
        if (!hSqA)    std::cerr << " - Missing: " << nameSqA.str()    << "\n";
        if (!hSqD)    std::cerr << " - Missing: " << nameSqD.str()    << "\n";
        if (!hCountA) std::cerr << " - Missing: " << nameCountA.str() << "\n";
        if (!hCountD) std::cerr << " - Missing: " << nameCountD.str() << "\n";
        return;
    }
    for (int ipt2 = 1; ipt2 <= 4; ++ipt2) {
        for (int iz = 1; iz <= 4; ++iz) {
            double wA = hA->GetBinContent(ipt2, iz);
            double wD = hD->GetBinContent(ipt2, iz);
            double sqA = hSqA->GetBinContent(ipt2, iz);
            double sqD = hSqD->GetBinContent(ipt2, iz);
            double nA = hCountA->GetBinContent(ipt2, iz);
            double nD = hCountD->GetBinContent(ipt2, iz);

            if (nA > 0 && nD > 0) {
                double meanA = wA / nA;
                double meanD = wD / nD;
                double errA = std::sqrt(std::abs(sqA / nA - meanA * meanA) / nA);
                double errD = std::sqrt(std::abs(sqD / nD - meanD * meanD) / nD);

                matrix.val[ipt2 - 1][iz - 1] = meanA / meanD;
                matrix.err[ipt2 - 1][iz - 1] = matrix.val[ipt2 - 1][iz - 1] * std::sqrt((errA / meanA) * (errA / meanA) + (errD / meanD) * (errD / meanD));
            } else {
                matrix.val[ipt2 - 1][iz - 1] = 0;
                matrix.err[ipt2 - 1][iz - 1] = 0;
            }
        }
    }
    f->Close();
}
void drawCratioRegionPdf(const CratioMatrixSimple& C2, const CratioMatrixSimple& Cu, const CratioMatrixSimple& Sn, const std::string& regionTag) {
    const std::string refTag = "C2";  // ← pick any reliable target
    const std::string refFilePath = "CratioInputFiles/CratioHist_" + refTag + "_RGD.root";

    TFile* f = TFile::Open(refFilePath.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "[drawCratioRegionPdf] Failed to open reference ROOT file: " << refFilePath << std::endl;
        return;
    }

    std::string refHistName = "region" + regionTag + "_D_" + refTag + "_RGD";
    TH2F* binref = (TH2F*) f->Get(refHistName.c_str());

    if (!binref) {
        std::cerr << "[drawCratioRegionPdf] Reference histogram not found: " << refHistName << std::endl;
        f->Close();
        return;
    }

    std::string outPdf = "Cratiovalues/CphiRegion_" + regionTag + ".pdf";
    TCanvas canvas("c", "CphiRatio", 1200, 800);
    canvas.Divide(2, 2);  // 4 pt2 bins per page

    const int Npt2 = 4;
    const int Nz = 4;

    for (int ipt2 = 0; ipt2 < Npt2; ++ipt2) {
        canvas.cd(ipt2 + 1);

        auto* gC2 = new TGraphErrors();
        auto* gCu = new TGraphErrors();
        auto* gSn = new TGraphErrors();

        for (int iz = 0; iz < Nz; ++iz) {
            const double zVal = binref->GetYaxis()->GetBinCenter(iz + 1);

            double valC2 = C2.val[ipt2][iz];
            double valCu = Cu.val[ipt2][iz];
            double valSn = Sn.val[ipt2][iz];
            double errC2 = C2.err[ipt2][iz];
            double errCu = Cu.err[ipt2][iz];
            double errSn = Sn.err[ipt2][iz];

            gSn->SetPoint(iz, zVal, valSn);        gSn->SetPointError(iz, 0.0, errSn);
            gCu->SetPoint(iz, zVal + 0.01, valCu); gCu->SetPointError(iz, 0.0, errCu);
            gC2->SetPoint(iz, zVal + 0.02, valC2); gC2->SetPointError(iz, 0.0, errC2);
        }

        gC2->SetMarkerStyle(20); gC2->SetMarkerColor(kBlack);
        gCu->SetMarkerStyle(20); gCu->SetMarkerColor(kGreen);
        gSn->SetMarkerStyle(20); gSn->SetMarkerColor(kOrange);

        TMultiGraph* mg = new TMultiGraph();
        mg->Add(gC2); mg->Add(gCu); mg->Add(gSn);
        mg->Draw("APE1");
        mg->GetXaxis()->SetTitle("z");
        mg->GetYaxis()->SetTitle("#LTcos(#phi_{h})#GT_{A}/#LTcos(#phi_{h})#GT_{D}");
        mg->GetYaxis()->SetRangeUser(0.5, 1.5);

        TLegend* legend = new TLegend(0.15, 0.15, 0.35, 0.30);
        legend->SetTextSize(0.035);
        legend->AddEntry(gC2, "C", "lp");
        legend->AddEntry(gCu, "Cu", "lp");
        legend->AddEntry(gSn, "Sn", "lp");
        legend->Draw("same");

        TLine* line = new TLine(0.0, 1.0, 1.0, 1.0);
        line->SetLineStyle(2);
        line->Draw("same");

        TLatex text;
        text.SetTextSize(0.045);
        text.DrawLatexNDC(0.45, 0.22, Form("%.2f < p_{T}^{2} < %.2f",
            binref->GetXaxis()->GetBinLowEdge(ipt2 + 1),
            binref->GetXaxis()->GetBinUpEdge(ipt2 + 1)));

        TLatex watermark;
        watermark.SetTextSize(0.08);
        watermark.SetTextAngle(45);
        watermark.SetTextColorAlpha(kGray + 1, 0.3);
        watermark.SetNDC();
        watermark.SetTextAlign(22);
        watermark.DrawLatex(0.5, 0.5, "very preliminary");
    }

    canvas.Print(outPdf.c_str());
    f->Close();
}


// end REGIONS //


int main() {
    std::cout << "[main] testing.\n";
//    CratioMatrixSimple mat_A;
//    TH2F* refHist = nullptr;
//
//    computeCratioRegion4x4("CratioInputFiles/CratioHist_Cu_RGD.root", "Cu_RGD", 'A', mat_A, refHist);
//    drawCratioRegionPdf(mat_A, refHist, "A");
CratioMatrixSimple mat_C2_A, mat_Cu_A, mat_Sn_A;
CratioMatrixSimple mat_C2_B, mat_Cu_B, mat_Sn_B;
CratioMatrixSimple mat_C2_C, mat_Cu_C, mat_Sn_C;
CratioMatrixSimple mat_C2_D, mat_Cu_D, mat_Sn_D;

// Region A
computeCratioRegion4x4("CratioInputFiles/CratioHist_C2_RGD.root", "C2_RGD", 'A', mat_C2_A);
computeCratioRegion4x4("CratioInputFiles/CratioHist_Cu_RGD.root", "Cu_RGD", 'A', mat_Cu_A);
computeCratioRegion4x4("CratioInputFiles/CratioHist_Sn_RGD.root", "Sn_RGD", 'A', mat_Sn_A);
drawCratioRegionPdf(mat_C2_A, mat_Cu_A, mat_Sn_A, "A");

// Region B
computeCratioRegion4x4("CratioInputFiles/CratioHist_C2_RGD.root", "C2_RGD", 'B', mat_C2_B);
computeCratioRegion4x4("CratioInputFiles/CratioHist_Cu_RGD.root", "Cu_RGD", 'B', mat_Cu_B);
computeCratioRegion4x4("CratioInputFiles/CratioHist_Sn_RGD.root", "Sn_RGD", 'B', mat_Sn_B);
drawCratioRegionPdf(mat_C2_B, mat_Cu_B, mat_Sn_B, "B");

// Region C
computeCratioRegion4x4("CratioInputFiles/CratioHist_C2_RGD.root", "C2_RGD", 'C', mat_C2_C);
computeCratioRegion4x4("CratioInputFiles/CratioHist_Cu_RGD.root", "Cu_RGD", 'C', mat_Cu_C);
computeCratioRegion4x4("CratioInputFiles/CratioHist_Sn_RGD.root", "Sn_RGD", 'C', mat_Sn_C);
drawCratioRegionPdf(mat_C2_C, mat_Cu_C, mat_Sn_C, "C");

// Region D
computeCratioRegion4x4("CratioInputFiles/CratioHist_C2_RGD.root", "C2_RGD", 'D', mat_C2_D);
computeCratioRegion4x4("CratioInputFiles/CratioHist_Cu_RGD.root", "Cu_RGD", 'D', mat_Cu_D);
computeCratioRegion4x4("CratioInputFiles/CratioHist_Sn_RGD.root", "Sn_RGD", 'D', mat_Sn_D);
drawCratioRegionPdf(mat_C2_D, mat_Cu_D, mat_Sn_D, "D");


    return 0;
}
