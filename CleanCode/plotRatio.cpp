#include <TFile.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "constants.h"
#include <THnSparse.h>


//compile and run: 
//g++ plotRatio.cpp -o plotRatio `root-config --cflags --libs`
//./plotRatio   -- (default) uses txt input --also w/ 1
//./plotRatio 2 -- uses root input
const int Rbin = Constants::Rbin_nu; //exporting value from csts
const int FIXED_SIZE = Constants::Rbin_nu;

//struct RatioMatrix {
//    double val[50][50][50];
//    double err[50][50][50];
//    //ratMatrix(Rbin, std::vector<std::vector<double>>(Rbin, std::vector<double>(Rbin, 0.0))),
//    //errorMatrix(Rbin, std::vector<std::vector<double>>(Rbin, std::vector<double>(Rbin, 0.0))),
//};
struct RatioMatrix {
    std::vector<std::vector<std::vector<double>>> val;
    std::vector<std::vector<std::vector<double>>> err;
    RatioMatrix() :
        val(FIXED_SIZE, std::vector<std::vector<double>>(FIXED_SIZE, std::vector<double>(FIXED_SIZE, 0.0))),
        err(FIXED_SIZE, std::vector<std::vector<double>>(FIXED_SIZE, std::vector<double>(FIXED_SIZE, 0.0))) {}
};

void computeRatio(const std::string& fileName, const std::string& tag, RatioMatrix& outMatrix, TH3F*& refHisto) {
    std::cout << "[computeRatio] Attempting to open file: " << fileName << std::endl;
    TFile* f = new TFile(fileName.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "[computeRatio] Failed to open file: " << fileName << std::endl;
        return;
    }

    std::string nameD   = "nu_z_pt2_D_" + tag;
    std::string nameA   = "nu_z_pt2_A_" + tag;
    std::string nameNuD = "nu_D_" + tag;
    std::string nameNuA = "nu_A_" + tag;

    std::cout << "\n[computeRatio] Looking for histograms:\n  " << nameD << "\n  " << nameA << "\n  " << nameNuD << "\n  " << nameNuA << std::endl;

    auto* hD   = dynamic_cast<TH3F*>(f->Get(nameD.c_str()));
    auto* hA   = dynamic_cast<TH3F*>(f->Get(nameA.c_str()));
    auto* hnuD = dynamic_cast<TH1F*>(f->Get(nameNuD.c_str()));
    auto* hnuA = dynamic_cast<TH1F*>(f->Get(nameNuA.c_str()));

    if (!hD) std::cerr << "[computeRatio] Histogram not found: " << nameD << std::endl;
    if (!hA) std::cerr << "[computeRatio] Histogram not found: " << nameA << std::endl;
    if (!hnuD) std::cerr << "[computeRatio] Histogram not found: " << nameNuD << std::endl;
    if (!hnuA) std::cerr << "[computeRatio] Histogram not found: " << nameNuA << std::endl;

    if (!hD || !hA || !hnuD || !hnuA) {
        std::cerr << "[computeRatio] One or more histograms missing. Skipping file.\n" << std::endl;
        f->Close();
        return;
    }

    if (!refHisto) {
        refHisto = hD;
        refHisto->SetDirectory(nullptr);  // Prevent deletion on file close
    }

    int NX = hD->GetNbinsX();
    int NY = hD->GetNbinsY();
    int NZ = hD->GetNbinsZ();

    for (int x = 1; x <= NX; ++x) {
        double nuD = hnuD->GetBinContent(x);
        double nuA = hnuA->GetBinContent(x);
        for (int y = 1; y <= NY; ++y) {
            for (int z = 1; z <= NZ; ++z) {
                double valD = hD->GetBinContent(x, y, z);
                double valA = hA->GetBinContent(x, y, z);

                double rD = (nuD > 0) ? valD / nuD : 0.0;
                double rA = (nuA > 0) ? valA / nuA : 0.0;
                double R = (rD > 0) ? rA / rD : 0.0;
                double err = (valD > 0 && valA > 0 && nuD > 0 && nuA > 0) ?
                    R * std::sqrt(1.0 / valA + 1.0 / valD + 1.0 / nuA + 1.0 / nuD) : 0.0;

                if (!std::isfinite(R) || !std::isfinite(err)) {
                    std::cout << "[NaN Check] R= " << R << ", err= " << err
                              << " | valA= " << valA << ", valD= " << valD
                              << ", nuA= " << nuA << ", nuD= " << nuD
                              << " | bin(x,y,z)= (" << x << ", " << y << ", " << z << ")" << std::endl;
                }

                outMatrix.val[x - 1][y - 1][z - 1] = R;
                outMatrix.err[x - 1][y - 1][z - 1] = err;
            }
        }
    }
    f->Close();
}

const int Rdim = 5;  // assuming your THnSparseD has 5 dimensions

void checkTHnSparseContents(const std::string& filePath, const std::string& histoName) {
    TFile* file = TFile::Open(filePath.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[checkTHnSparseContents] Failed to open file: " << filePath << std::endl;
        return;
    }

    auto* h = dynamic_cast<THnSparseD*>(file->Get(histoName.c_str()));
    if (!h) {
        std::cerr << "[checkTHnSparseContents] Histogram " << histoName << " not found in " << filePath << std::endl;
        file->Close();
        return;
    }

    std::cout << "[checkTHnSparseContents] Contents of histogram: " << histoName << std::endl;
    Long64_t nbins = h->GetNbins();
    int ndim = h->GetNdimensions();
    std::vector<Int_t> indices(ndim);

    for (Long64_t i = 0; i < nbins; ++i) {
        double content = h->GetBinContent(i);
        if (content == 0) continue;

        Long64_t local = i;
        for (int d = ndim - 1; d >= 0; --d) {
            Int_t nb = h->GetAxis(d)->GetNbins() + 2; // include under/overflow
            indices[d] = local % nb;
            local /= nb;
        }

        std::cout << "Bin " << i << " | Content = " << content << " | Indices = [ ";
        for (int d = 0; d < ndim; ++d)
            std::cout << indices[d] << " ";
        std::cout << "]" << std::endl;
    }

    file->Close();
}


void drawRatioPdf(const RatioMatrix& A, const RatioMatrix& B, const RatioMatrix& C, TH3F* binref) {
    if (!binref) {
        std::cerr << "[drawRatioPdf] No reference histogram provided. Cannot plot.\n";
        return;
    }

    std::cout << "[drawRatioPdf] Starting plot with reference histogram: " << binref->GetName() << std::endl;

    TCanvas canvas("c", "canvas", 1200, 800);
    std::string outPdf = "Rvalues/RatioTripleTarget_fromRoot.pdf";

    int NX = binref->GetNbinsX();
    int NY = binref->GetNbinsY();
    int NZ = binref->GetNbinsZ();

    for (int x = 0; x < NX; ++x) {
        canvas.Clear();
        canvas.Divide(3, 2);
        double nu = binref->GetXaxis()->GetBinCenter(x + 1);

        for (int z = 0; z < NZ; ++z) {
            canvas.cd(z + 1);
            auto* gC = new TGraphErrors();
            auto* gCu = new TGraphErrors();
            auto* gSn = new TGraphErrors();

            for (int y = 0; y < NY; ++y) {
                double zVal = binref->GetYaxis()->GetBinCenter(y + 1);

                double valC = A.val[x][y][z];
                double valCu = B.val[x][y][z];
                double valSn = C.val[x][y][z];

                double errC = A.err[x][y][z];
                double errCu = B.err[x][y][z];
                double errSn = C.err[x][y][z];

                std::cout << "[drawRatioPdf] Bin (" << x << ", " << y << ", " << z << ")"
                          << " | C= " << valC << "±" << errC
                          << ", Cu= " << valCu << "±" << errCu
                          << ", Sn= " << valSn << "±" << errSn << std::endl;

                gSn->SetPoint(y, zVal , valSn);
                gSn->SetPointError(y, 0.0, errSn);
                gCu->SetPoint(y, zVal + 0.01, valCu);
                gCu->SetPointError(y, 0.0, errCu);
                gC->SetPoint(y, zVal+ 0.02, valC);
                gC->SetPointError(y, 0.0, errC);
            }

            TMultiGraph* mg = new TMultiGraph();
            gC->SetMarkerStyle(20); gC->SetMarkerColor(kBlack);
            gCu->SetMarkerStyle(20); gCu->SetMarkerColor(kGreen);
            gSn->SetMarkerStyle(20); gSn->SetMarkerColor(kOrange);
            mg->Add(gC); mg->Add(gCu); mg->Add(gSn);

            mg->Draw("APE1");
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("R");
            mg->GetYaxis()->SetRangeUser(0, 1.5);

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
            //instead of of bin center we want ranges for variables, so lets do get bin edges 
            text.DrawLatexNDC(0.45, 0.22, Form( "%.2f < #nu < %.2f", binref->GetXaxis()->GetBinLowEdge(x + 1), binref->GetXaxis()->GetBinUpEdge(x + 1)));
            text.DrawLatexNDC(0.45, 0.17, Form("%.2f < p_{T}^{2} < %.2f", binref->GetZaxis()->GetBinLowEdge(z + 1), binref->GetZaxis()->GetBinUpEdge(z + 1)));

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
    RatioMatrix rC, rCu, rSn;
    TH3F *refHist = nullptr;

    computeRatio("RinputFiles/Rhist_C2_RGD.root", "C2_RGD", rC, refHist);
    computeRatio("RinputFiles/Rhist_Cu_RGD.root", "Cu_RGD", rCu, refHist);
    computeRatio("RinputFiles/Rhist_Sn_RGD.root", "Sn_RGD", rSn, refHist);
    if (refHist) {
        drawRatioPdf(rC, rCu, rSn, refHist);
        std::cout << "[plotRatio] R PDF generated in ../Rvalues/" << std::endl;
    } else {
        std::cerr << "[plotRatio] Aborting plot: no valid reference histogram loaded.\n";
    }
    checkTHnSparseContents("RinputFiles/Rhist_C2_RGD.root", "Q2_xB_nu_pt2_z_D_C2_RGD");
    checkTHnSparseContents("RinputFiles/Rhist_Cu_RGD.root", "Q2_xB_nu_pt2_z_D_Cu_RGD");
    checkTHnSparseContents("RinputFiles/Rhist_Sn_RGD.root", "Q2_xB_nu_pt2_z_D_Sn_RGD");

    return 0;
}
