#include <TFile.h>
#include <TH3D.h>
#include <TH2D.h>
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

// -----------------------------------------------------------------------------
//  D⟨pT2⟩ in 4D  :  Q2 -> xB -> phih ->  z
//  counts  :  Q_x_phih_z_D_(tag)      , Q_x_phih_z_A_(tag) // this is for average calc 
//  weight (pt2) = val   :  Q_x_phih_z_wD_(tag)     , Q_x_phih_z_wA_(tag)   //this is actal value collected
//  weighted sq (pt2 *pt2)   :  Q_x_phih_z_w2D_(tag)    , Q_x_phih_z_w2A_(tag) /this is for error calculation ( variance )
// -----------------------------------------------------------------------------
void computeDpt4D_LoopsPhih(const std::string& fileName,
                            const std::string& tag)
{
    std::cout << "[computeDpt4D] Opening file " << fileName << std::endl;
    TFile* f = TFile::Open(fileName.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[computeDpt4D] Cannot open file – abort\n";
        return;
    }
    //hitos attribution 
    const std::string hCntD_name  = "Q_x_phih_z_D_"  + tag;
    const std::string hCntA_name  = "Q_x_phih_z_A_"  + tag;
    const std::string hWPtD_name  = "Q_x_phih_z_wD_" + tag;   // w
    const std::string hWPtA_name  = "Q_x_phih_z_wA_" + tag;
    const std::string hW2PtD_name = "Q_x_phih_z_w2D_" + tag;  // w sq
    const std::string hW2PtA_name = "Q_x_phih_z_w2A_" + tag;

    auto* hCntD  = dynamic_cast<THnSparseD*>(f->Get(hCntD_name .c_str()));
    auto* hCntA  = dynamic_cast<THnSparseD*>(f->Get(hCntA_name .c_str()));
    auto* hWPtD  = dynamic_cast<THnSparseD*>(f->Get(hWPtD_name .c_str()));
    auto* hWPtA  = dynamic_cast<THnSparseD*>(f->Get(hWPtA_name .c_str()));
    auto* hW2PtD = dynamic_cast<THnSparseD*>(f->Get(hW2PtD_name.c_str()));
    auto* hW2PtA = dynamic_cast<THnSparseD*>(f->Get(hW2PtA_name.c_str()));
    //sanity check
    if (!hCntD || !hCntA || !hWPtD || !hWPtA || !hW2PtD || !hW2PtA) {
        std::cerr << "[computeDpt4D] One or more 4-D histograms missing\n";
        f->ls();
        f->Close();
        return;
    }

    // create empty output container (same binning as counts-D all same binnign) 
    // this for results values output 
    auto* hDpt = static_cast<THnSparseD*>(hCntD->Clone(("dpt4D_" + tag).c_str()));
    hDpt->Reset();

    // axis extents
    const int nQ   = hCntD->GetAxis(0)->GetNbins();
    const int nxB  = hCntD->GetAxis(1)->GetNbins();
    const int nPhi = hCntD->GetAxis(2)->GetNbins();
    const int nZ   = hCntD->GetAxis(3)->GetNbins();

    Int_t idx[4];
    Long64_t zeroCounter = 0, totCounter = 0;

    // main loop  (Q² → xB → φh → z)
    for (int iQ = 1;  iQ <= nQ;   ++iQ){
        for (int ix = 1;  ix <= nxB;  ++ix){
            for (int ip = 1;  ip <= nPhi; ++ip){
                for (int iz = 1;  iz <= nZ;   ++iz) {
                    idx[0] = iQ-1;  idx[1] = ix-1;  idx[2] = ip-1;  idx[3] = iz-1;
                    ++totCounter;

                    // -- D target --------------------------------------------------
                    const double cD   = hCntD ->GetBinContent(idx);
                    const double sD   = hWPtD ->GetBinContent(idx);
                    //std::cout << "[computeDpt4D] idx: " << idx[0] << ", " << idx[1] << ", " << idx[2] << ", " << idx[3] << std::endl;
                    //std::cout << "[computeDpt4D] cD: " << cD << ", sD: " << sD << std::endl;
                    const double s2D  = hW2PtD->GetBinContent(idx);
                    const double avgD = (cD>0) ? sD/cD : 0.0;
                    const double varD = (cD>0) ? s2D/cD - avgD*avgD : 0.0;
                    const double errD = (cD>0) ? std::sqrt(varD/cD) : 0.0;

                    // -- A target --------------------------------------------------
                    const double cA   = hCntA ->GetBinContent(idx);
                    const double sA   = hWPtA ->GetBinContent(idx);
                    //std::cout << "[computeDpt4D] cA: " << cA << ", sA: " << sA << std::endl;
                    const double s2A  = hW2PtA->GetBinContent(idx);
                    const double avgA = (cA>0) ? sA/cA : 0.0;
                    const double varA = (cA>0) ? s2A/cA - avgA*avgA : 0.0;
                    const double errA = (cA>0) ? std::sqrt(varA/cA) : 0.0;

                    // -- delta⟨pt2⟩ and uncertainty --------------------------------
                    const double dpt   = avgA - avgD;
                    const double err   = std::sqrt(errA*errA + errD*errD);

                    if (dpt==0.0) ++zeroCounter;
                    //std::cout << "Bin ["
                    //          << idx[0] << " " << idx[1] << " " << idx[2] << " " << idx[3] << "] | "
                    //          << "sD = " << sD << " | cD = " << cD << " || "
                    //          << "sA = " << sA << " | cA = " << cA << " -> Dpt2 = " << dpt
                    //          << std::endl;
//
                    hDpt->SetBinContent(idx, dpt);
                    hDpt->SetBinError  (idx, err);
                }
            }
        }
    }

    std::cout << "[computeDpt4D] Zero-valued bins: "  << zeroCounter << " / " << totCounter << '\n';

    //output 
    TFile* fout = new TFile(("dpt4D_"+tag+".root").c_str(),"RECREATE");
    hDpt->Write();
    fout->Close();
    f->Close();

    std::cout << "[computeDpt4D] Saved D⟨pt2⟩ histogram to dpt4D_" << tag << ".root\n";
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

///basic axis info (bins, range, titles)
void inspectTHnSparseAxes(const THnSparseD* h, const std::string& label)
{
    if (!h) { std::cerr << "[inspectTHnSparse] null pointer for " << label << '\n'; return; }

    std::cout << "\n[inspectTHnSparse]  " << label << '\n';
    const Int_t ndim = h->GetNdimensions();
    for (Int_t ax = 0; ax < ndim; ++ax) {
        auto* axis = h->GetAxis(ax);
        std::cout << "  axis " << ax          << "  (\"" << axis->GetTitle() << "\") : "
                  << axis->GetNbins()         << " bins  |  "
                  << "[" << axis->GetXmin()
                  << ", " << axis->GetXmax()   << "]\n";
    }
    std::cout << "  total sparse bins   : " << h->GetNbins()   << '\n'
          << "  filled entries      : " << h->GetEntries() << '\n';

    Double_t sumW = 0.0;
    Long64_t totalBins = h->GetNbins();
    for (Long64_t i = 0; i < totalBins; ++i) {
        double val = h->GetBinContent(i);
        if (std::isfinite(val)) sumW += val;
    }
    std::cout << "  sum of weights      : " << sumW << "\n" << std::endl;

}

void checkTHnSparseContents(const std::string& filePath,
                            const std::string& histoName,
                            bool printPerBin = true)
{
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
    std::vector<Int_t> coords(ndim);

    for (Long64_t g = 0; g < nbins; ++g) {
        double content = h->GetBinContent(g);
        if (content == 0.0) continue;

        // Decode global bin index to per-dimension indices manually
        Long64_t local = g;
        for (int d = ndim - 1; d >= 0; --d) {
            Int_t nb = h->GetAxis(d)->GetNbins() + 2;
            coords[d] = local % nb;
            local /= nb;
        }

        std::cout << "Bin " << g
                  << " | Content = " << content
                  << " | Indices = [ ";
        for (int i = 0; i < ndim; ++i)
            std::cout << coords[i] << " ";
        std::cout << "]" << std::endl;
    }

    file->Close();
}

// REGIONS //
//annoying but temp, modify when we want to change to 6 points!!
struct DptMatrixSimple {
    double val[4][4];  //phih, z
    double err[4][4];
};
void computeDptRegion4x4(const std::string& fileName, const std::string& tag, char region, DptMatrixSimple& matrix) {
    std::cout << "[computeDptRegion4x4] Opening file: " << fileName << std::endl;
    TFile* f = TFile::Open(fileName.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "   Cannot open file – skipped.\n";
        return;
    }
    
// Build histogram names
std::ostringstream name_ptD, name_ptA, name_sqD, name_sqA, name_countD, name_countA;
name_ptD    << "region" << region << "_D_" << tag <<"_RGD";
name_ptA    << "region" << region << "_A_" << tag <<"_RGD";
name_sqD    << "region" << region << "_w2D_" << tag<<"_RGD";
name_sqA    << "region" << region << "_w2A_" << tag<<"_RGD";
name_countD << "region" << region << "_count_D_" << tag<<"_RGD";
name_countA << "region" << region << "_count_A_" << tag<<"_RGD";

std::cout << "[computeDptRegion4x4] Looking for histograms:\n";
std::cout << "   " << name_ptD.str()    << '\n';
std::cout << "   " << name_ptA.str()    << '\n';
std::cout << "   " << name_sqD.str()    << '\n';
std::cout << "   " << name_sqA.str()    << '\n';
std::cout << "   " << name_countD.str() << '\n';
std::cout << "   " << name_countA.str() << '\n';

// Try to retrieve histograms (use TH2 instead of TH2F to be safe)
TH2* h_ptD    = (TH2*) f->Get(name_ptD.str().c_str());
TH2* h_ptA    = (TH2*) f->Get(name_ptA.str().c_str());
TH2* h_sqptD  = (TH2*) f->Get(name_sqD.str().c_str());
TH2* h_sqptA  = (TH2*) f->Get(name_sqA.str().c_str());
TH2* h_countD = (TH2*) f->Get(name_countD.str().c_str());
TH2* h_countA = (TH2*) f->Get(name_countA.str().c_str());

// Print existence check for each histogram
if (!h_ptD)    std::cerr << "   Missing: " << name_ptD.str()    << '\n';
if (!h_ptA)    std::cerr << "   Missing: " << name_ptA.str()    << '\n';
if (!h_sqptD)  std::cerr << "   Missing: " << name_sqD.str()    << '\n';
if (!h_sqptA)  std::cerr << "   Missing: " << name_sqA.str()    << '\n';
if (!h_countD) std::cerr << "   Missing: " << name_countD.str() << '\n';
if (!h_countA) std::cerr << "   Missing: " << name_countA.str() << '\n';

// Global missing check
if (!h_ptD || !h_ptA || !h_sqptD || !h_sqptA || !h_countD || !h_countA) {
    std::cerr << "[computeDptRegion4x4] Missing one or more histograms for region " << region << " and tag " << tag << '\n';
    f->Close();
    return;
}




    int Nphih = h_ptD->GetNbinsX();  // φh
    int Nz = h_ptD->GetNbinsY();    // z

    if (Nphih != 4 || Nz != 4) {
        std::cerr << "[computeDptRegion4x4] Unexpected binning: φh=" << Nphih << ", z=" << Nz << '\n';
    }

    for (int iphih = 1; iphih <= std::min(Nphih, 4); ++iphih) {
        for (int iz = 1; iz <= std::min(Nz, 4); ++iz) {
            double valD    = h_ptD->GetBinContent(iphih, iz);
            double valA    = h_ptA->GetBinContent(iphih, iz);
            double countD  = h_countD->GetBinContent(iphih, iz);
            double countA  = h_countA->GetBinContent(iphih, iz);
            double sqD     = h_sqptD->GetBinContent(iphih, iz);
            double sqA     = h_sqptA->GetBinContent(iphih, iz);

            double avgD = (countD > 0) ? valD / countD : 0.0;
            double avgA = (countA > 0) ? valA / countA : 0.0;
            double dpt = avgA - avgD;

            double varD = (countD > 0) ? sqD / countD - avgD * avgD : 0.0;
            double varA = (countA > 0) ? sqA / countA - avgA * avgA : 0.0;
            double errD = (countD > 0) ? std::sqrt(varD / countD) : 0.0;
            double errA = (countA > 0) ? std::sqrt(varA / countA) : 0.0;
            double errDpt = std::sqrt(errD * errD + errA * errA);
            matrix.val[iphih - 1][iz - 1] = dpt;
            matrix.err[iphih - 1][iz - 1] = errDpt;
        }
    }

    f->Close();
}

//drawing REGIONS//
void drawDptRegionPdf(const DptMatrixSimple& C2, const DptMatrixSimple& Cu, const DptMatrixSimple& Sn, const std::string& regionTag) {    const std::string refTag = "C2";  // ← pick any reliable target
    const std::string refFilePath = "DptinputFiles/Dhist_" + refTag + "_RGD.root";

    TFile* f = TFile::Open(refFilePath.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "[drawDptRegionPdf] Failed to open reference ROOT file: " << refFilePath << std::endl;
        return;
    }

    std::string refHistName = "region" + regionTag + "_D_" + refTag + "_RGD";
    TH2F* binref = (TH2F*) f->Get(refHistName.c_str());

    if (!binref) {
        std::cerr << "[drawDptRegionPdf] Reference histogram not found: " << refHistName << std::endl;
        f->Close();
        return;
    }

    std::cout << "[drawDptRegionPdf] Loaded reference histogram: " << refHistName << " from " << refFilePath << std::endl;


    std::string outPdf = "Dptvalues/DptRegion_" + regionTag + ".pdf";
    TCanvas canvas("c", "Dpt2triple-target", 1200, 800);
    canvas.Divide(2, 2);  // 4 phih bins per page

    const int Nphih = 4;  // fixed structure
    const int Nz = 4;

    for (int iphih = 0; iphih < Nphih; ++iphih) {
        canvas.cd(iphih + 1);

        auto* gC2 = new TGraphErrors();
        auto* gCu = new TGraphErrors();
        auto* gSn = new TGraphErrors();

        for (int iz = 0; iz < Nz; ++iz) {
            const double zVal = binref->GetYaxis()->GetBinCenter(iz + 1);  // Y axis = z

            double valC2 = C2.val[iphih][iz];
            double valCu = Cu.val[iphih][iz];
            double valSn = Sn.val[iphih][iz];
            double errC2 = C2.err[iphih][iz];
            double errCu = Cu.err[iphih][iz];
            double errSn = Sn.err[iphih][iz];

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
        mg->GetYaxis()->SetTitle("#Delta#LTp_{T}^{2}#GT");
        mg->GetYaxis()->SetRangeUser(-0.05, 0.15);

        TLegend* legend = new TLegend(0.15, 0.15, 0.35, 0.30);
        legend->SetTextSize(0.035);
        legend->AddEntry(gC2, "C", "lp");
        legend->AddEntry(gCu, "Cu", "lp");
        legend->AddEntry(gSn, "Sn", "lp");
        legend->Draw("same");

        TLine* line = new TLine(0.0, 0.0, 1.0, 0.0);  // horizontal 0-line
        line->SetLineStyle(2);
        line->Draw("same");

        TLatex text;
        text.SetTextSize(0.045);
        text.DrawLatexNDC(0.45, 0.22, Form("%.2f < #phi < %.2f",
            binref->GetXaxis()->GetBinLowEdge(iphih + 1),
            binref->GetXaxis()->GetBinUpEdge(iphih + 1)));

        TLatex watermark;
        watermark.SetTextSize(0.08);
        watermark.SetTextAngle(45);
        watermark.SetTextColorAlpha(kGray + 1, 0.3);
        watermark.SetNDC();
        watermark.SetTextAlign(22);
        watermark.DrawLatex(0.5, 0.5, "very preliminary");
    }

    canvas.Print((outPdf).c_str());
}


int main() {
    DptMatrixStruct dptC, dptCu, dptSn;
    TH3F* refHist = nullptr;

    computeDpt("DptinputFiles/Dhist_C2_RGD.root", "C2_RGD", dptC, refHist);
    computeDpt("DptinputFiles/Dhist_Cu_RGD.root", "Cu_RGD", dptCu, refHist);
    computeDpt("DptinputFiles/Dhist_Sn_RGD.root", "Sn_RGD", dptSn, refHist);
    TFile* f = TFile::Open("DptinputFiles/Dhist_Cu_RGD.root");
if (!f || f->IsZombie()) { std::cerr << "cannot open file\n"; return 1; }

THnSparseD *hCntD{nullptr}, *hCntA{nullptr};
f->GetObject("Q_x_phih_z_D_Cu_RGD",  hCntD);   // counts (D)
f->GetObject("Q_x_phih_z_A_Cu_RGD",  hCntA);   // counts (A)

inspectTHnSparseAxes(hCntD, "Counts-D  (Cu)");

inspectTHnSparseAxes(hCntA, "Counts-A  (Cu)");

checkTHnSparseContents("DptinputFiles/Dhist_Cu_RGD.root","Q_x_phih_z_D_Cu_RGD");

checkTHnSparseContents("DptinputFiles/Dhist_C2_RGD.root","Q_x_phih_z_A_C2_RGD");
checkTHnSparseContents("DptinputFiles/Dhist_Sn_RGD.root","Q_x_phih_z_A_Sn_RGD");

//if (!okD || !okA) {
//    std::cerr << "Histogram sanity check failed — aborting analysis\n";
//    return 1;
//}


//==REGIONS==//
DptMatrixSimple mat_C2_A, mat_Cu_A, mat_Sn_A;  // for region A
DptMatrixSimple mat_C2_B, mat_Cu_B, mat_Sn_B;  // for region B
DptMatrixSimple mat_C2_C, mat_Cu_C, mat_Sn_C;  // for region C
DptMatrixSimple mat_C2_D, mat_Cu_D, mat_Sn_D;  // for region D

computeDptRegion4x4("DptinputFiles/Dhist_C2_RGD.root", "C2", 'A', mat_C2_A);
computeDptRegion4x4("DptinputFiles/Dhist_Cu_RGD.root", "Cu", 'A', mat_Cu_A);
computeDptRegion4x4("DptinputFiles/Dhist_Sn_RGD.root", "Sn", 'A', mat_Sn_A);
computeDptRegion4x4("DptinputFiles/Dhist_C2_RGD.root", "C2", 'B', mat_C2_B);
computeDptRegion4x4("DptinputFiles/Dhist_Cu_RGD.root", "Cu", 'B', mat_Cu_B);
computeDptRegion4x4("DptinputFiles/Dhist_Sn_RGD.root", "Sn", 'B', mat_Sn_B);
computeDptRegion4x4("DptinputFiles/Dhist_C2_RGD.root", "C2", 'C', mat_C2_C);
computeDptRegion4x4("DptinputFiles/Dhist_Cu_RGD.root", "Cu", 'C', mat_Cu_C);
computeDptRegion4x4("DptinputFiles/Dhist_Sn_RGD.root", "Sn", 'C', mat_Sn_C);
computeDptRegion4x4("DptinputFiles/Dhist_C2_RGD.root", "C2", 'D', mat_C2_D);
computeDptRegion4x4("DptinputFiles/Dhist_Cu_RGD.root", "Cu", 'D', mat_Cu_D);
computeDptRegion4x4("DptinputFiles/Dhist_Sn_RGD.root", "Sn", 'D', mat_Sn_D);
    // Draw regions
drawDptRegionPdf(mat_C2_A, mat_Cu_A, mat_Sn_A, "A");
drawDptRegionPdf(mat_C2_B, mat_Cu_B, mat_Sn_B, "B");
drawDptRegionPdf(mat_C2_C, mat_Cu_C, mat_Sn_C, "C");
drawDptRegionPdf(mat_C2_D, mat_Cu_D, mat_Sn_D, "D");


// END REGIONS //


    //computeDpt4D_LoopsPhih("DptinputFiles/Dhist_C2_RGD.root", "C2_RGD");
    //computeDpt4D_LoopsPhih("DptinputFiles/Dhist_Cu_RGD.root", "Cu_RGD");
    //computeDpt4D_LoopsPhih("DptinputFiles/Dhist_Sn_RGD.root", "Sn_RGD");

    if (refHist) {
        drawDptPdf(dptC, dptCu, dptSn, refHist);
        //std::cout << "[plotDpt] Dpt PDF generated in ../Dptvalues/" << std::endl;
    } else {
        //std::cerr << "[plotDpt] No valid reference histogram.\n";
    }
    return 0;
}