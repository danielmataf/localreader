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
#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "constants.h"
#include <THnSparse.h>
#include <TStyle.h>


//compile and run: 
//g++ plotRatio.cpp -o plotRatio `root-config --cflags --libs`
//./plotRatio   -- (default) uses txt input --also w/ 1
//./plotRatio 2 -- uses root input
const int Rbin = Constants::Rbin_nu; //exporting value from csts
const int FIXED_SIZE = Constants::Rbin_nu;

// == REMEMBER ==
// 1. variables in 3D are nu, z, pt2 (last 2 are hadron only)
// 2. variabnles in 5D are xb, Q2, nu, z, pt2 (last two are hadron only )


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

//this is for 1D
struct RzData {
    std::vector<double> bins;
    std::vector<double> values;
    std::vector<double> errors;
};

RzData calcRz(const std::string& filePath, const std::string& tag) {
    RzData data;

    TFile* f = TFile::Open(filePath.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[calcRz] Cannot open file: " << filePath << std::endl;
        return data;
    }

    auto* h_z_D  = dynamic_cast<TH1F*>(f->Get(("z_D_new" + tag).c_str()));
    auto* h_z_A  = dynamic_cast<TH1F*>(f->Get(("z_A_new" + tag).c_str()));
    auto* h_nu_D = dynamic_cast<TH1F*>(f->Get(("nu_D_" + tag).c_str()));
    auto* h_nu_A = dynamic_cast<TH1F*>(f->Get(("nu_A_" + tag).c_str()));

    if (!h_z_D || !h_z_A || !h_nu_D || !h_nu_A) {
        std::cerr << "[calcRz] Missing histograms in file: " << filePath << std::endl;
        f->Close();
        return data;
    }

    double normD = h_nu_D->GetEntries();
    double normA = h_nu_A->GetEntries();
    int nbins = h_z_D->GetNbinsX();

    data.bins.resize(nbins);
    data.values.resize(nbins);
    data.errors.resize(nbins);

    for (int i = 1; i <= nbins; ++i) {
        double z = h_z_D->GetBinCenter(i);
        double d = h_z_D->GetBinContent(i);
        double a = h_z_A->GetBinContent(i);

        double ratio = 0.0;
        double error = 0.0;

        if (d > 0 && a > 0 && normD > 0 && normA > 0) {
            double normedD = d / normD;
            double normedA = a / normA;
            ratio = normedA / normedD;
            error = ratio * sqrt(1.0 / a + 1.0 / d + 1.0 / normA + 1.0 / normD);
        }

        data.bins[i - 1] = z;
        data.values[i - 1] = ratio;
        data.errors[i - 1] = error;
    }

    f->Close();
    return data;
}


void plotRzComparison(const RzData& rC, const RzData& rSn, const RzData& rCu, const std::string& outNameBase) {
    TCanvas* c = new TCanvas("cRz", "R(z) Comparison", 800, 600);

    TGraphErrors* gC  = new TGraphErrors(rC.bins.size());
    TGraphErrors* gSn = new TGraphErrors(rSn.bins.size());
    TGraphErrors* gCu = new TGraphErrors(rCu.bins.size());

    for (size_t i = 0; i < rC.bins.size(); ++i) {
        gC->SetPoint(i, rC.bins[i]+0.02, rC.values[i]);
        gC->SetPointError(i, 0.0, rC.errors[i]);
    }
    for (size_t i = 0; i < rSn.bins.size(); ++i) {
        gSn->SetPoint(i, rSn.bins[i], rSn.values[i]);
        gSn->SetPointError(i, 0.0, rSn.errors[i]);
    }
    for (size_t i = 0; i < rCu.bins.size(); ++i) {
        gCu->SetPoint(i, rCu.bins[i]+0.01, rCu.values[i]);
        gCu->SetPointError(i, 0.0, rCu.errors[i]);
    }

    gC->SetMarkerStyle(20);  gC->SetMarkerColor(kBlack);
    gSn->SetMarkerStyle(20); gSn->SetMarkerColor(kOrange);
    gCu->SetMarkerStyle(20); gCu->SetMarkerColor(kGreen);

    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gC, "P");
    mg->Add(gSn, "P");
    mg->Add(gCu, "P");

    mg->SetTitle("R(z) comparison");
    mg->GetXaxis()->SetTitle("z");
    mg->GetYaxis()->SetTitle("R(z)");
    mg->GetYaxis()->SetRangeUser(0, 2.0);
    mg->Draw("APE1");

    // Proper dotted line
    double xmin = mg->GetXaxis()->GetXmin();
    double xmax = mg->GetXaxis()->GetXmax();
    TLine* line = new TLine(xmin, 1.0, xmax, 1.0);
    line->SetLineStyle(2);
    line->Draw("same");

    TLegend* legend = new TLegend(0.7, 0.7, 0.88, 0.88);
    legend->AddEntry(gC, "C", "lp");
    legend->AddEntry(gSn, "Sn", "lp");
    legend->AddEntry(gCu, "Cu", "lp");
    legend->Draw("same");

    TLatex watermark;
    watermark.SetTextSize(0.08);
    watermark.SetTextAngle(45);
    watermark.SetTextColorAlpha(kGray + 1, 0.3);
    watermark.SetNDC();
    watermark.SetTextAlign(22);
    watermark.DrawLatex(0.5, 0.5, "preliminary pass0v11");

    //TLatex tag;
    //tag.SetTextSize(0.06);
    //tag.SetTextAngle(45);
    //tag.SetTextColorAlpha(kGray + 2, 0.3);
    //tag.SetNDC();
    //tag.SetTextAlign(22);
    //tag.DrawLatex(0.5, 0.4, "preliminary pass0v11");

    c->SaveAs((outNameBase + ".pdf").c_str());
    c->SaveAs((outNameBase + ".png").c_str());

    delete c;
}

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
        refHisto->SetDirectory(nullptr);  //to avoid deletion when file is closed
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


void computeSelfRatio3D(const std::string& fileC1, const std::string& fileC2, const std::string& tagC1, const std::string& tagC2, RatioMatrix& outMatrix, TH3F*& refHisto) {
    std::cout << "[computeSelfRatio3D] opening files..." << std::endl;

    TFile* f1 = TFile::Open(fileC1.c_str(), "READ");
    TFile* f2 = TFile::Open(fileC2.c_str(), "READ");

    if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie()) {
        std::cerr << "[computeSelfRatio3D] failed to open input files." << std::endl;
        return;
    }

    std::string nameH1 = "nu_z_pt2_A_" + tagC1;
    std::string nameH2 = "nu_z_pt2_A_" + tagC2;
    std::string nameNu1 = "nu_A_" + tagC1;
    std::string nameNu2 = "nu_A_" + tagC2;

    auto* hC1 = dynamic_cast<TH3F*>(f1->Get(nameH1.c_str()));
    auto* hC2 = dynamic_cast<TH3F*>(f2->Get(nameH2.c_str()));
    auto* hNuC1 = dynamic_cast<TH1F*>(f1->Get(nameNu1.c_str()));
    auto* hNuC2 = dynamic_cast<TH1F*>(f2->Get(nameNu2.c_str()));

    if (!hC1 || !hC2 || !hNuC1 || !hNuC2) {
        std::cerr << "[computeSelfRatio3D] histograms missing!" << std::endl;
        f1->Close(); f2->Close();
        return;
    }

    if (!refHisto) {
        refHisto = hC1;
        refHisto->SetDirectory(nullptr);
    }

    int NX = hC1->GetNbinsX();
    int NY = hC1->GetNbinsY();
    int NZ = hC1->GetNbinsZ();

    for (int x = 1; x <= NX; ++x) {
        double nuC1 = hNuC1->GetBinContent(x);
        double nuC2 = hNuC2->GetBinContent(x);
        for (int y = 1; y <= NY; ++y) {
            for (int z = 1; z <= NZ; ++z) {
                double valC1 = hC1->GetBinContent(x, y, z);
                double valC2 = hC2->GetBinContent(x, y, z);

                double normedC1 = (nuC1 > 0) ? valC1 / nuC1 : 0.0;
                double normedC2 = (nuC2 > 0) ? valC2 / nuC2 : 0.0;
                double R = (normedC2 > 0) ? normedC1 / normedC2 : 0.0;
                double err = (valC1 > 0 && valC2 > 0 && nuC1 > 0 && nuC2 > 0) ? R * std::sqrt(1.0 / valC1 + 1.0 / valC2 + 1.0 / nuC1 + 1.0 / nuC2) : 0.0;
                //std::cout << "[computeSelfRatio3D] R=" << R << ", err=" << err << std::endl;
                outMatrix.val[x - 1][y - 1][z - 1] = R;
                outMatrix.err[x - 1][y - 1][z - 1] = err;
            }
        }
    }

    f1->Close();
    f2->Close();
    std::cout << "[computeSelfRatio3D] Completed self-ratio computation.\n";
}


const int Rdim = 5;  //setting THnSparseD 5 dim



void compute5DRatio_withLoops(const std::string& fileName, const std::string& tag) {
    std::cout << "[compute5DRatio_withLoops] Attempting to open file: " << fileName << std::endl;
    TFile* f = new TFile(fileName.c_str(), "READ");
    //saniry checks 
    if (!f || f->IsZombie()) {
        std::cerr << "[compute5DRatio_withLoops] Failed to open file: " << fileName << std::endl;
        return;
    }

    std::string nameD5D = "Q2_xB_nu_pt2_z_D_" + tag;
    std::string nameA5D = "Q2_xB_nu_pt2_z_A_" + tag;
    std::string nameD3D = "Q2_xB_nu_D_" + tag;
    std::string nameA3D = "Q2_xB_nu_A_" + tag;

    auto* hD_5D = dynamic_cast<THnSparseD*>(f->Get(nameD5D.c_str()));
    auto* hA_5D = dynamic_cast<THnSparseD*>(f->Get(nameA5D.c_str()));
    auto* hD_3D = dynamic_cast<TH3F*>(f->Get(nameD3D.c_str()));
    auto* hA_3D = dynamic_cast<TH3F*>(f->Get(nameA3D.c_str()));

    if (!hD_5D || !hA_5D || !hD_3D || !hA_3D) {
        std::cerr << "[compute5DRatio_withLoops] One or more histograms missing. Exiting." << std::endl;
        f->Close();
        return;
    }

    // cloning a histo to replicate binning... but will be used for point value storage
    auto* hRatio = static_cast<THnSparseD*>(hD_5D->Clone(("h_5D_RatioLoop_" + tag).c_str()));
    hRatio->Reset();
    //getting all bins from each dimension
    const int nbinsQ2  = hD_5D->GetAxis(0)->GetNbins();
    const int nbinsxB  = hD_5D->GetAxis(1)->GetNbins();
    const int nbinsNu  = hD_5D->GetAxis(2)->GetNbins();
    const int nbinsPt2 = hD_5D->GetAxis(4)->GetNbins();
    const int nbinsZ   = hD_5D->GetAxis(3)->GetNbins();

    Int_t indices[5];  //indexing in 5dim   I guess prefered used is long64
    int Riszero = 0;
    int nbofvalues = 0;
    //nested looping on each dim
    for (int iQ2 = 1; iQ2 <= nbinsQ2; ++iQ2) {
        for (int ixB = 1; ixB <= nbinsxB; ++ixB) {
            for (int iNu = 1; iNu <= nbinsNu; ++iNu) {
                //stopped at 3dim for electrons (3D histograms) values for counts 
                double normD = hD_3D->GetBinContent(iQ2, ixB, iNu);
                double normA = hA_3D->GetBinContent(iQ2, ixB, iNu);
                //proceeding with hadron variables 
                for (int iZ = 1; iZ <= nbinsZ; ++iZ) {
                    for (int iPt2 = 1; iPt2 <= nbinsPt2; ++iPt2) {
                        nbofvalues++;
                        //substract 1 because index has a 0-start 
                        indices[0] = iQ2 - 1;
                        indices[1] = ixB - 1;
                        indices[2] = iNu - 1;
                        indices[3] = iZ  - 1;
                        indices[4] = iPt2 - 1;

                        double valD = hD_5D->GetBinContent(indices);
                        double valA = hA_5D->GetBinContent(indices);

                        double rD = (normD > 0) ? valD / normD : 0.0;
                        double rA = (normA > 0) ? valA / normA : 0.0;
                        double R  = (rD > 0) ? rA / rD : 0.0;
                        double err = (valD > 0 && valA > 0 && normD > 0 && normA > 0) ? R * std::sqrt(1.0 / valA + 1.0 / valD + 1.0 / normA + 1.0 / normD) : 0.0;
                        if (R == 0 ) {
                            Riszero++;
                        }
                        hRatio->SetBinContent(indices, R);
                        hRatio->SetBinError(indices, err);
                    }
                }
            }
        }
    }
    std::cout << "[compute5DRatio_withLoops] Number of zero ratios: " << Riszero << "out of " << nbofvalues << std::endl;
    TFile* outFile = new TFile(("ratio_loop_" + tag + ".root").c_str(), "RECREATE");
    hRatio->Write();
    outFile->Close();
    f->Close();

    std::cout << "[compute5DRatio_withLoops] Ratio histogram written to: ratio_loop_" << tag << ".root" << std::endl;
}


void compute5D_LoopsPhih(const std::string& fileName, const std::string& tag) {
    std::cout << "[compute5D_LoopsPhih] Attempting to open file: " << fileName << std::endl;
    TFile* f = new TFile(fileName.c_str(), "READ");
    
    if (!f || f->IsZombie()) {
        std::cerr << "[compute5D_LoopsPhih] Failed to open file: " << fileName << std::endl;
        return;
    }

    // Histogram names
    std::string nameD5D = "Q2_xB_phih_z_pt2_D_" + tag;
    std::string nameA5D = "Q2_xB_phih_z_pt2_A_" + tag;
    std::string nameD2D = "xB_Q2_D" + tag;
    std::string nameA2D = "xB_Q2_A" + tag;

    // Load histograms
    auto* hD_5D = dynamic_cast<THnSparseD*>(f->Get(nameD5D.c_str()));
    auto* hA_5D = dynamic_cast<THnSparseD*>(f->Get(nameA5D.c_str()));
    auto* hD_2D = dynamic_cast<TH2F*>(f->Get(nameD2D.c_str()));
    auto* hA_2D = dynamic_cast<TH2F*>(f->Get(nameA2D.c_str()));

    if (!hD_5D || !hA_5D || !hD_2D || !hA_2D) {
        std::cerr << "[compute5D_LoopsPhih] One or more histograms missing. Exiting." << std::endl;
        f->Close();
        return;
    }

    // Prepare output ratio histogram
    auto* hRatio = static_cast<THnSparseD*>(hD_5D->Clone(("h_5D_RatioLoop_" + tag).c_str()));
    hRatio->Reset();

    // Get bin numbers
    const int nbinsQ2   = hD_5D->GetAxis(0)->GetNbins();
    const int nbinsxB   = hD_5D->GetAxis(1)->GetNbins();
    const int nbinsPhih = hD_5D->GetAxis(2)->GetNbins();
    const int nbinsZ    = hD_5D->GetAxis(3)->GetNbins();
    const int nbinsPt2  = hD_5D->GetAxis(4)->GetNbins();

    Int_t indices[5];
    int Riszero = 0;
    int nbofvalues = 0;

    for (int iQ2 = 1; iQ2 <= nbinsQ2; ++iQ2) {
        for (int ixB = 1; ixB <= nbinsxB; ++ixB) {
            double normD = hD_2D->GetBinContent(ixB, iQ2);  // TH2F: x = xB, y = Q²
            double normA = hA_2D->GetBinContent(ixB, iQ2);

            for (int iPhih = 1; iPhih <= nbinsPhih; ++iPhih) {
                for (int iZ = 1; iZ <= nbinsZ; ++iZ) {
                    for (int iPt2 = 1; iPt2 <= nbinsPt2; ++iPt2) {
                        nbofvalues++;
                        indices[0] = iQ2   - 1;
                        indices[1] = ixB   - 1;
                        indices[2] = iPhih - 1;
                        indices[3] = iZ    - 1;
                        indices[4] = iPt2  - 1;

                        double valD = hD_5D->GetBinContent(indices);
                        double valA = hA_5D->GetBinContent(indices);

                        double rD = (normD > 0) ? valD / normD : 0.0;
                        double rA = (normA > 0) ? valA / normA : 0.0;
                        double R  = (rD > 0) ? rA / rD : 0.0;
                        double err = (valD > 0 && valA > 0 && normD > 0 && normA > 0)
                            ? R * std::sqrt(1.0 / valA + 1.0 / valD + 1.0 / normA + 1.0 / normD)
                            : 0.0;

                        if (R == 0.0) ++Riszero;

                        hRatio->SetBinContent(indices, R);
                        hRatio->SetBinError(indices, err);
                    }
                }
            }
        }
    }

    std::cout << "[compute5D_LoopsPhih] Number of zero ratios: " << Riszero << " out of " << nbofvalues << std::endl;

    TFile* outFile = new TFile(("ratio_loop_" + tag + ".root").c_str(), "RECREATE");
    hRatio->Write();
    outFile->Close();
    f->Close();

    std::cout << "[compute5D_LoopsPhih] Ratio histogram written to: ratio_loop_" << tag << ".root" << std::endl;
}


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
            Int_t nb = h->GetAxis(d)->GetNbins() + 2; 
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
void inspectTHnSparseAxes(const char* filename, const char* histoname) {
    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    auto* h = dynamic_cast<THnSparseD*>(f->Get(histoname));
    if (!h) {
        std::cerr << "Histogram " << histoname << " not found." << std::endl;
        f->Close();
        return;
    }

    int ndim = h->GetNdimensions();
    std::cout << "THnSparseD: " << histoname << " has " << ndim << " dimensions." << std::endl;

    for (int i = 0; i < ndim; ++i) {
        TAxis* axis = h->GetAxis(i);
        std::cout << "Axis " << i << ": " << axis->GetTitle()
                  << " | Nbins = " << axis->GetNbins()
                  << " | Range = [" << axis->GetXmin()
                  << ", " << axis->GetXmax() << "]" << std::endl;
    }

    f->Close();
}

void compute5DRatio(const std::string& fileName, const std::string& tag) {
    std::cout << "[compute5DRatio] Attempting to open file: " << fileName << std::endl;
    TFile* f = new TFile(fileName.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[compute5DRatio] Failed to open file: " << fileName << std::endl;
        return;
    }

std::string nameD5D = "Q2_xB_nu_pt2_z_D_" + tag;
std::string nameA5D = "Q2_xB_nu_pt2_z_A_" + tag;
std::string nameD3D = "Q2_xB_nu_D_" + tag;
std::string nameA3D = "Q2_xB_nu_A_" + tag;

    std::cout << "[compute5DRatio] Looking for histograms:" << std::endl;
    std::cout << "  " << nameD5D << std::endl;
    std::cout << "  " << nameA5D << std::endl;
    std::cout << "  " << nameD3D << std::endl;
    std::cout << "  " << nameA3D << std::endl;

    auto* hD_5D = dynamic_cast<THnSparseD*>(f->Get(nameD5D.c_str()));
    auto* hA_5D = dynamic_cast<THnSparseD*>(f->Get(nameA5D.c_str()));
    auto* hD_3D = dynamic_cast<TH3F*>(f->Get(nameD3D.c_str()));
    auto* hA_3D = dynamic_cast<TH3F*>(f->Get(nameA3D.c_str()));

    if (!hD_5D) std::cerr << "[compute5DRatio] Missing histogram: " << nameD5D << std::endl;
    if (!hA_5D) std::cerr << "[compute5DRatio] Missing histogram: " << nameA5D << std::endl;
    if (!hD_3D) std::cerr << "[compute5DRatio] Missing histogram: " << nameD3D << std::endl;
    if (!hA_3D) std::cerr << "[compute5DRatio] Missing histogram: " << nameA3D << std::endl;

    if (!hD_5D || !hA_5D || !hD_3D || !hA_3D) {
        std::cerr << "[compute5DRatio] One or more histograms missing. Exiting." << std::endl;
        f->Close();
        return;
    }

    auto* hRatio = static_cast<THnSparseD*>(hD_5D->Clone(("h_5D_Ratio_" + tag).c_str()));
    hRatio->Reset();
    int Riszero = 0;
    Long64_t nbins = hD_5D->GetNbins();
    for (Long64_t i = 0; i < nbins; ++i) {
        double valD = hD_5D->GetBinContent(i);
        double valA = hA_5D->GetBinContent(i);
        std::cout << "[compute5DRatio] Bin " << i << " | valD = " << valD << ", valA = " << valA << std::endl;
        if (valD == 0 && valA == 0) continue;

        Int_t indices[5];
        hD_5D->GetBinContent(i, indices); //not proper use of this but idea.  needs custom handling for bin indices
//
        int binQ2 = indices[0] + 1;
        int binxB = indices[1] + 1;
        int binnu = indices[2] + 1;
        std::cout << "[compute5DRatio] Bin indices: " << binQ2 << ", " << binxB << ", " << binnu << std::endl;

        double nuD = hD_3D->GetBinContent(binQ2, binxB, binnu);
        double nuA = hA_3D->GetBinContent(binQ2, binxB, binnu);
//
        double rD = (nuD > 0) ? valD / nuD : 0.0;
        double rA = (nuA > 0) ? valA / nuA : 0.0;
        double R = (rD > 0) ? rA / rD : 0.0;
        double err = (valD > 0 && valA > 0 && nuD > 0 && nuA > 0) ? R * std::sqrt(1.0 / valA + 1.0 / valD + 1.0 / nuA + 1.0 / nuD) : 0.0;
//
        std::cout << "[compute5DRatio] R = " << R << ", err = " << err << std::endl;
        if (R == 0) {
            Riszero++;
        }
        //hRatio->SetBinContent(i, R);
        //hRatio->SetBinError(i, err);
    }
    std::cout << "[compute5DRatio] Number of bins with R = 0: " << Riszero << "out of " << nbins << std::endl;
    //currently round half of the points are =0. This probably due to bad handling of binning.... Need s to be plotted to check. 
    //1500 points are enough for 1500/6 250 plots (?)  
//
    TFile* outFile = new TFile(("ratio_" + tag + ".root").c_str(), "RECREATE");
    //hRatio->Write();
    //hRatio handling for plotting is not properly working. Do not use a THnSparseD for this.Use instead a map/matrix?
    outFile->Close();
    f->Close();
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

void drawRatio5DfromTHnSparse(std::vector<THnSparseD*> ratios, const std::string& outDir) {
    if (ratios.size() != 3) {
        std::cerr << "[drawRatio5DfromTHnSparse] Expected 3 THnSparseD histograms (C, Cu, Sn)." << std::endl;
        return;
    }

    auto* hC  = ratios[0];
    auto* hCu = ratios[1];
    auto* hSn = ratios[2];

    const int nQ2  = hC->GetAxis(0)->GetNbins();
    const int nxB  = hC->GetAxis(1)->GetNbins();
    const int nNu  = hC->GetAxis(2)->GetNbins();
    const int nPt2 = hC->GetAxis(3)->GetNbins();
    const int nZ   = hC->GetAxis(4)->GetNbins();

    gStyle->SetOptStat(0);

    for (int iQ2 = 1; iQ2 <= nQ2; ++iQ2) {
        std::string outPdf = outDir + "/Ratio_Q2bin" + std::to_string(iQ2) + ".pdf";
        std::cout << "[drawRatio5DfromTHnSparse] Processing Q2 bin " << iQ2 << " with output PDF: " << outPdf << std::endl;
        for (int ixB = 1; ixB <= nxB; ++ixB) {
            for (int iNu = 1; iNu <= nNu; ++iNu) {
                TCanvas* canvas = new TCanvas("canvas", "5D Ratio Plots", 1200, 800);
                canvas->Divide(3, 2);

                for (int iPt2 = 1; iPt2 <= nPt2; ++iPt2) {
                    canvas->cd(iPt2);

                    TMultiGraph* mg = new TMultiGraph();
                    TGraphErrors* gC  = new TGraphErrors();
                    TGraphErrors* gCu = new TGraphErrors();
                    TGraphErrors* gSn = new TGraphErrors();

                    for (int iZ = 1; iZ <= nZ; ++iZ) {

                        Int_t idx[5] = {iQ2 - 1, ixB - 1, iNu - 1,  iZ-1 , iPt2-1 };
                    std::cout << "[drawRatio5DfromTHnSparse] : z bin value =" << hC->GetAxis(3)->GetBinCenter(iZ) << std::endl;

                        double zVal = hC->GetAxis(3)->GetBinCenter(iZ);
                        double valC  = hC->GetBinContent(idx);
                        double errC  = hC->GetBinError(idx);
                        double valCu = hCu->GetBinContent(idx);
                        double errCu = hCu->GetBinError(idx);
                        double valSn = hSn->GetBinContent(idx);
                        double errSn = hSn->GetBinError(idx);
                        std::cout << "[drawRatio5DfromTHnSparse] Bin (Q2=" << iQ2
                        << ", xB=" << ixB
                        << ", Nu=" << iNu
                        << ", pt2=" << iPt2
                        << ", z=" << iZ
                        << ") | C= " << valC << " ± " << errC
                        << ", Cu= " << valCu << " ± " << errCu
                        << ", Sn= " << valSn << " ± " << errSn << std::endl;

                        //ward for skipping empty points (useful? ) commented
                        //if (valC == 0 && valCu == 0 && valSn == 0) continue;    //ward

                        int point = iZ - 1;
                        gC->SetPoint(point, zVal, valC);
                        gC->SetPointError(point, 0.0, errC);
                        gCu->SetPoint(point, zVal + 0.01, valCu);
                        gCu->SetPointError(point, 0.0, errCu);
                        gSn->SetPoint(point, zVal + 0.02, valSn);
                        gSn->SetPointError(point, 0.0, errSn);
                    }

                    gC->SetMarkerStyle(20);  gC->SetMarkerColor(kBlack);
                    gCu->SetMarkerStyle(20); gCu->SetMarkerColor(kGreen);
                    gSn->SetMarkerStyle(20); gSn->SetMarkerColor(kOrange);
                    mg->Add(gC); mg->Add(gCu); mg->Add(gSn);
                    mg->Draw("APE1");
                    mg->GetXaxis()->SetTitle("z");
                    mg->GetYaxis()->SetTitle("R");
                    mg->GetYaxis()->SetRangeUser(0.0, 1.5);

                    // Annotate bin ranges
                    TLatex label;
                    label.SetTextSize(0.03);
                    label.DrawLatexNDC(0.12, 0.88, Form("%.2f < Q^{2} < %.2f", hC->GetAxis(0)->GetBinLowEdge(iQ2), hC->GetAxis(0)->GetBinUpEdge(iQ2)));
                    label.DrawLatexNDC(0.12, 0.83, Form("%.2f < x_{B} < %.2f", hC->GetAxis(1)->GetBinLowEdge(ixB), hC->GetAxis(1)->GetBinUpEdge(ixB)));
                    label.DrawLatexNDC(0.12, 0.78, Form("%.2f < #nu < %.2f", hC->GetAxis(2)->GetBinLowEdge(iNu), hC->GetAxis(2)->GetBinUpEdge(iNu)));
                    label.DrawLatexNDC(0.12, 0.73, Form("%.2f < p_{T}^{2} < %.2f", hC->GetAxis(4)->GetBinLowEdge(iPt2), hC->GetAxis(4)->GetBinUpEdge(iPt2)));

                    // Reference line
                    TLine* refLine = new TLine(hC->GetAxis(3)->GetXmin(), 1.0, hC->GetAxis(3)->GetXmax(), 1.0);
                    refLine->SetLineStyle(2);
                    refLine->SetLineColor(kGray + 2);
                    refLine->Draw("same");

                    // Legend
                    TLegend* legend = new TLegend(0.12, 0.12, 0.35, 0.27);
                    legend->SetTextSize(0.025);
                    legend->AddEntry(gC, "C", "lp");
                    legend->AddEntry(gCu, "Cu", "lp");
                    legend->AddEntry(gSn, "Sn", "lp");
                    legend->Draw("same");

                    // Watermark
                    TLatex watermark;
                    watermark.SetTextSize(0.06);
                    watermark.SetTextAngle(45);
                    watermark.SetTextColorAlpha(kGray + 1, 0.25);
                    watermark.SetNDC();
                    watermark.SetTextAlign(22);
                    watermark.DrawLatex(0.5, 0.5, "very preliminary");
                }

                // Save canvas to PDF (handle open/close)
                if (ixB == 1 && iNu == 1)
                    canvas->Print((outPdf + "(").c_str());  // Open multi-page
                else if (ixB == nxB && iNu == nNu)
                    canvas->Print((outPdf + ")").c_str());  // Close multi-page
                else
                    canvas->Print(outPdf.c_str());          // Add page

                delete canvas;
            }
        }
    }

    std::cout << "[drawRatio5DfromTHnSparse] PDF output: 1 file per Q² bin, 36 pages each (6 xB × 6 ν).\n";
}


int main() {
    RzData rzC  = calcRz("RinputFiles/Rhist_C2_RGD.root", "C2_RGD");
    RzData rzSn = calcRz("RinputFiles/Rhist_Sn_RGD.root", "Sn_RGD");
    RzData rzCu = calcRz("RinputFiles/Rhist_Cu_RGD.root", "Cu_RGD");

    plotRzComparison(rzC, rzSn, rzCu, "Rvalues/Rz_comparison");

    RatioMatrix rC, rCu, rSn;
    TH3F *refHist = nullptr;
    computeSelfRatio3D("RinputFiles/Rhist_C1_RGD.root", "RinputFiles/Rhist_C2_RGD.root", "C1_RGD", "C2_RGD", rC, refHist); //change this to C1 once you have C1 file 
    computeRatio("RinputFiles/Rhist_C2_RGD.root", "C2_RGD", rC, refHist);
    computeRatio("RinputFiles/Rhist_Cu_RGD.root", "Cu_RGD", rCu, refHist);
    computeRatio("RinputFiles/Rhist_Sn_RGD.root", "Sn_RGD", rSn, refHist);
    if (refHist) {
        drawRatioPdf(rC, rCu, rSn, refHist);
        std::cout << "[plotRatio] R PDF generated in ../Rvalues/" << std::endl;
    } else {
        std::cerr << "[plotRatio] Aborting plot: no valid reference histogram loaded.\n";
    }
//    checkTHnSparseContents("RinputFiles/Rhist_C2_RGD.root", "Q2_xB_nu_pt2_z_D_C2_RGD");
//    checkTHnSparseContents("RinputFiles/Rhist_Cu_RGD.root", "Q2_xB_nu_pt2_z_D_Cu_RGD");
//    checkTHnSparseContents("RinputFiles/Rhist_Sn_RGD.root", "Q2_xB_nu_pt2_z_D_Sn_RGD");
//    inspectTHnSparseAxes("RinputFiles/Rhist_C2_RGD.root", "Q2_xB_nu_pt2_z_D_C2_RGD");
    //compute5DRatio("RinputFiles/Rhist_C2_RGD.root", "C2_RGD");
    //compute5DRatio("RinputFiles/Rhist_Cu_RGD.root", "Cu_RGD");
    //compute5DRatio("RinputFiles/Rhist_Sn_RGD.root", "Sn_RGD");
    
    
    //compute5DRatio_withLoops("RinputFiles/Rhist_C2_RGD.root", "C2_RGD");
    //compute5DRatio_withLoops("RinputFiles/Rhist_Cu_RGD.root", "Cu_RGD");
    //compute5DRatio_withLoops("RinputFiles/Rhist_Sn_RGD.root", "Sn_RGD");
    compute5D_LoopsPhih("RinputFiles/Rhist_C2_RGD.root", "C2_RGD");
    compute5D_LoopsPhih("RinputFiles/Rhist_Cu_RGD.root", "Cu_RGD");
    compute5D_LoopsPhih("RinputFiles/Rhist_Sn_RGD.root", "Sn_RGD");

    //TFile* fC  = TFile::Open("ratio_loop_C2_RGD.root");
    //TFile* fCu = TFile::Open("ratio_loop_Cu_RGD.root");
    //TFile* fSn = TFile::Open("ratio_loop_Sn_RGD.root");
//
    //auto* hC  = (THnSparseD*)fC->Get("h_5D_RatioLoop_C2_RGD");
    //auto* hCu = (THnSparseD*)fCu->Get("h_5D_RatioLoop_Cu_RGD");
    //auto* hSn = (THnSparseD*)fSn->Get("h_5D_RatioLoop_Sn_RGD");
//
    //if (!hC || !hCu || !hSn) {
    //    std::cerr << "ERROR: One or more THnSparse histograms could not be loaded. Check names!" << std::endl;
    //    return 1;
    //}
//
    //std::vector<THnSparseD*> ratios = {hC, hCu, hSn};
    //drawRatio5DfromTHnSparse(ratios, "Rvalues");
//



    return 0;
}
