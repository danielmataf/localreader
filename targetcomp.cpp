#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TPad.h>

#include <unordered_map>
#include <iostream>
#include <vector>
#include <tuple>
#include <string>

// g++ targetcomp.cpp $(root-config --cflags --libs) -o targetcomp
// ./targetcomp

// -------- cut values (edit as you like) --------
double cut_loQ2       = 1.0;   double cut_hiQ2       = 100.0;
double cut_loW2       = 4.0;   double cut_hiW2       = 100.0;
double cut_lonu       = 100.0; double cut_hinu       = 7.0;
double cut_lophih     = -100.0; double cut_hiphih    = -100.0;
double cut_loxB       = 100.0; double cut_hixB       = 100.0;
double cut_loy        = 100.0; double cut_hiy        = .7;
double cut_loz        = .3;    double cut_hiz        = .7;
double cut_lotargetVz = 100.0; double cut_hitargetVz = 100.0;
double cut_lopt2      = 100.0; double cut_hipt2      = 1.2;

double cut_loptot_el  = 100.0; double cut_hiptot_el  = 100.0;
double cut_lopx_el    = 100.0; double cut_hipx_el    = 100.0;
double cut_lopy_el    = 100.0; double cut_hipy_el    = 100.0;
double cut_lopz_el    = 100.0; double cut_hipz_el    = 100.0;
double cut_loE_el     = 100.0; double cut_hiE_el     = 100.0;
double cut_loE_pi     = 100.0; double cut_hiE_pi     = 100.0;
double cut_lotheta_el = 100.0; double cut_hitheta_el = 100.0;
double cut_lophi_el   = 100.0; double cut_hiphi_el   = 100.0;
double cut_lotheta_pi = 100.0; double cut_hitheta_pi = 100.0;
double cut_lophi_pi   = 100.0; double cut_hiphi_pi   = 100.0;
double cut_lochi2_el  = 100.0; double cut_hichi2_el  = 100.0;
double cut_lochi2_pi  = 100.0; double cut_hichi2_pi  = 100.0;

static void SetNiceStyle() {
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
}

void CompareHistograms() {
    // map histogram base -> low/high cut value
    const std::unordered_map<std::string, double> cutLo = {
        {"Q2_", cut_loQ2}, {"W2_", cut_loW2}, {"nu_", cut_lonu}, {"phih_", cut_lophih},
        {"xb_", cut_loxB}, {"y_", cut_loy}, {"z_", cut_loz}, {"targetVz_", cut_lotargetVz},
        {"pt2_", cut_lopt2},
        {"ptot_ele_", cut_loptot_el}, {"px_ele_", cut_lopx_el}, {"py_ele_", cut_lopy_el},
        {"pz_ele_", cut_lopz_el}, {"E_el", cut_loE_el}, {"E_pi", cut_loE_pi},
        {"theta_el", cut_lotheta_el}, {"phi_el", cut_lophi_el},
        {"theta_pi", cut_lotheta_pi}, {"phi_pi", cut_lophi_pi},
        {"chi2_el_", cut_lochi2_el}, {"chi2_pi_", cut_lochi2_pi}
    };
    const std::unordered_map<std::string, double> cutHi = {
        {"Q2_", cut_hiQ2}, {"W2_", cut_hiW2}, {"nu_", cut_hinu}, {"phih_", cut_hiphih},
        {"xb_", cut_hixB}, {"y_", cut_hiy}, {"z_", cut_hiz}, {"targetVz_", cut_hitargetVz},
        {"pt2_", cut_hipt2},
        {"ptot_ele_", cut_hiptot_el}, {"px_ele_", cut_hipx_el}, {"py_ele_", cut_hipy_el},
        {"pz_ele_", cut_hipz_el}, {"E_el", cut_hiE_el}, {"E_pi", cut_hiE_pi},
        {"theta_el", cut_hitheta_el}, {"phi_el", cut_hiphi_el},
        {"theta_pi", cut_hitheta_pi}, {"phi_pi", cut_hiphi_pi},
        {"chi2_el_", cut_hichi2_el}, {"chi2_pi_", cut_hichi2_pi}
    };

    // (tag, path)  — order sets legend order
    std::vector<std::tuple<std::string, std::string>> targets = {
        {"C2",  "~/dump/pass1C2.root"},
        {"Sn",  "~/dump/pass1Sn.root"},
        {"Cu",  "~/dump/pass1Cu.root"},
        {"LD2", "~/dump/pass1LD2.root"}
    };

    // open files
    std::vector<TFile*> rootFiles;
    for (const auto& [target, filePath] : targets) {
        TFile* file = new TFile(filePath.c_str(), "READ");
        if (!file || !file->IsOpen()) {
            std::cerr << "Error: Cannot open file for " << target << std::endl;
            return;
        }
        rootFiles.push_back(file);
    }

    // -------------- 1D overlays (unchanged in spirit, styled to match finComp) --------------
    SetNiceStyle();

    const std::vector<std::tuple<std::string, std::string>> histogramPairs = {
        {"Q2_", "Q^{2} [GeV^{2}]"}, {"W2_", "W^{2} [GeV^{2}]"}, {"nu_", "#nu [GeV]"},
        {"phih_", "#phi_{h} [deg]"}, {"xb_", "x_{B}"}, {"y_", "y"}, {"z_", "z"},
        {"targetVz_", "V_{z} [cm]"}, {"pt2_", "p_{T}^{2} [GeV^{2}]"},
        {"ptot_ele_", "p_{tot} (e^{-}) [GeV]"}, {"px_ele_", "p_{x} (e^{-}) [GeV]"},
        {"py_ele_", "p_{y} (e^{-}) [GeV]"}, {"pz_ele_", "p_{z} (e^{-}) [GeV]"},
        {"E_el", "E (e^{-}) [GeV]"}, {"E_pi", "E (#pi^{+}) [GeV]"},
        {"theta_el", "#theta (e^{-}) [deg]"}, {"phi_el", "#phi (e^{-}) [deg]"},
        {"theta_pi", "#theta (#pi^{+}) [deg]"}, {"phi_pi", "#phi (#pi^{+}) [deg]"},
        {"chi2_el_", "#chi^{2} (e^{-})"}, {"chi2_pi_", "#chi^{2} (#pi^{+})"}
    };

    // color order similar to finComp overlays
    std::vector<int> colors = { kBlack, kOrange+7, kGreen+2, kBlue+1 };

    TCanvas* pdfCanvas = new TCanvas("pdfCanvas", "Combined Histogram Comparison", 1000, 800);
    pdfCanvas->Divide(3, 3);
    int canvasIndex = 1, histCount = 0;
    pdfCanvas->Print("targetcomp.pdf[");

    for (const auto& [histBaseName, xAxisTitle] : histogramPairs) {
        std::vector<TH1F*> histograms;

        int idx = 0;
        for (const auto& [target, filePath] : targets) {
            TFile* file = rootFiles[idx];
            std::string histName = histBaseName + target + "_RGD";
            TH1F* hist = dynamic_cast<TH1F*>(file->Get(histName.c_str()));
            if (!hist) {
                std::cerr << "Error: Cannot retrieve histogram " << histName << std::endl;
                ++idx; continue;
            }
            hist->SetLineColor(colors[idx]);
            hist->SetStats(0);
            histograms.push_back(hist);
            ++idx;
        }
        if (histograms.size() < targets.size()) continue;

        // Normalize & max
        double maxVal = 0;
        for (auto* h : histograms) {
            double integral = h->Integral();
            if (integral != 0) h->Scale(1.0 / integral);
            if (h->GetMaximum() > maxVal) maxVal = h->GetMaximum();
        }

        // fetch cuts
        double xlo = 100.0, xhi = 100.0;
        if (auto it = cutLo.find(histBaseName); it != cutLo.end()) xlo = it->second;
        if (auto it = cutHi.find(histBaseName); it != cutHi.end()) xhi = it->second;

        // ---- PNG canvas ----
        TCanvas* canvas = new TCanvas(Form("ComparisonCanvas_%s", histBaseName.c_str()), "", 1000, 700);
        gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.12);

        histograms[0]->SetMaximum(maxVal * 1.15);
        histograms[0]->SetMinimum(0);
        histograms[0]->GetXaxis()->SetTitle(xAxisTitle.c_str());
        histograms[0]->GetYaxis()->SetTitle("Normalized events");
        histograms[0]->SetLineWidth(2);
        histograms[0]->Draw("hist");
        for (size_t i = 1; i < histograms.size(); ++i) {
            histograms[i]->SetLineWidth(2);
            histograms[i]->Draw("hist same");
        }

        // dashed vertical lines
        if (xlo != 100.0 && xlo != -100.0) { TLine* l = new TLine(xlo, 0.0, xlo, maxVal*1.15); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw("same"); }
        if (xhi != 100.0 && xhi != -100.0) { TLine* r = new TLine(xhi, 0.0, xhi, maxVal*1.15); r->SetLineStyle(2); r->SetLineWidth(3); r->Draw("same"); }

        TLegend* legend = new TLegend(0.70, 0.70, 0.93, 0.90);
        legend->SetTextSize(0.045);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->AddEntry(histograms[0], "C",   "l");
        legend->AddEntry(histograms[1], "Sn",  "l");
        legend->AddEntry(histograms[2], "Cu",  "l");
        legend->AddEntry(histograms[3], "LD2", "l");
        legend->Draw();

        canvas->SaveAs(Form("targetcomp_%s.png", histBaseName.c_str()));

        // ---- PDF pad ----
        pdfCanvas->cd(canvasIndex);
        gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.12);
        histograms[0]->Draw("hist");
        for (size_t i = 1; i < histograms.size(); ++i) histograms[i]->Draw("hist same");
        if (xlo != 100.0 && xlo != -100.0) { TLine l(xlo, 0.0, xlo, maxVal*1.15); l.SetLineStyle(2); l.SetLineWidth(3); l.Draw(); }
        if (xhi != 100.0 && xhi != -100.0) { TLine r(xhi, 0.0, xhi, maxVal*1.15); r.SetLineStyle(2); r.SetLineWidth(3); r.Draw(); }
        legend->Draw();

        ++canvasIndex; ++histCount;
        if (canvasIndex > 9) {
            pdfCanvas->Print("targetcomp.pdf");
            pdfCanvas->Clear();
            pdfCanvas->Divide(3, 3);
            canvasIndex = 1;
        }
        delete canvas;
    }

    if (histCount % 9 != 0) pdfCanvas->Print("targetcomp.pdf");
    pdfCanvas->Print("targetcomp.pdf]");
    delete pdfCanvas;

    // -------------- NEW: 2D plots (one PNG per target for each base) --------------
    // List the 2D bases you want to draw. Names are "<base><TAG>_RGD" (e.g., "pt2z_Cu_RGD").
    struct Var2D { std::string base; std::string title; std::string xlab; std::string ylab; };
    const std::vector<Var2D> vars2D = {
        {"pt2z_", "p_{T}^{2} vs z",       "z",                 "p_{T}^{2} [GeV^{2}]"},
        {"xQ2_",  "Q^{2} vs x_{B}",       "x_{B}",             "Q^{2} [GeV^{2}]"}
    };

    for (size_t i = 0; i < targets.size(); ++i) {
        const auto& [tag, path] = targets[i];
        TFile* file = rootFiles[i];

        for (const auto& v2 : vars2D) {
            const std::string hname = v2.base + tag + "_RGD";
            TH2F* h2 = dynamic_cast<TH2F*>(file->Get(hname.c_str()));
            if (!h2) {
                std::cerr << "Warning: missing 2D histogram " << hname << " in " << path << "\n";
                continue;
            }

            TCanvas c2(Form("c2_%s_%s", v2.base.c_str(), tag.c_str()), "", 1100, 800);
            gPad->SetRightMargin(0.12);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);

            // Title shows variable + target (ASCII hyphen to avoid encoding issues)
            h2->SetTitle(Form("%s - %s", v2.title.c_str(), tag.c_str()));
            h2->GetXaxis()->SetTitle(v2.xlab.c_str());
            h2->GetYaxis()->SetTitle(v2.ylab.c_str());
            h2->SetContour(60);
            h2->SetStats(0);
            h2->Draw("COLZ");
            //gPad->Update();
            //if (auto* st = dynamic_cast<TPaveStats*>(gPad->GetPrimitive("stats"))) st->Delete();
            c2.SaveAs(Form("targetcomp2D_%s%s.png", v2.base.c_str(), tag.c_str()));
        }
    }

    // close files
    for (auto* f : rootFiles) { f->Close(); delete f; }
}

int main() {
    CompareHistograms();
    return 0;
}
