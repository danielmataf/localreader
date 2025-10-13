#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <tuple>

// g++ targetcomp.cpp $(root-config --cflags --libs) -o targetcomp
// ./targetcomp

// -------- cut values (edit as you like) --------
double cut_loQ2       = 1.0;   double cut_hiQ2       = 100.0;
double cut_loW2       = 4.0; double cut_hiW2       = 100.0;
double cut_lonu       = 100.0; double cut_hinu       = 7.0;
double cut_lophih     = -100.0; double cut_hiphih     = -100.0;
double cut_loxB       = 100.0; double cut_hixB       = 100.0;
double cut_loy        = 100.0; double cut_hiy        = .7;
double cut_loz        = .3; double cut_hiz        = .7;
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

    std::vector<std::tuple<std::string, std::string>> targets = {
        //{"C2",  "~/Desktop/rootJUL8/REcutC2_test.root"},
        //{"Sn",  "~/Desktop/rootJUL8/REcutSn_test.root"},
        //{"Cu",  "~/Desktop/rootJUL8/REcutCu_test.root"},
        //{"LD2", "~/Desktop/rootJUL8/REcutLD2_test.root"}
        {"CxC",  "~/pass1CC.root"},
        {"Sn",  "~/pass1Sn.root"},
        {"Cu",  "~/pass1Cu.root"},
        {"LD2", "~/pass1LD2.root"}
   
    };

    std::vector<TFile*> rootFiles;
    for (const auto& [target, filePath] : targets) {
        TFile* file = new TFile(filePath.c_str(), "READ");
        if (!file->IsOpen()) {
            std::cerr << "Error: Cannot open file for " << target << std::endl;
            return;
        }
        rootFiles.push_back(file);
    }

    const std::vector<std::tuple<std::string, std::string>> histogramPairs = {
        {"Q2_", "Q^{2} [GeV^{2}]"}, {"W2_", "W^{2} [GeV^{2}]"}, {"nu_", "#nu [GeV]"},
        {"phih_", "#phi_{h} [deg]"}, {"xb_", "x_{B}"}, {"y_", "y"}, {"z_", "z"},
        {"targetVz_", "Target V_{z} [cm]"}, {"pt2_", "p_{T}^{2} [GeV^{2}]"},
        {"ptot_ele_", "p_{tot} Electron [GeV]"}, {"px_ele_", "p_{x} Electron [GeV]"},
        {"py_ele_", "p_{y} Electron [GeV]"}, {"pz_ele_", "p_{z} Electron [GeV]"},
        {"E_el", "E Electron [GeV]"}, {"E_pi", "E Pion [GeV]"},
        {"theta_el", "#theta Electron [deg]"}, {"phi_el", "#phi Electron [deg]"},
        {"theta_pi", "#theta Pion [deg]"}, {"phi_pi", "#phi Pion [deg]"},
        {"chi2_el_", "#chi^{2} Electron"}, {"chi2_pi_", "#chi^{2} Pion"}
    };

    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(0);

    TCanvas* pdfCanvas = new TCanvas("pdfCanvas", "Combined Histogram Comparison", 1000, 800);
    pdfCanvas->Divide(3, 3);
    int canvasIndex = 1, histCount = 0;
    pdfCanvas->Print("targetcomp.pdf[");

    for (const auto& [histBaseName, xAxisTitle] : histogramPairs) {
        std::vector<TH1F*> histograms;
        std::vector<int> colors = {kBlack, kOrange, kGreen, kBlue};

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
        TCanvas* canvas = new TCanvas(Form("ComparisonCanvas_%s", histBaseName.c_str()), "", 800, 600);
        histograms[0]->SetMaximum(maxVal * 1.1);
        histograms[0]->SetMinimum(0);
        histograms[0]->GetXaxis()->SetTitle(xAxisTitle.c_str());
        histograms[0]->Draw("hist");
        for (size_t i = 1; i < histograms.size(); ++i) histograms[i]->Draw("hist same");

        // two vertical lines
        TLine vlo_png(xlo, 0.0, xlo, maxVal * 1.1);
        vlo_png.SetLineStyle(3);   // dotted
        vlo_png.SetLineWidth(2);
        // vlo_png.SetLineColor(kGray+2); // optional

        TLine vhi_png(xhi, 0.0, xhi, maxVal * 1.1);
        vhi_png.SetLineStyle(2);   // dashed
        vhi_png.SetLineWidth(2);
        // vhi_png.SetLineColor(kGray+1); // optional

        vlo_png.Draw();
        vhi_png.Draw();

        TLegend* legend = new TLegend(0.72, 0.70, 0.95, 0.88);
        legend->SetTextSize(0.05);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->AddEntry(histograms[0], "C", "l");
        legend->AddEntry(histograms[1], "Sn", "l");
        legend->AddEntry(histograms[2], "Cu", "l");
        legend->AddEntry(histograms[3], "LD2", "l");
        legend->Draw();

        canvas->SaveAs(Form("targetcomp_%s.png", histBaseName.c_str()));

        // ---- PDF pad ----
        pdfCanvas->cd(canvasIndex);
        histograms[0]->Draw("hist");
        for (size_t i = 1; i < histograms.size(); ++i) histograms[i]->Draw("hist same");

        TLine vlo_pdf(xlo, 0.0, xlo, maxVal * 1.1); vlo_pdf.SetLineStyle(3); vlo_pdf.SetLineWidth(2); vlo_pdf.Draw();
        TLine vhi_pdf(xhi, 0.0, xhi, maxVal * 1.1); vhi_pdf.SetLineStyle(2); vhi_pdf.SetLineWidth(2); vhi_pdf.Draw();

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
    for (auto* f : rootFiles) f->Close();
}

int main() {
    CompareHistograms();
    return 0;
}
