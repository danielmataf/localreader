// emc_compare.cpp
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TAxis.h"

struct Group {
    std::string title;    // e.g. "Regions A–B"
    std::string tag;      // e.g. "AB"
    std::vector<int> idx; // 0-based indices into A..R
};

// ---- Fixed Y range via SetRangeUser
static const bool   kFixY = true;   // true -> lock Y to [0,1]
static const double kYMin = 0.0;
static const double kYMax = 1.0;

// Q2 LO/HI per region A..R (index 0..17). Use hi=-1 for open upper bound.
static const std::vector<double> q2Lo = {
    1.045, 1.000, 1.309, 1.104, 1.000, 1.490, 1.231, 1.000, 1.672,
    1.352, 1.000, 1.870, 1.479, 0.987, 2.065, 1.599, 1.233, 1.753
};
static const std::vector<double> q2Hi = {
    1.388, 1.309, 1.719, 1.490, 1.231, 2.115, 1.672, 1.352, 2.644,
    1.870, 1.479, 3.305, 2.065, 1.599, 4.759, 2.369, 1.775, -1.0
};

static std::string q2LabelForIndex(int idx) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3)
        << q2Lo[idx] << " < Q2 < ";
    if (q2Hi[idx] > 0.0) oss << q2Hi[idx];
    else                 oss << "...";
    return oss.str();
}

static TGraph* makeGraph(const std::vector<double>& data, const std::vector<int>& idxs, double xShift,
                         int color, int mstyle, double msize)
{
    std::vector<double> xs; xs.reserve(idxs.size());
    std::vector<double> ys; ys.reserve(idxs.size());
    for (size_t j = 0; j < idxs.size(); ++j) {
        int k = idxs[j];
        if (k >= 0 && k < (int)data.size() && data[k] >= 0.0) {
            xs.push_back((double)(j+1) + xShift);
            ys.push_back(data[k]);
        }
    }
    if (xs.empty()) return nullptr;

    auto* g = new TGraph((int)xs.size());
    for (int i = 0; i < (int)xs.size(); ++i) g->SetPoint(i, xs[i], ys[i]);
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(mstyle);
    g->SetMarkerSize(msize);
    g->SetLineWidth(0); // markers only (no connecting lines)
    return g;
}

static std::string xbTextForTag(const std::string& tag)
{
    if (tag == "AB")  return "x_{B} \\in [0.07, 0.10]";
    if (tag == "CDE") return "x_{B} \\in [0.10, 0.13]";
    if (tag == "FGH") return "x_{B} \\in [0.13, 0.16]";
    if (tag == "IJK") return "x_{B} \\in [0.16, 0.20]";
    if (tag == "LMN") return "x_{B} \\in [0.20, 0.25]";
    if (tag == "OPQ") return "x_{B} \\in [0.25, 0.36]";
    if (tag == "R")   return "x_{B} \\in [0.36, \\ldots]";
    return "";
}

static void drawGroup(const Group& G,
                      const std::vector<double>& rCxC,
                      const std::vector<double>& rSn,
                      const std::vector<double>& rCu,
                      bool firstPage, bool lastPage,
                      const std::string& pdfName)
{
    const int N = (int)G.idx.size();

    TCanvas* c = new TCanvas(("c_"+G.tag).c_str(), G.title.c_str(), 900, 600);
    gPad->SetGridx(true);
    gPad->SetGridy(true);

    // Empty title, no X-axis title, Y-axis title set
    TH1F* frame = new TH1F(("frame_"+G.tag).c_str(),
                           ";;Electron ratio (Target/LD2)",
                           N, 0.5, N + 0.5);
    frame->SetStats(0);

    // X-axis bin labels: Q2 intervals for each region in the group
    for (int i = 1; i <= N; ++i) {
        int k = G.idx[i-1]; // region index
        frame->GetXaxis()->SetBinLabel(i, q2LabelForIndex(k).c_str());
    }

    // Y-axis range
    if (kFixY) {
        frame->GetYaxis()->SetRangeUser(kYMin, kYMax);
    } else {
        double ymin = +1e9, ymax = -1e9;
        auto upd = [&](const std::vector<double>& v){
            for (int i = 0; i < N; ++i) {
                int k = G.idx[i];
                if (k >= 0 && k < (int)v.size() && v[k] >= 0.0) {
                    ymin = std::min(ymin, v[k]);
                    ymax = std::max(ymax, v[k]);
                }
            }
        };
        upd(rCxC); upd(rSn); upd(rCu);
        if (!(ymax > ymin)) { ymin = 0.4; ymax = 0.8; }
        const double pad = std::max(0.02, 0.10*(ymax - ymin));
        frame->GetYaxis()->SetRangeUser(std::max(0.0, ymin - pad), std::min(1.0, ymax + pad));
    }

    frame->GetXaxis()->SetLabelSize(0.045); // slightly smaller for longer labels
    frame->GetXaxis()->SetLabelOffset(0.01);
    frame->GetYaxis()->SetTitleOffset(1.1);
    frame->Draw();

    // Slight x-offsets so markers don't overlap
    TGraph* gCxC = makeGraph(rCxC, G.idx, -0.15, kBlue+1, 20, 1.3);
    TGraph* gSn  = makeGraph(rSn,  G.idx,  0.00, kRed+1,  21, 1.3);
    TGraph* gCu  = makeGraph(rCu,  G.idx, +0.15, kGreen+2,22, 1.3);

    if (gCxC) gCxC->Draw("P same");
    if (gSn)  gSn ->Draw("P same");
    if (gCu)  gCu ->Draw("P same");

    // Legend
    TLegend* leg = new TLegend(0.15, 0.78, 0.40, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    if (gCxC) leg->AddEntry(gCxC, "CxC / LD2", "p");
    if (gSn)  leg->AddEntry(gSn,  "Cu / LD2",  "p");
    if (gCu)  leg->AddEntry(gCu,  "Sn / LD2",  "p");
    leg->Draw();

    // x_B window label (top-right, via TLatex)
    const std::string xb = xbTextForTag(G.tag);
    if (!xb.empty()) {
        TLatex* lab = new TLatex();
        lab->SetNDC(true);
        lab->SetTextSize(0.045);
        lab->DrawLatex(0.62, 0.93, xb.c_str());
    }

    // Outputs
    c->SaveAs((std::string("EMC_") + G.tag + ".png").c_str());
    if (firstPage) c->SaveAs((pdfName + "(").c_str());
    else if (lastPage) c->SaveAs((pdfName + ")").c_str());
    else c->SaveAs(pdfName.c_str());
}

int main()
{
    gStyle->SetOptStat(0);

    // Ratios (A..R)
    std::vector<double> rCxC = {
        0.7753,0.7748,0.7952,0.7853,0.7688,0.7861,0.7756,0.7677,0.7846,
        0.7729,0.7648,0.7846,0.7733,0.7631,0.7805,0.7623,0.7585,0.7572
    };
    std::vector<double> rSn = {
        0.4791,0.4751,0.4849,0.4743,0.4594,0.4774,0.4693,0.4581,0.4750,
        0.4633,0.4513,0.4664,0.4555,0.4461,0.4515,0.4425,0.4391,0.4283
    };
    std::vector<double> rCu = {
        0.3260,0.3255,0.3367,0.3328,0.3250,0.3324,0.3311,0.3268,0.3341,
        0.3306,0.3310,0.3337,0.3309,0.3325,0.3258,0.3274,0.3262,0.3207
    };

    // Groupings
    std::vector<Group> groups = {
        {"Regions A–B", "AB",  {0,1}},
        {"Regions C–E", "CDE", {2,3,4}},
        {"Regions F–H", "FGH", {5,6,7}},
        {"Regions I–K", "IJK", {8,9,10}},
        {"Regions L–N", "LMN", {11,12,13}},
        {"Regions O–Q", "OPQ", {14,15,16}},
        {"Region R",    "R",   {17}}
    };

    const std::string pdfName = "EMC_Ratios_Compare.pdf";
    for (size_t i = 0; i < groups.size(); ++i) {
        const bool first = (i == 0);
        const bool last  = (i + 1 == groups.size());
        drawGroup(groups[i], rCxC, rSn, rCu, first, last, pdfName);
    }

    std::cout << "Wrote " << pdfName << " and PNGs for each group.\n";
    if (kFixY) std::cout << "Y-axis fixed via SetRangeUser(" << kYMin << ", " << kYMax << ").\n";
    return 0;
}
