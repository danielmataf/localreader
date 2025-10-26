// finComp.cpp
// g++ finComp.cpp $(root-config --cflags --libs) -o finComp && ./finComp
// 1D: overlays across targets (normalized) with dashed cut lines.
// 2D: both a 2x2 panel (Cu, CxC, LD2, Sn) and separate full-size plots per target.

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TClass.h>
#include <TColor.h>
#include <TPad.h>

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <array>
#include <algorithm>

// -------------------- CUTS (100.0 means "don't draw") --------------------
// Physics cuts for guide lines; adjust as needed.
double cut_loQ2       = 1.0;   double cut_hiQ2       = 100.0;
double cut_loW2       = 4.0;   double cut_hiW2       = 100.0;
double cut_lonu       = 100.0; double cut_hinu       = 7.0;

// DO NOT draw Vz window (per request)
double cut_lotargetVz = 100.0; double cut_hitargetVz = 100.0;

// pt2 upper guide
double cut_lopt2      = 100.0; double cut_hipt2      = 1.2;

// NEW: y upper guide at 0.75
double cut_loy        = 0.25; double cut_hiy        = 0.75;

// NEW: z two-sided window 0.30–0.70
double cut_loz        = 0.30;  double cut_hiz        = 0.70;

// -------------------- TARGETS (legend order & colors) --------------------
struct Target {
    std::string tag;     // Cu, CxC, LD2, Sn
    std::string label;   // legend label
    int         color;   // ROOT color
    std::string path;    // file path
};

static const std::vector<Target> kTargets = {
    {"Cu",  "Cu",        kGreen+2,  "~/fin_Cu_test.root"},
    {"CxC", "C (CxC)",   kBlack,    "~/fin_CxC_test.root"},
    {"LD2", "LD2",       kBlue+1,   "~/fin_LD2_test.root"},
    {"Sn",  "Sn",        kOrange+7, "~/fin_Sn_test.root"}
};

// -------------------- 1D variables --------------------
struct Var { std::string base; std::string xtitle; std::string title; };
static const std::vector<Var> kVars = {
    {"Q2",  "Q^{2} [GeV^{2}]",     "Q^{2}"},
    {"Vz",  "V_{z} [cm]",          "V_{z}"},
    {"W2",  "W^{2} [GeV^{2}]",     "W^{2}"},
    {"nu",  "#nu [GeV]",           "#nu"},
    {"pt2", "p_{T}^{2} [GeV^{2}]", "p_{T}^{2}"},
    {"y",   "y",                   "y"},            // with line at 0.75
    {"z",   "z",                   "z"}             // with lines at 0.30 and 0.70
    // If you add xB later:
    // {"xb", "x_{B}",               "x_{B}"}
};

// -------------------- 2D variables (unchanged) --------------------
struct Var2D { std::string base; std::string title; };
static const std::vector<Var2D> kVars2D = {
    {"sampl_el",  "Sampling fraction vs p (e-)"},
    {"beta_p",    "#beta vs p"},
    {"m2_p",      "m^{2} vs p"},
    {"pt2_vz_pi", "p_{T}^{2} vs V_{z} (#pi)"},
};
static const std::vector<std::string> kStates2D = {"pre", "post", "posPID", "posKIN"}; // tried in order

// -------------------- helpers --------------------
static std::pair<double,double> cutsFor(const std::string& base) {
    if (base == "Q2")  return {cut_loQ2,       cut_hiQ2};
    if (base == "W2")  return {cut_loW2,       cut_hiW2};
    if (base == "nu")  return {cut_lonu,       cut_hinu};
    if (base == "Vz")  return {cut_lotargetVz, cut_hitargetVz}; // disabled (100,100)
    if (base == "pt2") return {cut_lopt2,      cut_hipt2};
    if (base == "y")   return {cut_loy,        cut_hiy};        // 0.75 upper line
    if (base == "z")   return {cut_loz,        cut_hiz};        // 0.30 & 0.70
    return {100.0, 100.0};
}

static void setStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
}

// -------------------- HISTO LISTING (stdout) --------------------
static void ListHistosRecursive(TDirectory* dir,
                                std::vector<std::pair<std::string,std::string>>& out,
                                const std::string& prefix = "")
{
    TIter next(dir->GetListOfKeys());
    while (TKey* key = (TKey*)next()) {
        const char* cname = key->GetClassName();
        const char* kname = key->GetName();
        TClass* cls = gROOT->GetClass(cname);
        if (!cls) continue;

        if (cls->InheritsFrom(TDirectory::Class())) {
            if (TDirectory* sub = dir->GetDirectory(kname)) {
                ListHistosRecursive(sub, out, prefix + kname + "/");
            }
        } else if (cls->InheritsFrom(TH1::Class())) {
            std::string kind = "TH1";
            if (cls->InheritsFrom(TH2::Class())) kind = "TH2";
            out.emplace_back(kind, prefix + kname);
        }
    }
}

static void DumpHistListToStdout(TFile* f, const std::string& tag, const std::string& path) {
    std::vector<std::pair<std::string,std::string>> names; // (kind, name)
    ListHistosRecursive(f, names, "");
    std::sort(names.begin(), names.end(),
              [](auto& a, auto& b){ return a.second < b.second; });

    std::cout << "\n===== " << tag << " : " << path << " =====\n";
    std::cout << "Found " << names.size() << " histograms (TH1/TH2):\n";
    for (auto& kv : names) std::cout << "  [" << kv.first << "] " << kv.second << "\n";
    if (names.empty()) std::cout << "  (no histograms found)\n";
}

// -------------------- fetchers --------------------
static TH1* fetchClone(TFile* f, const std::string& name) {
    TObject* obj = f->Get(name.c_str());
    if (!obj || !obj->InheritsFrom(TH1::Class())) return nullptr;
    TH1* h = static_cast<TH1*>(obj);
    TH1* c = static_cast<TH1*>(h->Clone((name + "_cl").c_str()));
    c->SetDirectory(nullptr);
    return c;
}

static TH1* getHistFlexible(TFile* f,
                            const std::string& base,
                            const std::string& state,
                            const std::string& tag)
{
    // synonyms for base
    std::vector<std::string> bases = { base };
    if (base == "Vz")  bases.push_back("targetVz");
    if (base == "pt2") bases.push_back("pT2");
    if (base == "xb")  bases.push_back("xB");

    for (const auto& b : bases) {
        std::array<std::string,4> patterns = {
            b + "_" + state + "_" + tag + "_RGD",
            b + "_" + state + "_" + tag,
            b + "_" + state,
            "h_" + b + "_" + state
        };
        for (const auto& n : patterns) {
            if (TH1* h = fetchClone(f, n)) return h;
        }
    }
    return nullptr;
}

// TH2
static TH2* fetchClone2D(TFile* f, const std::string& name) {
    TObject* obj = f->Get(name.c_str());
    if (!obj || !obj->InheritsFrom(TH2::Class())) return nullptr;
    TH2* h = static_cast<TH2*>(obj);
    TH2* c = static_cast<TH2*>(h->Clone((name + "_cl2").c_str()));
    c->SetDirectory(nullptr);
    return c;
}

static TH2* getHist2DFlexible(TFile* f,
                              const std::string& base,
                              const std::string& state,
                              const std::string& tag)
{
    // with tag
    std::array<std::string,4> withTag = {
        base + "_" + state + "_" + tag + "_RGD",
        base + "_" + state + "_" + tag,
        "h_"  + base + "_" + state + "_" + tag + "_RGD",
        "h_"  + base + "_" + state + "_" + tag
    };
    for (const auto& n : withTag) if (TH2* h = fetchClone2D(f, n)) return h;

    // without tag
    std::array<std::string,4> noTag = {
        base + "_" + state + "_RGD",
        base + "_" + state,
        "h_"  + base + "_" + state + "_RGD",
        "h_"  + base + "_" + state
    };
    for (const auto& n : noTag) if (TH2* h = fetchClone2D(f, n)) return h;

    return nullptr;
}

// -------------------- utils --------------------
static void normalize(TH1* h) {
    if (!h) return;
    const double I = h->Integral();
    if (I > 0) h->Scale(1.0 / I);
}

static double maxOf(const std::vector<TH1*>& hs) {
    double m = 0.0;
    for (auto* h : hs) if (h && h->GetMaximum() > m) m = h->GetMaximum();
    return m;
}

// -------------------- drawers --------------------
// 1D overlay (pre/post). Titles do NOT include the state.
static void drawOne(const Var& v,
                    const std::string& state,
                    const std::map<std::string, std::unique_ptr<TFile>>& files,
                    const std::string& pdfName)
{
    const std::string& base   = v.base;
    const std::string& xtitle = v.xtitle;

    std::vector<TH1*> hs;
    hs.reserve(kTargets.size());

    const std::string title = v.title; // no pre/post in title

    for (const auto& t : kTargets) {
        TH1* h = getHistFlexible(files.at(t.tag).get(), base, state, t.tag);
        if (!h) {
            std::cerr << "[WARN] Missing " << base << "_" << state
                      << " in " << t.path << " (tried common variants)\n";
        }
        if (h) {
            h->SetLineColor(t.color);
            h->SetLineWidth(2);
            h->SetTitle(title.c_str());
            h->GetXaxis()->SetTitle(xtitle.c_str());
            h->GetYaxis()->SetTitle("Normalized events");
            normalize(h);
            hs.push_back(h);
        }
    }

    if (hs.size() < 2) {
        std::cerr << "[WARN] Not enough histograms to draw " << base << "_" << state << "\n";
        return;
    }

    const double ymax = maxOf(hs) * 1.15;
    std::unique_ptr<TCanvas> c(new TCanvas(
        ("c_" + base + "_" + state).c_str(), "", 1000, 700));

    hs.front()->SetMaximum(ymax);
    hs.front()->SetMinimum(0);
    hs.front()->Draw("hist");
    for (size_t i = 1; i < hs.size(); ++i) hs[i]->Draw("hist same");

    // legend
    std::unique_ptr<TLegend> leg(new TLegend(0.70, 0.70, 0.93, 0.90));
    leg->SetTextSize(0.045);
    for (const auto& t : kTargets) {
        for (auto* h : hs) {
            if (h->GetLineColor() == t.color) { leg->AddEntry(h, t.label.c_str(), "l"); break; }
        }
    }
    leg->Draw();

    // dashed cut lines
    auto [xlo, xhi] = cutsFor(base);
    if (xlo != 100.0) {
        auto* l = new TLine(xlo, 0.0, xlo, ymax);
        l->SetLineStyle(2);   // dashed
        l->SetLineWidth(3);
        l->Draw("same");
    }
    if (xhi != 100.0) {
        auto* r = new TLine(xhi, 0.0, xhi, ymax);
        r->SetLineStyle(2);   // dashed
        r->SetLineWidth(3);
        r->Draw("same");
    }
    gPad->Modified();
    gPad->Update();

    // save
    c->SaveAs(Form("finComp_%s_%s.png", base.c_str(), state.c_str()));
    c->Print(pdfName.c_str());
}

// 2D panel (2x2, one pad per target), shared Z-scale, titles without pre/post
static void drawOne2DPanel(const std::string& base,
                           const std::string& prettyTitle,
                           const std::string& state,
                           const std::map<std::string, std::unique_ptr<TFile>>& files,
                           const std::string& pdfName)
{
    struct Item { const Target* tgt; TH2* h; };
    std::vector<Item> items; items.reserve(kTargets.size());
    double zmax = 0.0;

    for (const auto& t : kTargets) {
        TH2* h2 = getHist2DFlexible(files.at(t.tag).get(), base, state, t.tag);
        if (!h2) {
            std::cerr << "[WARN] Missing TH2 " << base << "_" << state
                      << " in " << t.path << " (tried with/without tag)\n";
            continue;
        }
        zmax = std::max(zmax, h2->GetMaximum());
        items.push_back({&t, h2});
    }
    if (items.empty()) return;

    std::unique_ptr<TCanvas> c(new TCanvas(
        ("c2_" + base + "_" + state).c_str(), "", 1200, 900));
    c->Divide(2, 2, 0.001, 0.001);

    int pad = 1;
    for (const auto& t : kTargets) {
        c->cd(pad++);
        gPad->SetRightMargin(0.12);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);

        TH2* h2 = nullptr;
        for (const auto& it : items) if (it.tgt->tag == t.tag) { h2 = it.h; break; }
        if (!h2) { gPad->Clear(); continue; }

        // Title: no "pre/post", keep target name; ASCII hyphen
        h2->SetTitle(Form("%s - %s", prettyTitle.c_str(), t.label.c_str()));
        h2->SetContour(60);
        h2->SetMaximum(zmax);
        h2->Draw("COLZ");
    }

    c->SaveAs(Form("finComp2D_%s_%s.png", base.c_str(), state.c_str()));
    c->Print(pdfName.c_str());
}

// 2D single-plot per target (full-size), titles without pre/post
static void drawOne2DSoloPerTarget(const std::string& base,
                                   const std::string& prettyTitle,
                                   const std::string& state,
                                   const std::map<std::string, std::unique_ptr<TFile>>& files,
                                   const std::string& pdfName)
{
    for (const auto& t : kTargets) {
        TH2* h2 = getHist2DFlexible(files.at(t.tag).get(), base, state, t.tag);
        if (!h2) continue;

        std::unique_ptr<TCanvas> c(new TCanvas(
            ("c2solo_" + base + "_" + state + "_" + t.tag).c_str(), "", 1100, 800));
        gPad->SetRightMargin(0.12);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);

        h2->SetTitle(Form("%s - %s", prettyTitle.c_str(), t.label.c_str()));
        h2->SetContour(60);
        h2->Draw("COLZ");

        c->SaveAs(Form("finComp2D_%s_%s_%s.png", base.c_str(), state.c_str(), t.tag.c_str()));
        c->Print(pdfName.c_str());
    }
}

// -------------------- main --------------------
int main() {
    setStyle();
    // gROOT->SetBatch(true); // uncomment to suppress GUI windows if desired

    // open files
    std::map<std::string, std::unique_ptr<TFile>> files;
    for (const auto& t : kTargets) {
        std::unique_ptr<TFile> f(TFile::Open(t.path.c_str(), "READ"));
        if (!f || f->IsZombie() || !f->IsOpen()) {
            std::cerr << "[ERROR] Cannot open " << t.path << "\n";
            return 1;
        }
        files[t.tag] = std::move(f);
    }

    // list contents to stdout
    for (const auto& t : kTargets)
        DumpHistListToStdout(files[t.tag].get(), t.tag, t.path);

    const std::string pdfName = "finComp.pdf";
    { TCanvas opener("opener","",10,10); opener.Print((pdfName + "[").c_str()); }

    // 1D overlays (pre & post). Titles do not include state.
    for (const auto& v : kVars) {
        drawOne(v, "pre",  files, pdfName);
        drawOne(v, "post", files, pdfName);
    }

    // 2D: panel + solo-per-target (titles without pre/post)
    for (const auto& v2 : kVars2D) {
        for (const auto& st : kStates2D) {
            drawOne2DPanel(v2.base, v2.title, st, files, pdfName);
            drawOne2DSoloPerTarget(v2.base, v2.title, st, files, pdfName);
        }
    }

    { TCanvas closer("closer","",10,10); closer.Print((pdfName + "]").c_str()); }
    std::cout << "\nWrote finComp.pdf, finComp_<var>_<state>.png,"
                 " finComp2D_<base>_<state>.png, and finComp2D_<base>_<state>_<TAG>.png\n";
    return 0;
}
