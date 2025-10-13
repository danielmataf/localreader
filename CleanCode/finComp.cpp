// finComp.cpp
// g++ finComp.cpp $(root-config --cflags --libs) -o finComp && ./finComp
// this compares pre & pos cuts
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

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <array>
#include <algorithm>

// -------- Cut values (100.0 sentinel = "don't draw") --------
double cut_loQ2       = 1.0;   double cut_hiQ2       = 100.0;
double cut_loW2       = 4.0;   double cut_hiW2       = 100.0;
double cut_lonu       = 100.0; double cut_hinu       = 7.0;
double cut_lotargetVz = 100.0; double cut_hitargetVz = 100.0;
double cut_lopt2      = 100.0; double cut_hipt2      = 1.2;

// -------- Targets (order defines legend order & colors) --------
struct Target {
    std::string tag;     // Cu, CxC, LD2, Sn
    std::string label;   // legend label
    int         color;   // ROOT color
    std::string path;    // file path
};

static const std::vector<Target> kTargets = {
    {"Cu",  "Cu",        kGreen+2,  "build/fin_Cu_test.root"},
    {"CxC", "C (CxC)",   kBlack,    "build/fin_CxC_test.root"},
    {"LD2", "LD2",       kBlue,     "build/fin_LD2_test.root"},
    {"Sn",  "Sn",        kOrange+7, "build/fin_Sn_test.root"}
};

// -------- Variables to draw (exact set you requested) --------
struct Var { std::string base; std::string xtitle; };
static const std::vector<Var> kVars = {
    {"Q2",  "Q^{2} [GeV^{2}]"},
    {"Vz",  "Target V_{z} [cm]"},
    {"W2",  "W^{2} [GeV^{2}]"},
    {"nu",  "#nu [GeV]"},
    {"pt2", "p_{T}^{2} [GeV^{2}]"}
};

static std::pair<double,double> cutsFor(const std::string& base) {
    if (base == "Q2")  return {cut_loQ2,       cut_hiQ2};
    if (base == "W2")  return {cut_loW2,       cut_hiW2};
    if (base == "nu")  return {cut_lonu,       cut_hinu};
    if (base == "Vz")  return {cut_lotargetVz, cut_hitargetVz};
    if (base == "pt2") return {cut_lopt2,      cut_hipt2};
    return {100.0, 100.0};
}

static void setStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
}

// ------------------ HISTO LISTING (stdout only) ------------------
static void ListHistosRecursive(TDirectory* dir,
                                std::vector<std::pair<std::string,std::string>>& out,
                                const std::string& prefix = "")
{
    TIter next(dir->GetListOfKeys());
    while (TKey* key = (TKey*)next()) {
        const char* cname = key->GetClassName();
        const char* kname = key->GetName();
        TClass*     cls   = gROOT->GetClass(cname);
        if (!cls) continue;

        if (cls->InheritsFrom(TDirectory::Class())) {
            if (TDirectory* sub = dir->GetDirectory(kname)) {
                ListHistosRecursive(sub, out, prefix + kname + "/");
            }
        } else if (cls->InheritsFrom(TH1::Class())) {
            // mark dimensionality for clarity
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
    for (auto& kv : names) {
        std::cout << "  [" << kv.first << "] " << kv.second << "\n";
    }
    if (names.empty()) {
        std::cout << "  (no histograms found)\n";
    }
}
// ----------------------------------------------------------------

// Fetch and clone a TH1 if present
static TH1* fetchClone(TFile* f, const std::string& name) {
    TObject* obj = f->Get(name.c_str());
    if (!obj || !obj->InheritsFrom(TH1::Class())) return nullptr;
    TH1* h = static_cast<TH1*>(obj);
    TH1* c = static_cast<TH1*>(h->Clone((name + "_cl").c_str()));
    c->SetDirectory(nullptr);
    return c;
}

// Try multiple naming patterns and synonyms (handles per-file vs per-hist target tags)
static TH1* getHistFlexible(TFile* f,
                            const std::string& base,
                            const std::string& state,
                            const std::string& tag)
{
    // Optional synonyms for base
    std::vector<std::string> bases = { base };
    if (base == "Vz")  bases.push_back("targetVz");
    if (base == "pt2") bases.push_back("pT2");

    // Name patterns to try, in order:
    //  base_state_tag_RGD, base_state_tag, base_state, h_base_state
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

// Draw one overlay for given variable & state ("pre" or "post")
static void drawOne(const std::string& base,
                    const std::string& xtitle,
                    const std::string& state,
                    const std::map<std::string, std::unique_ptr<TFile>>& files,
                    const std::string& pdfName)
{
    // Collect the four hists (one per target)
    std::vector<TH1*> hs;
    hs.reserve(kTargets.size());

    const std::string title = base + " (" + state + " " + base + " cut)";

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
    hs.front()->GetXaxis()->SetTitle(xtitle.c_str());
    hs.front()->GetYaxis()->SetTitle("Normalized events");
    hs.front()->Draw("hist");
    for (size_t i = 1; i < hs.size(); ++i) hs[i]->Draw("hist same");

    // Legend in fixed order (Cu, CxC, LD2, Sn)
    std::unique_ptr<TLegend> leg(new TLegend(0.70, 0.70, 0.93, 0.90));
    leg->SetTextSize(0.045);
    for (const auto& t : kTargets) {
        for (auto* h : hs) {
            if (h->GetLineColor() == t.color) { leg->AddEntry(h, t.label.c_str(), "l"); break; }
        }
    }
    leg->Draw();

    // Optional cut lines
    auto [xlo, xhi] = cutsFor(base);
    if (xlo != 100.0) { TLine l(xlo, 0.0, xlo, ymax); l.SetLineStyle(3); l.SetLineWidth(2); l.Draw(); }
    if (xhi != 100.0) { TLine r(xhi, 0.0, xhi, ymax); r.SetLineStyle(2); r.SetLineWidth(2); r.Draw(); }

    // Save PNG and append to PDF
    c->SaveAs(Form("finComp_%s_%s.png", base.c_str(), state.c_str()));
    c->Print(pdfName.c_str());
}

int main() {
    setStyle();
    // gROOT->SetBatch(true); // uncomment to suppress GUI windows if desired

    // Open files
    std::map<std::string, std::unique_ptr<TFile>> files; // tag -> file
    for (const auto& t : kTargets) {
        std::unique_ptr<TFile> f(TFile::Open(t.path.c_str(), "READ"));
        if (!f || f->IsZombie() || !f->IsOpen()) {
            std::cerr << "[ERROR] Cannot open " << t.path << "\n";
            return 1;
        }
        files[t.tag] = std::move(f);
    }

    // ---- LIST HISTOGRAMS TO STDOUT (no files written) ----
    for (const auto& t : kTargets) {
        DumpHistListToStdout(files[t.tag].get(), t.tag, t.path);
    }

    // Open combined PDF
    const std::string pdfName = "finComp.pdf";
    { TCanvas opener("opener","",10,10); opener.Print((pdfName + "[").c_str()); }

    // Draw pre & post for each requested variable
    for (const auto& v : kVars) {
        drawOne(v.base, v.xtitle, "pre",  files, pdfName);
        drawOne(v.base, v.xtitle, "post", files, pdfName);
    }

    // Close PDF
    { TCanvas closer("closer","",10,10); closer.Print((pdfName + "]").c_str()); }

    std::cout << "\nWrote finComp.pdf and finComp_<var>_<state>.png\n";
    return 0;
}
