#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <map>

//compile and run 
//g++ -o comphist compHistoJun34.cpp root-config --cflags --libs
//g++ -o comphist compHistoJun34.cpp $(root-config --cflags --libs)
//./comphist

void CompareHistograms(const char* target, const std::vector<std::string>& plotTitles, const std::vector<std::string>& xTitles) {
    //std::string file1 = std::string("/home/matamoros/sept") + target + "_test.root";
    //std::string file2 = std::string("/home/matamoros/sept") + target + "_sim.root";
    std::string file1 = std::string("/home/matamoros/elev") + target + "_test.root";
    std::string file2 = std::string("/home/matamoros/elev") + target + "_sim.root";

    TFile* rootFile1 = new TFile(file1.c_str(), "READ");
    if (!rootFile1->IsOpen()) {
        std::cerr << "Error: Cannot open file1!" << std::endl;
        return;
    }
    TFile* rootFile2 = new TFile(file2.c_str(), "READ");
    if (!rootFile2->IsOpen()) {
        std::cerr << "Error: Cannot open file2!" << std::endl;
        return;
    }

    const std::vector<std::pair<std::string, std::string>> histogramPairs = {
        {"Q2_" + std::string(target) + "_true", "Q2_" + std::string(target) + "_sim"},
        {"W2_" + std::string(target) + "_true", "W2_" + std::string(target) + "_sim"},
        {"nu_" + std::string(target) + "_true", "nu_" + std::string(target) + "_sim"},
        {"phih_" + std::string(target) + "_true", "phih_" + std::string(target) + "_sim"},
        {"xb_" + std::string(target) + "_true", "xb_" + std::string(target) + "_sim"},
        {"y_" + std::string(target) + "_true", "y_" + std::string(target) + "_sim"},
        {"z_" + std::string(target) + "_true", "z_" + std::string(target) + "_sim"},
        {"targetVz_" + std::string(target) + "_true", "targetVz_" + std::string(target) + "_sim"},
        {"pt2_" + std::string(target) + "_true", "pt2_" + std::string(target) + "_sim"},
        {"ptot_ele_" + std::string(target) + "_true", "ptot_ele_" + std::string(target) + "_sim"},
        {"px_ele_" + std::string(target) + "_true", "px_ele_" + std::string(target) + "_sim"},
        {"py_ele_" + std::string(target) + "_true", "py_ele_" + std::string(target) + "_sim"},
        {"pz_ele_" + std::string(target) + "_true", "pz_ele_" + std::string(target) + "_sim"},
        {"E_el" + std::string(target) + "_true", "E_el" + std::string(target) + "_sim"},
        {"E_pi" + std::string(target) + "_true", "E_pi" + std::string(target) + "_sim"},
        {"theta_el" + std::string(target) + "_true", "theta_el" + std::string(target) + "_sim"},
        {"phi_el" + std::string(target) + "_true", "phi_el" + std::string(target) + "_sim"},
        {"ptot_pro_" + std::string(target) + "_true", "ptot_pro_" + std::string(target) + "_sim"},
        {"px_pro_" + std::string(target) + "_true", "px_pro_" + std::string(target) + "_sim"},
        {"py_pro_" + std::string(target) + "_true", "py_pro_" + std::string(target) + "_sim"},
        {"pz_pro_" + std::string(target) + "_true", "pz_pro_" + std::string(target) + "_sim"},
        {"theta_pi" + std::string(target) + "_true", "theta_pi" + std::string(target) + "_sim"},
        {"phi_pi" + std::string(target) + "_true", "phi_pi" + std::string(target) + "_sim"},
        {"lu_el"+ std::string(target) + "_true", "lu_el" + std::string(target) + "_sim"},
        {"lv_el"+ std::string(target) + "_true", "lv_el" + std::string(target) + "_sim"},
        {"lw_el"+ std::string(target) + "_true", "lw_el" + std::string(target) + "_sim"},
        {"epcal_el"+ std::string(target) + "_true", "epcal_el" + std::string(target) + "_sim"},
        {"Nphe15_"+ std::string(target) + "_true", "Nphe15_" + std::string(target) + "_sim"},
        {"Nphe16_"+ std::string(target) + "_true", "Nphe16_" + std::string(target) + "_sim"}
    };

    std::map<std::string, int> targetMap = {
        {"LD2", 1},
        {"Sn", 2},
        {"C1", 3},
        {"Cu", 4},
        {"C2", 5}
    };

    double lowvz, highvz;

    switch (targetMap[target]) {
        case 1: // LD2
            lowvz = -7.5;
            highvz = -2.5;
            break;
        case 2: // Sn
        case 3: // C1
            lowvz = -3.5;
            highvz = -1.5;
            break;
        case 4: // Cu
        case 5: // C2
            lowvz = -8.5;
            highvz = -6.5;
            break;
        default:
            std::cerr << "Error: Unknown target!" << std::endl;
            return;
    }

    const std::map<std::string, std::vector<double>> xStartValues = {
        {"Q2_" + std::string(target) + "_true", {1.5, 100.0}},
        {"W2_" + std::string(target) + "_true", {0.0, 100.0}},
        {"nu_" + std::string(target) + "_true", {100.0, 100.0}}, // No lines for "nu"
        {"phih_" + std::string(target) + "_true", {100.0, 100.0}}, // No lines for "phih"
        {"xb_" + std::string(target) + "_true", {100.0, 100.0}}, // No lines for "xb"
        {"y_" + std::string(target) + "_true", {100.0, 100.0}}, // No lines for "y2"
        {"z_" + std::string(target) + "_true", {0.3, 0.7}},
        {"targetVz_" + std::string(target) + "_true", {lowvz, highvz}},
        {"pt2_" + std::string(target) + "_true", {-10.0, 1.5}},
        {"ptot_ele_" + std::string(target) + "_true", {100.0, 100.0}},  // exaggerating values because no lines here
        {"px_ele_" + std::string(target) + "_true", {100.0, 100.0}}, 
        {"py_ele_" + std::string(target) + "_true", {100.0, 100.0}}, 
        {"pz_ele_" + std::string(target) + "_true", {100.0, 100.0}},
        {"E_el" + std::string(target) + "_true", {100.0, 100.0}},
        {"E_pi" + std::string(target) + "_true", {100.0, 100.0}},
        {"theta_el" + std::string(target) + "_true", {100.0, 100.0}}, 
        {"phi_el" + std::string(target) + "_true", {100.0, 100.0}}, 
        {"ptot_pro_" + std::string(target) + "_true", {100.0, 100.0}}, 
        {"px_pro_" + std::string(target) + "_true", {100.0, 100.0}}, 
        {"py_pro_" + std::string(target) + "_true", {100.0, 100.0}}, 
        {"pz_pro_" + std::string(target) + "_true", {100.0, 100.0}},
        {"theta_pi" + std::string(target) + "_true", {100.0, 100.0}}, 
        {"phi_pi" + std::string(target) + "_true", {100.0, 100.0}}, 
        {"lu_el"+ std::string(target) + "_true", {100.0, 100.0}}, 
        {"lv_el"+ std::string(target) + "_true", {100.0, 100.0}}, 
        {"lw_el"+ std::string(target) + "_true", {100.0, 100.0}},
        {"epcal_el"+ std::string(target) + "_true", {100.0, 100.0}}, 
        {"Nphe15_"+ std::string(target) + "_true", {100.0, 100.0}}, 
        {"Nphe16_"+ std::string(target) + "_true", {100.0, 100.0}} 
    };

    //map defining x-axis titles and plot titles & association to histos  
    std::map<std::string, std::pair<std::string, std::string>> titleMap = {
        {"pt2_" + std::string(target) + "_true", {"Comparison of p_{t}^{2}", "p_{t}^{2} (GeV^{2})"}},
        {"Q2_" + std::string(target) + "_true", {"Comparison of Q^{2}", "Q^{2} (GeV^{2})"}},
        {"W2_" + std::string(target) + "_true", {"Comparison of W^{2}", "W^{2} (GeV^{2})"}},
        {"xb_" + std::string(target) + "_true", {"Comparison of x_{B}", "x_{B}"}},
        {"y_" + std::string(target) + "_true", {"Comparison of y", "y"}},
        {"z_" + std::string(target) + "_true", {"Comparison of z", "z"}},
        {"targetVz_" + std::string(target) + "_true", {"Comparison of V_{z}", "V_{z} (cm)"}},
        {"phih_" + std::string(target) + "_true", {"Comparison of \\phi_{h}", "\\phi_{h} (deg)"}},
        {"ptot_ele_" + std::string(target) + "_true", {"Comparison of p_{tot, e}", "p_{tot, e} (GeV)"}},
        {"px_ele_" + std::string(target) + "_true", {"Comparison of p_{x, e}", "p_{x, e} (GeV)"}},
        {"py_ele_" + std::string(target) + "_true", {"Comparison of p_{y, e}", "p_{y, e} (GeV)"}},
        {"pz_ele_" + std::string(target) + "_true", {"Comparison of p_{z, e}", "p_{z, e} (GeV)"}},
        {"E_el" + std::string(target) + "_true", {"Comparison of E_{e}", "E_{e} (GeV)"}},
        {"E_pi" + std::string(target) + "_true", {"Comparison of E_{\\pi}", "E_{\\pi} (GeV)"}},
        {"theta_el" + std::string(target) + "_true", {"Comparison of \\theta_{e}", "\\theta_{e} (deg)"}},
        {"phi_el" + std::string(target) + "_true", {"Comparison of \\phi_{e}", "\\phi_{e} (deg)"}},
        {"ptot_pro_" + std::string(target) + "_true", {"Comparison of p_{tot, \\pi}", "p_{tot, \\pi} (GeV)"}},
        {"px_pro_" + std::string(target) + "_true", {"Comparison of p_{x, \\pi}", "p_{x, \\pi} (GeV)"}},
        {"py_pro_" + std::string(target) + "_true", {"Comparison of p_{y, \\pi}", "p_{y, \\pi} (GeV)"}},
        {"pz_pro_" + std::string(target) + "_true", {"Comparison of p_{z, \\pi}", "p_{z, \\pi} (GeV)"}},
        {"theta_pi" + std::string(target) + "_true", {"Comparison of \\theta_{\\pi}", "\\theta_{\\pi} (deg)"}},
        {"phi_pi" + std::string(target) + "_true", {"Comparison of \\phi_{\\pi}", "\\phi_{\\pi} (deg)"}},
        {"lu_el"+ std::string(target) + "_true", {"Comparison of l_{u}", "l_{u} (cm)"}},
        {"lv_el"+ std::string(target) + "_true", {"Comparison of l_{v}", "l_{v} (cm)"}},
        {"lw_el"+ std::string(target) + "_true", {"Comparison of l_{w}", "l_{w} (cm)"}},
        {"epcal_el"+ std::string(target) + "_true", {"Comparison of E_{pcal e}", "E_{pcal e} (GeV)"}},
        {"Nphe15_"+ std::string(target) + "_true", {"Comparison of Nphe15", "Nphe15"}},
        {"Nphe16_"+ std::string(target) + "_true", {"Comparison of Nphe16", "Nphe16"}}

    };

    for (size_t i = 0; i < histogramPairs.size(); ++i) {
        const auto& pair = histogramPairs[i];
        TCanvas* canvas = new TCanvas(Form("ComparisonCanvas_%s_%s", pair.first.c_str(), pair.second.c_str()), Form("Histogram Comparison %s vs %s", pair.first.c_str(), pair.second.c_str()), 800, 600);
        
        TH1F* h1 = dynamic_cast<TH1F*>(rootFile1->Get(pair.first.c_str()));     //data
        TH1F* h2 = dynamic_cast<TH1F*>(rootFile2->Get(pair.second.c_str()));    //sim

        if (!h1 || !h2) {
            std::cerr << "Error: Cannot retrieve histograms from files!" << std::endl;
            return;
        }
        if (pair.first == "pt2_" + std::string(target) + "_sim") { // || pair.first == "U_z_" + std::string(target) + "_data") {
            canvas->SetLogy();
            //maybe change this if statement for a switch case with only one case with pt2, then a case default where we do nothing (for the rest of the histos)
        }
        h1->SetLineColor(kBlue);    //data
        h2->SetLineColor(kRed);     //sim

        std::string plotTitle;
        std::string xTitle;

        
        auto it = titleMap.find(pair.first); // using find to get title and x-axis label based on the histogram name
        if (it != titleMap.end()) {
            plotTitle = it->second.first;
            xTitle = it->second.second;
        } else {
            // if there is no match for a title use a default title "comparison"
            plotTitle = "Comparison Sim vs RGD";
            xTitle = "";    // no x-axis label
        }

        h1->SetTitle(plotTitle.c_str());
        h1->GetXaxis()->SetTitle(xTitle.c_str());

        // Scaling
        double integral1 = h1->Integral();
        double integral2 = h2->Integral();
        if (integral1 != 0) h1->Scale(1.0 / integral1);
        if (integral2 != 0) h2->Scale(1.0 / integral2);

        // retrieving max post scaling 
        float max1 = h1->GetMaximum();
        float max2 = h2->GetMaximum();
        // getting the max of 2 histos. useful to define a proper yaxis range
        float max = (max1 > max2) ? max1 : max2;

        // yaxis range
        h1->SetMaximum(max * 1.1);  // multiplying max to leave some space in the axis
        h1->SetMinimum(0);          // forcing start at 0

        h1->Draw("hist");
        h2->Draw("hist same");

        TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.69);  // Adjusted position
        legend->AddEntry(h1, "Data", "l");
        legend->AddEntry(h2, "Simulation", "l");
        legend->Draw();

        // Drawing lines for specific histograms based on xStartValues map
        auto xStartIt = xStartValues.find(pair.first);
        if (xStartIt != xStartValues.end()) {
            for (double xStart : xStartIt->second) {
                TLine* line = new TLine(xStart, 0, xStart, max * 1.1);
                line->SetLineColor(kBlack);
                line->SetLineStyle(2); // dashed line
                line->Draw();
            }
        }

        canvas->SaveAs(Form("septComparison_%s_vs_%s.png", pair.first.c_str(), pair.second.c_str()));
    }

    rootFile1->Close();
    rootFile2->Close();
}

int main() {
    const char* target = "C2"; // Change this to the desired target
    CompareHistograms(target, {}, {});
    return 0;
}
