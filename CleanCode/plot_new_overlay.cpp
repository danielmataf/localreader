#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLine.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "constants.h"

// Compile & run:
// g++ plot_new_overlay.cpp -o plotbinsxQ `root-config --cflags --glibs`
// ./plotbinsxQ

constexpr double M = 0.938272;
constexpr double E = 10.4;

constexpr double THETA_HI  = 27.0;
constexpr double THETA_MID = 11;
constexpr double THETA_LOW = 8.8;

inline double deg2rad(double d) { return d * M_PI / 180.0; }

double Q2_of_x_theta(double x, double thetaDeg) {
    double c = 1.0 - std::cos(deg2rad(thetaDeg));
    double A = 2.0 * M * E * E * c;
    double B = E * c;
    return A * x / (B + M * x);
}
double Q2_of_x_y(double x, double y) {
    double s = 2.0 * M * E;
    return x * y * s;
}

double Q2_of_x_W2(double x, double W2) {
    return (W2 - M * M) * x / (1.0 - x);
}

int main() {
    //TFile* fin = TFile::Open("DptinputFiles/REcutC2_test.root", "READ");
    TFile* fin = TFile::Open("~/pass0v11C2.root", "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open input ROOT file.\n";
        return 1;
    }

    auto* h = dynamic_cast<TH2F*>(fin->Get("xQ2_C2_RGD"));
    if (!h) {
        std::cerr << "ERROR: histogram xQ2_C2_RGD not found.\n";
        return 2;
    }

    // ===== Five vertical xB lines (5 edges → 4 closed ranges + tail) =====
    std::vector<double> xBins = {0.075, 0.11, 0.15, 0.19, 0.29}; // << exactly 5 lines
    if (xBins.size() != 5) {
        std::cerr << "ERROR: need exactly 5 x-bin edges (5 vertical lines); got "
                  << xBins.size() << "\n";
        return 3;   
    }
        const double regABx_LO  = xBins.at(0); // 0.06
        const double regABx_HI  = xBins.at(1); // 0.13
        const double regCDEx_LO = xBins.at(1); // 0.13
        const double regCDEx_HI = xBins.at(2); // 0.17
        const double regFGHx_LO = xBins.at(2); // 0.17
        const double regFGHx_HI = xBins.at(3); // 0.26
        const double regIJKx_LO = xBins.at(3); // 0.26
        const double regIJKx_HI = xBins.at(4); // 0.32
        const double regLMNx_LO = xBins.at(4); // 0.32
    // R1 ≈ old AB+CDE: [0.06, 0.13]
    // R2 ≈ old FGH:    [0.13, 0.17]
    // R3 ≈ old IJK+LMN:[0.17, 0.26]
    // R4 ≈ old OPQ:    [0.26, 0.36]
    // R5 ≈ old RS:     [0.36, +∞)
    const double R1_LO = xBins.at(0), R1_HI = xBins.at(1);
    const double R2_LO = xBins.at(1), R2_HI = xBins.at(2);
    const double R3_LO = xBins.at(2), R3_HI = xBins.at(3);
    const double R4_LO = xBins.at(3), R4_HI = xBins.at(4);
    const double R5_LO = xBins.at(4); // tail

    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c", "xQ2 with overlays", 1400, 800);
    c->SetRightMargin(0.15);

    // θ curves
    auto make_theta_curve = [](double thetaDeg, const char* name, int color) {
        double c = 1.0 - std::cos(deg2rad(thetaDeg));
        double A = 2.0 * M * E * E * c;
        double B = E * c;
        TF1* f = new TF1(name, "[0]*x/([1]+0.9396*x)", 0.0, 0.7);
        f->SetParameters(A, B);
        f->SetLineWidth(2);
        f->SetLineColor(color);
        f->SetLineStyle(1);
        f->SetRange(0.0, 0.7);
        f->SetMinimum(1.0);
        return f;
    };

    TF1* fThetaHi  = make_theta_curve(THETA_HI , "fThetaHi" , kRed);
    TF1* fThetaMid = make_theta_curve(THETA_MID, "fThetaMid", kRed + 1);
    TF1* fThetaLow = make_theta_curve(THETA_LOW, "fThetaLow", kRed + 2);
// y = 0.75 line    
    double y_fixed = 0.67;  //change y value here 
    double s = 2.0 * M * E;
    double C_y = y_fixed * s;

        // y = 0.25 line    (this is from low cut=0.25)
    double y_fixedlow = 0.25;  //change y value here 
    //double s = 2.0 * M * E;
    double C_y2 = y_fixedlow * s;


    TF1* fy = new TF1("fy", Form("%.5f * x", C_y), 0.0, 0.7);
    fy->SetLineColor(kMagenta + 1);
    fy->SetLineWidth(2);
    fy->SetLineStyle(2);  // dashed
    
    TF1* fylow = new TF1("fy2", Form("%.5f * x", C_y2), 0.0, 0.7);
    fylow->SetLineColor(kMagenta + 1);
    fylow->SetLineWidth(2);
    fylow->SetLineStyle(2);  // dashed

    std::vector<double> yVals = {0.75, 0.25};  // cuts in y 

    // W² = 4 GeV² line
    TF1* fW2 = new TF1("fW2", "(4 - 0.938272*0.938272)*x/(1 - x)", 0.0, 0.7);
    fW2->SetLineColor(kOrange + 1);
    fW2->SetLineWidth(2);
    fW2->SetLineStyle(1);
    fW2->SetMinimum(1.0);

    // Vertical xB lines (now 5)
    std::vector<TLine*> xLines;
    for (double x : xBins) {
        TLine* lx = new TLine(x, 1.0, x, 5.0);
        lx->SetLineColor(kBlue + 2);
        lx->SetLineWidth(2);
        lx->SetLineStyle(1);
        xLines.push_back(lx);
    }

    // Horizontal Q² lines at 1 and 5
    std::vector<double> Q2LinesToDraw = {1.0, 5.0};
    std::vector<TLine*> Q2Lines;
    double x_min = 0.05, x_max = 0.6;
    for (double Q2val : Q2LinesToDraw) {
        TLine* lQ = new TLine(x_min, Q2val, x_max, Q2val);
        lQ->SetLineColor(kGreen + 2);
        lQ->SetLineWidth(2);
        lQ->SetLineStyle(1);
        Q2Lines.push_back(lQ);
    }

    // Draw
    h->GetXaxis()->SetTitle("x_{B}");
    h->GetYaxis()->SetTitle("Q^{2}  [GeV^{2}]");
    h->Draw("colz");

    // Fix vertical line top to match the pad range
    double yTop = gPad->GetUymax();
    //for (auto* lx : xLines) { lx->SetY2(yTop); lx->Draw("same"); }

    //fThetaHi ->Draw("same");
    fThetaMid->Draw("same");
    fThetaLow->Draw("same");
    fy       ->Draw("same");
    fylow    ->Draw("same");
    fW2      ->Draw("same");
    for (auto* l : Q2Lines) l->Draw("same");
    for (auto* l : xLines)  l->Draw("same");


    TLegend* leg = new TLegend(0.70, 0.60, 0.93, 0.87);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    //leg->AddEntry(fThetaHi , Form("#theta_{e}=%.1f^{#circ}", THETA_HI ), "l");
    leg->AddEntry(fThetaMid, Form("#theta_{e}=%.1f^{#circ}", THETA_MID), "l");
    leg->AddEntry(fThetaLow, Form("#theta_{e}=%.1f^{#circ}", THETA_LOW), "l");
    leg->AddEntry(fW2, "W^{2} = 4 GeV^{2}", "l");
    leg->AddEntry(fy, "y = 0.75", "l");
    leg->AddEntry(fylow, "y = 0.25", "l");

    leg->Draw();

    c->SaveAs("xQ2_overlay_bis.png");

    // ================== Console output (reduced to R1..R5) ==================
    auto q2th = [](double x, double th){ return Q2_of_x_theta(x, th); };
    auto q2w  = [](double x){ return Q2_of_x_W2(x, 4.0); };

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nθ_HI=" << THETA_HI << " ; θ_MID=" << THETA_MID
              << " ; θ_LOW=" << THETA_LOW << " ;\n";

    std::cout << std::setprecision(2)
              << "R1: [" << R1_LO << ", " << R1_HI << "]\n"
              << "R2: [" << R2_LO << ", " << R2_HI << "]\n"
              << "R3: [" << R3_LO << ", " << R3_HI << "]\n"
              << "R4: [" << R4_LO << ", " << R4_HI << "]\n"
              << "R5: [" << R5_LO << ", ...]\n";

    std::cout << std::setprecision(3);

    auto dump_range = [&](const char* name, double xlo, double xhi, bool useWmin){
        const double q2_lo_mid = q2th(xlo, THETA_MID);
        const double q2_hi_hi  = q2th(xhi, THETA_HI);
        const double q2_lo_low = q2th(xlo, THETA_LOW);
        const double q2_hi_mid = q2th(xhi, THETA_MID);
        const double q2_w_lo   = std::max(1.0, q2w(xlo));

        std::cout << name << ": "
                  << "Q2(" << xlo << ", θ_mid)=" << q2_lo_mid << " ; "
                  << "Q2(" << xhi << ", θ_hi)="  << q2_hi_hi  << "\n"
                  << "    Q2(" << xlo << ", θ_low)=" << q2_lo_low << " ; "
                  << "Q2(" << xhi << ", θ_mid)="    << q2_hi_mid << "\n";
        if (useWmin) {
            std::cout << "    min Q2 via W2 at x_lo (clipped at 1): " << q2_w_lo << "\n";
        } else {
            std::cout << "    min Q2 fixed at 1.000\n";
        }
    };

    std::cout << "\n===== Region [FEW] Definitions  =====\n";
    std::cout << "A LO = Q2(" << regABx_LO << ", " << THETA_MID << ") = " << Q2_of_x_theta(regABx_LO,THETA_MID) << "  // loQregA\n";
    std::cout << "A HI = Q2(" << regABx_HI << ", y=" << y_fixed << ") = "<< Q2_of_x_y(regABx_HI, y_fixed) << "  // hiQregA\n";
    std::cout << std::endl;
    std::cout << "B LO = Q2(fixed) = 1     // loQregB \n";
    std::cout << "B HI = Q2(" << regABx_HI << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regABx_HI, THETA_MID) << "  // hiQregB\n";
    std::cout << std::endl;
    std::cout << "C LO = Q2(" << regCDEx_LO << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regCDEx_LO,THETA_MID  ) << "  // loQregC\n";
    std::cout << "C HI = Q2(" << regCDEx_HI << ", y=" << y_fixed << ") = "<< Q2_of_x_y(regCDEx_HI,y_fixed  ) << "  // hiQregC\n";
    std::cout << std::endl;
    std::cout << "D LO = Q2(" << regCDEx_LO << ", " << THETA_LOW << ") = "<< Q2_of_x_theta(regCDEx_LO,THETA_LOW  ) << "  // loQregD\n";
    std::cout << "D HI = Q2(" << regCDEx_HI << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regCDEx_HI,THETA_MID  ) << "  // hiQregD\n";
    std::cout << std::endl;
    std::cout << "E LO = Q2(fixed) = 1     //loQregE \n";
    std::cout << "E HI = Q2(" << regCDEx_HI << ", " << THETA_LOW << ") = "<< Q2_of_x_theta(regCDEx_HI,THETA_LOW  ) << "  // hiQregE\n";
    std::cout << std::endl;
    std::cout << "F LO = Q2(" << regFGHx_LO << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regFGHx_LO,THETA_MID  ) << "  // loQregF\n";
    std::cout << "F HI = Q2(" << regFGHx_HI << ", y=" << y_fixed << ") = "<< Q2_of_x_y(regFGHx_HI,y_fixed  ) << "  // hiQregF\n";
    std::cout << std::endl;
    std::cout << "G LO = Q2(" << regFGHx_LO << ", " << THETA_LOW << ") = "<< Q2_of_x_theta(regFGHx_LO,THETA_LOW  ) << "  // loQregG\n";
    std::cout << "G HI = Q2(" << regFGHx_HI << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regFGHx_HI,THETA_MID  ) << "  // hiQregG\n";
    std::cout << std::endl;
    std::cout << "H LO = Q2(" << regFGHx_LO << ", y=" << y_fixedlow << ") = "<< Q2_of_x_theta(regFGHx_LO,y_fixedlow  ) << " //loQregH \n";
    std::cout << "H HI = Q2(" << regFGHx_HI << ", " << THETA_LOW << ") = "<< Q2_of_x_theta(regFGHx_HI,THETA_LOW  ) << "  // hiQregH\n";
    std::cout << std::endl;
    std::cout << "I LO = Q2(" << regIJKx_LO << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regIJKx_LO,THETA_MID  ) << "  // loQregI\n";
    std::cout << "I HI = Q2(" << regIJKx_HI << ", y=" << y_fixed << ") = "<< Q2_of_x_y(regIJKx_HI,  y_fixed) << "  // hiQregI\n";
    std::cout << std::endl;
    std::cout << "J LO = Q2(" << regIJKx_LO << ", " << THETA_LOW << ") = "<< Q2_of_x_theta(regIJKx_LO,THETA_LOW  ) << "  // loQregJ\n";
    std::cout << "J HI = Q2(" << regIJKx_HI << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regIJKx_HI,THETA_MID  ) << "  // hiQregJ\n";
    std::cout << std::endl;
    std::cout << "K LO = Q2(" << regIJKx_LO << ", y=" << y_fixedlow << ") = "<< Q2_of_x_y(regIJKx_HI,y_fixedlow  ) << "//loQregK \n";
    std::cout << "K HI = Q2(" << regIJKx_HI << ", " << THETA_LOW << ") = "<< Q2_of_x_theta(regIJKx_HI,THETA_LOW  ) << "  // hiQregK\n";
    std::cout << std::endl;

    std::cout << "L LO = Q2(" << regLMNx_LO << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regLMNx_LO,THETA_MID  ) << "  // loQregL\n";
    std::cout << "non Q2 limit fixed for L-M region    TBD  //loQregRS \n";


    fin->Close();
    return 0;
}


