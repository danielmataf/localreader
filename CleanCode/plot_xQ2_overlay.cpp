#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLine.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <cmath>
#include <iomanip> // put this at the top of the file
#include <iostream>
#include <vector>

// Compile & run:
// g++ plot_xQ2_overlay.cpp -o plot_xQ2_overlay `root-config --cflags --glibs`
// ./plot_xQ2_overlay

constexpr double M = 0.9396;
constexpr double E = 10.5;

constexpr double THETA_HI  = 27.0;
constexpr double THETA_MID = 10.3;
constexpr double THETA_LOW = 8.4;

inline double deg2rad(double d) { return d * M_PI / 180.0; }

double Q2_of_x_theta(double x, double thetaDeg) {
    double c = 1.0 - std::cos(deg2rad(thetaDeg)); //this is sin sq (val/2) c= 2*sin2(thet/2)
    double A = 2.0 * M * E * E * c; // =2*M*E*E*c
    double B = E * c; // = 2*E*sin2(theta/2)
    //x is tha value of xb as input 
    return A * x / (B + M * x);
}

double Q2_of_x_y(double x, double y) {
    double s = 2.0 * M * E;
    return x * y * s;
}

double x_of_Q2_theta(double Q2, double thetaDeg) {
    double c = 1.0 - std::cos(deg2rad(thetaDeg));
    return (E * c * Q2) / (M * (2.0 * E * E * c - Q2));
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
    //auto* h = dynamic_cast<TH2F*>(fin->Get("xQ2_Sn_RGD"));
    if (!h) {
        std::cerr << "ERROR: histogram xQ2_C2_RGD not found.\n";
        return 2;
    }

    // Final xB lines (user defined + two chosen): total 4
    //std::vector<double> xBins = {0.075, 0.1, 0.13, 0.17, 0.21, 0.26, 0.34};  // 6 bin edges → 5 xB regions
    std::vector<double> xBins = {0.075, 0.105, 0.13, 0.16, 0.20, 0.25, 0.36};  // 6 bin edges → 5 xB regions
        std::vector<double> Q2LinesToDraw = {1.0, 5.0};        // keep only Q² = 1 and Q² = 5
        if (xBins.size() < 7) {
            std::cerr << "ERROR: need at least 7 x-bin edges; got " << xBins.size() << "\n";
            return 3;
        }

        const double regABx_LO  = xBins.at(0); // 0.06
        const double regABx_HI  = xBins.at(1); // 0.09
        const double regCDEx_LO = xBins.at(1); // 0.09
        const double regCDEx_HI = xBins.at(2); // 0.13
        const double regFGHx_LO = xBins.at(2); // 0.13
        const double regFGHx_HI = xBins.at(3); // 0.17
        const double regIJKx_LO = xBins.at(3); // 0.17
        const double regIJKx_HI = xBins.at(4); // 0.22
        const double regLMNx_LO = xBins.at(4); // 0.22
        const double regLMNx_HI = xBins.at(5); // 0.28
        const double regOPQx_LO = xBins.at(5); // 0.28
        const double regOPQx_HI = xBins.at(6); // 0.37
        const double regRSx_LO  = xBins.at(6); // 0.37
        // no regRSx_HI

    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c", "xQ2 with overlays", 1400, 800);
    c->SetRightMargin(0.15);

    // Create θ curves
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
    TF1* fW2 = new TF1("fW2", "(4 - 0.9396*0.9396)*x/(1 - x)", 0.0, 0.7);
    fW2->SetLineColor(kOrange + 1);
    fW2->SetLineWidth(2);
    fW2->SetLineStyle(1);
    fW2->SetMinimum(1.0);

    

    // Vertical lines: full height
    std::vector<TLine*> xLines;
    for (double x : xBins) {
        TLine* lx = new TLine(x, 1.0, x, 5.0);
        lx->SetLineColor(kBlue + 2);
        lx->SetLineWidth(2);
        lx->SetLineStyle(1);
        xLines.push_back(lx);
    }

    // Horizontal lines: only Q² = 1 and Q² = 5
    std::vector<TLine*> Q2Lines;
    double x_min = 0.05;
    double x_max = 0.6;
    for (double Q2val : Q2LinesToDraw) {
        TLine* lQ = new TLine(x_min, Q2val, x_max, Q2val);
        lQ->SetLineColor(kGreen + 2);
        lQ->SetLineWidth(2);
        lQ->SetLineStyle(1);
        Q2Lines.push_back(lQ);
    }

    // Draw everything
    h->GetXaxis()->SetTitle("x_{B}");
    h->GetYaxis()->SetTitle("Q^{2}  [GeV^{2}]");
    h->Draw("colz");

    //fThetaHi ->Draw("same");
    fThetaMid->Draw("same");
    fThetaLow->Draw("same");
    fW2      ->Draw("same");
    fy       ->Draw("same");
    fylow    ->Draw("same");


    for (auto* l : xLines)  l->Draw("same");
    for (auto* l : Q2Lines) l->Draw("same");

    TLegend* leg = new TLegend(0.70, 0.60, 0.93, 0.87);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    //leg->AddEntry(fThetaHi , Form("#theta_{e}=%.1f^{#circ}", THETA_HI ), "l");
    leg->AddEntry(fThetaMid, Form("#theta_{e}=%.1f^{#circ}", THETA_MID), "l");
    leg->AddEntry(fThetaLow, Form("#theta_{e}=%.1f^{#circ}", THETA_LOW), "l");
    leg->AddEntry(fW2, "W^{2} = 4 GeV^{2}", "l");
    leg->AddEntry(fy, "y = 0.75", "l");
    leg->AddEntry(fylow, "y = 0.25", "l");
    leg->Draw();

    c->SaveAs("xQ2_overlay_final.png");

    ///print Q2 outputs 
    // ================== check ouput ==================
    
    std::cout << std::fixed << std::setprecision(3);

    std::cout << "\ntheta_HI = "  << THETA_HI
              << " ; theta_mid = " << THETA_MID
              << " ; theta_lo = "  << THETA_LOW << " ;\n";
              const double loQregN  = std::max(1.0, Q2_of_x_W2(regLMNx_LO, 4)); // x in [0.22, 0.28]
const double loQregQ  = std::max(1.0, Q2_of_x_W2(regOPQx_LO, 4)); // x in [0.28, 0.37]
const double loQregRS =    Q2_of_x_W2(regRSx_LO, 4); // starts at x=0.37 and grows with x
              std::cout << std::fixed << std::setprecision(2)
          << "AB: ["  << regABx_LO  << ", " << regABx_HI  << "]\n"
          << "CDE: [" << regCDEx_LO << ", " << regCDEx_HI << "]\n"
          << "FGH: [" << regFGHx_LO << ", " << regFGHx_HI << "]\n"
          << "IJK: [" << regIJKx_LO << ", " << regIJKx_HI << "]\n"
          << "LMN: [" << regLMNx_LO << ", " << regLMNx_HI << "]\n"
          << "OPQ: [" << regOPQx_LO << ", " << regOPQx_HI << "]\n"
          << "RS:  [" << regRSx_LO  << ", ...]\n";  
const double loQregA = Q2_of_x_theta(regABx_LO, THETA_MID);
const double hiQregA = Q2_of_x_theta(regABx_HI, THETA_HI);
const double loQregB = 1.0; // fixed
const double hiQregB = Q2_of_x_theta(regABx_HI, THETA_MID);

std::cout << std::fixed << std::setprecision(3);
std::cout << "\n===== Region [MANY] Definitions  =====\n";

std::cout << "A LO = Q2(" << regABx_LO << ", " << THETA_MID << ") = " << Q2_of_x_theta(regABx_LO,THETA_MID  ) << "  // loQregA\n";
std::cout << "A HI = Q2(" << regABx_HI << ", y= " << y_fixed << ") = "<< Q2_of_x_y(regABx_HI,y_fixed  ) << "  // hiQregA\n";
std::cout << std::endl;
std::cout << "B LO = Q2(fixed) = 1     // loQregB \n";
std::cout << "B HI = Q2(" << regABx_HI << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regABx_HI,THETA_MID  ) << "  // hiQregB\n";
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
std::cout << "H LO = Q2(fixed) = 1     //loQregH \n";
std::cout << "H HI = Q2(" << regFGHx_HI << ", " << THETA_LOW << ") = "<< Q2_of_x_theta(regFGHx_HI,THETA_LOW  ) << "  // =1hiQregH\n";
std::cout << std::endl;


std::cout << "I LO = Q2(" << regIJKx_LO << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regIJKx_LO,THETA_MID  ) << "  // loQregI\n";
std::cout << "I HI = Q2(" << regIJKx_HI << ", y=" << y_fixed << ") = "<< Q2_of_x_y(regIJKx_HI,y_fixed  ) << "  // hiQregI\n";
std::cout << std::endl;
std::cout << "J LO = Q2(" << regIJKx_LO << ", " << THETA_LOW << ") = "<< Q2_of_x_theta(regIJKx_LO,THETA_LOW  ) << "  // loQregJ\n";
std::cout << "J HI = Q2(" << regIJKx_HI << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regIJKx_HI,THETA_MID  ) << "  // hiQregJ\n";
std::cout << std::endl;
std::cout << "K LO = Q2(fixed) = 1     //loQregK \n";
std::cout << "K HI = Q2(" << regIJKx_HI << ", " << THETA_LOW << ") = "<< Q2_of_x_theta(regIJKx_HI,THETA_LOW  ) << "  // hiQregK\n";
std::cout << std::endl;

std::cout << "L LO = Q2(" << regLMNx_LO << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regLMNx_LO,THETA_MID  ) << "  // loQregL\n";
std::cout << "L HI = Q2(" << regLMNx_HI << ", y=" << y_fixed << ") = "<< Q2_of_x_y(regLMNx_HI,y_fixed  ) << "  // hiQregL\n";
std::cout << std::endl;
std::cout << "M LO = Q2(" << regLMNx_LO << ", " << THETA_LOW << ") = "<< Q2_of_x_theta(regLMNx_LO,THETA_LOW  ) << "  // loQregM\n";
std::cout << "M HI = Q2(" << regLMNx_HI << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regLMNx_HI,THETA_MID  ) << "  // hiQregM\n";
std::cout << std::endl;
std::cout << "N LO = Q2(" << regLMNx_LO << ", y=" << y_fixedlow << ") = " << Q2_of_x_y(regLMNx_LO,y_fixedlow ) << "  // loQregN \n";
std::cout << "N HI = Q2(" << regLMNx_HI << ", " << THETA_LOW << ") = "<< Q2_of_x_theta(regLMNx_HI,THETA_LOW  ) << "  // hiQregN\n";
std::cout << std::endl;

std::cout << "O LO = Q2(" << regOPQx_LO << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regOPQx_LO,THETA_MID  ) << "  // loQregO\n";
std::cout << "O HI = Q2(" << regOPQx_HI << ", y=" << y_fixed << ") = "<< Q2_of_x_y(regOPQx_HI,y_fixed  ) << "  // hiQregO\n";
std::cout << std::endl;
std::cout << "P LO = Q2(" << regOPQx_LO << ", " << THETA_LOW << ") = "<< Q2_of_x_theta(regOPQx_LO,THETA_LOW  ) << "  // loQregP\n";
std::cout << "P HI = Q2(" << regOPQx_HI << ", " << THETA_MID << ") = "<< Q2_of_x_theta(regOPQx_HI,THETA_MID  ) << "  // hiQregP\n";
std::cout << std::endl;
std::cout << "Q LO = Q2(" << regOPQx_LO << ", y=" << y_fixedlow<< ") = " << Q2_of_x_y(regOPQx_LO,y_fixedlow  ) << "  // loQregQ (W-assisted; max with 1)\n";
std::cout << "Q HI = Q2(" << regOPQx_HI << ", " << THETA_LOW << ") = "<< Q2_of_x_theta(regOPQx_HI,THETA_LOW  ) << "  // hiQregQ\n";
std::cout << std::endl;

std::cout << "R LO = Q2(" << regRSx_LO  << ", W2=" << 4 << ") = " << Q2_of_x_W2(regRSx_LO, 4) << "  // loQregRS (W-assisted at x=RS_LO; increases with x)\n";
std::cout << "R HI = non Q2 limit fixed for R-S region    TBD  //loQregRS \n";


    fin->Close();
    return 0;
}



