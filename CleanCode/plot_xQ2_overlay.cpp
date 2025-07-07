#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLine.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <cmath>
#include <iostream>

constexpr double M = 0.9396;        //nucleon M  [GeV]
constexpr double E = 10.5;          //beam

constexpr double THETA_HI  = 28.0;  //top line
constexpr double THETA_MID = 12.0;  //middle line

inline double deg2rad(double d) { return d*M_PI/180.0; }
//defining fct w/ x and th
double Q2_of_x_theta(double x, double thetaDeg)
{
    double c = 1.0 - std::cos(deg2rad(thetaDeg));
    double A = 2.0 * M * E * E * c;
    double B = E * c;
    return A * x / (B + M * x);
}

//inverse function relation   
double x_of_Q2_theta(double Q2, double thetaDeg)
{
    double c = 1.0 - std::cos(deg2rad(thetaDeg));
    return (E * c * Q2) / (M * (2.0 * E * E * c - Q2));
}

// Q2 with x and W2
double Q2_of_x_W2(double x, double W2)
{
    return (W2 - M*M) * x / (1.0 - x);
}




int main()
{
    //sanity checj
    TFile* fin = TFile::Open("DptinputFiles/REcutC2_test.root", "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open input ROOT file.\n";
        return 1;
    }

    auto* h = dynamic_cast<TH2F*>(fin->Get("xQ2_C2_RGD"));
    if (!h) {
        std::cerr << "ERROR: histogram xQ2_C2_RGD not found.\n";
        return 2;
    }

    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c","xQ2 with overlays",1400,800);
    c->SetRightMargin(0.15);

    auto make_theta_curve = [](double thetaDeg, const char* name) {
        double c = 1.0 - std::cos(deg2rad(thetaDeg));
        double A = 2.0 * M * E * E * c;
        double B = E * c;
        TF1* f = new TF1(name, "[0]*x/([1]+0.9396*x)", 0.0, 0.7);
        f->SetParameters(A,B);
        f->SetLineWidth(2);
        f->SetLineColor(kRed);
        f->SetRange(0.0, 0.7);  // x range
        f->SetMinimum(1.0);     // Q2 min
        return f;
    };

    TF1* fThetaHi  = make_theta_curve(THETA_HI , "fThetaHi");
    TF1* fThetaMid = make_theta_curve(THETA_MID, "fThetaMid");
    fThetaMid->SetLineStyle(2);      

    // W2 = 4 curve
    TF1* fW2 = new TF1("fW2", "(4 - 0.9396*0.9396)*x/(1 - x)", 0.0, 0.7);
    fW2->SetLineColor(kOrange+1);
    fW2->SetLineWidth(2);
    fW2->SetMinimum(1.0);  //starts at Q2=1


    //vertical xB = 0.20 and Q2= 1 to Q2 = Q2(x=0.20, θ=THETA_HI) using fct
    double Q2_stop = Q2_of_x_theta(0.20, THETA_HI);
    TLine* line_x02 = new TLine(0.20, 1.0, 0.20, Q2_stop);
    line_x02->SetLineColor(kRed);
    line_x02->SetLineWidth(2);

    //vertical at xB = 0.50 to Q2 = 5
    TLine* line_x05 = new TLine(0.50,Q2_of_x_W2(0.50,4.0),0.50,5.0);
    line_x05->SetLineColor(kRed);
    line_x05->SetLineWidth(2);

    //horizontal Q2 = 1   from x(Q²,THETA_HI) up to 0.20
    double x1_start = x_of_Q2_theta(1.0, THETA_HI);
    TLine* line_Q2_1 = new TLine(x1_start, 1.0, 0.20, 1.0);
    line_Q2_1->SetLineColor(kViolet);
    line_Q2_1->SetLineWidth(2);

    //horizontal Q2 = 5   from xB 0.30 0.50, arbitrary
    TLine* line_Q2_5 = new TLine(0.30, 5.0, 0.50, 5.0);
    line_Q2_5->SetLineColor(kViolet);
    line_Q2_5->SetLineWidth(2);

    h->GetXaxis()->SetTitle("x_{B}");
    h->GetYaxis()->SetTitle("Q^{2}  [GeV^{2}]");
    h->Draw("colz");

    fThetaHi ->Draw("same");
    fThetaMid->Draw("same");
    fW2      ->Draw("same");
    line_x02 ->Draw("same");
    line_x05 ->Draw("same");
    line_Q2_1->Draw("same");
    line_Q2_5->Draw("same");

    TLegend* leg = new TLegend(0.71,0.66,0.93,0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(fThetaHi ,Form("#theta_{e}=%.0f^{#circ}",THETA_HI ),"l");
    leg->AddEntry(fThetaMid,Form("#theta_{e}=%.0f^{#circ}",THETA_MID),"l");
    leg->AddEntry(fW2,"W^{2}=4 GeV^{2}","l");
    leg->Draw();

    c->SaveAs("xQ2_overlay.png");
    fin->Close();
    return 0;
}
