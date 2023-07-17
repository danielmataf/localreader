// Include necessary ROOT headers
//Run and compile with 
//g++ -o testsratio sratio.cpp `root-config --cflags --libs`
// ./testsratio

////////////////////////////////////////////////////////////////
//Code separated from the first one. 
// June 2023 
// Need to test " plot sRatio.pdf" and "testsRatio.pdf" recovered from recoveryposCLAS.cpp
//  sratio = <sphi>_Sn/ <sphi>_D        //sRatio(Q2,v,z,pt)
//Warning: for it to work properly you will need to change the binning of the histogram found in the root file; 
//  remember binning = 10 in order to have only 10 points (ideal)
////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TH1.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <TCanvas.h>
#include <TH2.h>
#include <TF1.h>
#include <TRatioPlot.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TLine.h>
#include <TArrow.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#include "histograms.hpp"
using namespace std;


            //w_sQ_Sn->Fill(Q2,sin(phih));


int main() {
    h_sratio hist_sratio;
    TCanvas* sR=new TCanvas("ratio","ratio"); //sReate new canvas
    sR->Divide(2,2);
    TGraphErrors *sin_Q = new TGraphErrors();
    TGraphErrors *sin_v = new TGraphErrors();
    TGraphErrors *sin_z = new TGraphErrors();
    TGraphErrors *sin_pt = new TGraphErrors();
    // Open the first ROOT file
    //TFile* file1 = new TFile("build/output1.root", "READ");
    TFile* fileD = new TFile("output_D.root", "READ");
    TFile* fileSn = new TFile("output_Sn.root", "READ");

    // Retrieve the first histogram from the ROOT file
    //TH1F* h_Q1test = dynamic_cast<TH1F*>(file1->Get("Q2"));
    TH1F* h_QtestD = dynamic_cast<TH1F*>(fileD->Get("Q2"));
    TH1F* h_QtestSn = dynamic_cast<TH1F*>(fileSn->Get("Q2"));
    TH1F* h_vtestD = dynamic_cast<TH1F*>(fileD->Get("gamnu"));
    TH1F* h_vtestSn = dynamic_cast<TH1F*>(fileSn->Get("gamnu"));
    TH1F* h_ztestD = dynamic_cast<TH1F*>(fileD->Get("z"));
    TH1F* h_ztestSn = dynamic_cast<TH1F*>(fileSn->Get("z"));
    TH1F* h_pttestD = dynamic_cast<TH1F*>(fileD->Get("P_t"));
    TH1F* h_pttestSn = dynamic_cast<TH1F*>(fileSn->Get("P_t"));
    TH1F* h_phiSn = dynamic_cast<TH1F*>(fileSn->Get("phih"));
    TH1F* h_phiD = dynamic_cast<TH1F*>(fileD->Get("phih"));
    TH2F* h2DQ_D = dynamic_cast<TH2F*>(fileD->Get("sinQ"));  // Retrieve the TH2F histogram
    TH2F* h2Dv_D = dynamic_cast<TH2F*>(fileD->Get("sinv"));  
    TH2F* h2Dz_D = dynamic_cast<TH2F*>(fileD->Get("sinz"));
    TH2F* h2Dp_D = dynamic_cast<TH2F*>(fileD->Get("sinp"));
    TH2F* h2DQ_Sn = dynamic_cast<TH2F*>(fileSn->Get("sinQ"));
    TH2F* h2Dv_Sn = dynamic_cast<TH2F*>(fileSn->Get("sinv"));
    TH2F* h2Dz_Sn = dynamic_cast<TH2F*>(fileSn->Get("sinz"));
    TH2F* h2Dp_Sn = dynamic_cast<TH2F*>(fileSn->Get("sinp"));

    int numBinsQ = h2DQ_D->GetNbinsY();
    int numBinsc = h2DQ_D->GetNbinsX();

    for (int binQ = 1; binQ <= numBinsQ; binQ++) {
        double sums = 0.0;
        double sumWeights = 0.0;
        double sums_Sn = 0.0;
        double sumWeights_Sn = 0.0;

        for (int binc = 1; binc <= numBinsc; binc++) {
            double s = h2DQ_D->GetXaxis()->GetBinCenter(binc);
            double weight = h2DQ_D->GetBinContent(binc, binQ);
            double s_Sn = h2DQ_Sn->GetXaxis()->GetBinCenter(binc);
            double weight_Sn = h2DQ_Sn->GetBinContent(binc, binQ);
            sums += s * weight;
            sumWeights += weight;
            sums_Sn += s_Sn * weight_Sn;
            sumWeights_Sn += weight_Sn;

                

        }
        double averages = sums / sumWeights;
        double averages_Sn = sums_Sn / sumWeights_Sn;
    }



    int numBinsv = h2Dv_D->GetNbinsY();

    for (int binv = 1; binv <= numBinsv; binv++) {
        double sums = 0.0;
        double sumWeightsv = 0.0;
        double sums_Sn = 0.0;
        double sumWeightsv_Sn = 0.0;

        for (int binc = 1; binc <= numBinsc; binc++) {
            double c = h2Dv_D->GetXaxis()->GetBinCenter(binc);
            double weightv = h2Dv_D->GetBinContent(binc, binv);
            double c_Sn = h2Dv_Sn->GetXaxis()->GetBinCenter(binv);
            double weightv_Sn = h2Dv_Sn->GetBinContent(binc, binv);
            sums += c * weightv;
            sumWeightsv += weightv;
            sums_Sn += c_Sn * weightv_Sn;
            sumWeightsv_Sn += weightv_Sn;

                

        }
        double averagesv = sums / sumWeightsv;
        double averagesv_Sn = sums_Sn / sumWeightsv_Sn;
    }


    int numBinsz = h2Dz_D->GetNbinsY();

    for (int binz = 1; binz <= numBinsz; binz++) {
        double sums = 0.0;
        double sumWeightsz = 0.0;
        double sums_Sn = 0.0;
        double sumWeightsz_Sn = 0.0;

        for (int binc = 1; binc <= numBinsc; binc++) {
            double c = h2Dz_D->GetXaxis()->GetBinCenter(binc);
            double weightz = h2Dz_D->GetBinContent(binc, binz);
            double c_Sn = h2Dz_Sn->GetXaxis()->GetBinCenter(binz);
            double weightz_Sn = h2Dz_Sn->GetBinContent(binc, binz);
            sums += c * weightz;
            sumWeightsz += weightz;
            sums_Sn += c_Sn * weightz_Sn;
            sumWeightsz_Sn += weightz_Sn;

                

        }
        double averagesz = sums / sumWeightsz;
        double averagesz_Sn = sums_Sn / sumWeightsz_Sn;
        cout<<averagesz_Sn<< " / " <<averagesz<<" = "<<averagesz_Sn / averagesz<< endl;
    }

    int numBinsp = h2Dp_D->GetNbinsY();

    for (int binp = 1; binp <= numBinsp; binp++) {
        double sums = 0.0;
        double sumWeightsp = 0.0;
        double sums_Sn = 0.0;
        double sumWeightsp_Sn = 0.0;

        for (int binc = 1; binc <= numBinsc; binc++) {
            double c = h2Dp_D->GetXaxis()->GetBinCenter(binc);
            double weightp = h2Dp_D->GetBinContent(binc, binp);
            double c_Sn = h2Dp_Sn->GetXaxis()->GetBinCenter(binp);
            double weightp_Sn = h2Dp_Sn->GetBinContent(binc, binp);
            sums += c * weightp;
            sumWeightsp += weightp;
            sums_Sn += c_Sn * weightp_Sn;
            sumWeightsp_Sn += weightp_Sn;

                

        }
        double averagesp = sums / sumWeightsp;
        double averagesp_Sn = sums_Sn / sumWeightsp_Sn;
        //cout<<averagesp_Sn<< " / " <<averagesp<<" = "<<averagesp_Sn / averagesp<< endl;
    }






    // Close the first ROOT file
    fileD->Close();
    fileSn->Close();
    /*
sR->cd(1);
sin_Q->GetXaxis()->SetTitleSize(0.05);
sin_Q->GetYaxis()->SetTitleSize(0.05);
sin_Q->SetMarkerSize(0.5);
sin_Q->SetMarkerStyle(21);
sin_Q->GetXaxis()->SetRangeUser(1,6.5);
sin_Q->GetXaxis()->SetTitle("Q^{2} " );
sin_Q->GetYaxis()->SetTitle("<sin(#Phi)>_{Sn} / <sin(#Phi)>_{D}");
sin_Q->Draw("AP");

sR->cd(2);
sin_v->GetXaxis()->SetTitleSize(0.05);
sin_v->GetYaxis()->SetTitleSize(0.05);
sin_v->GetXaxis()->SetRangeUser(2,9);
sin_v->SetMarkerSize(0.5);
sin_v->SetMarkerStyle(21);
sin_v->GetXaxis()->SetTitle("#nu " );
sin_v->GetYaxis()->SetTitle("<sin(#Phi)>_{Sn} / <sin(#Phi)>_{D}");
//sin_v->Draw("AP");

sR->cd(3);
sin_z->GetXaxis()->SetTitleSize(0.05);
sin_z->GetYaxis()->SetTitleSize(0.05);
sin_z->SetMarkerSize(0.5);
sin_z->SetMarkerStyle(21);
sin_z->GetXaxis()->SetTitle("p_{t}^{2} " );
sin_z->GetYaxis()->SetTitle("<sin(#Phi)>_{Sn} / <sin(#Phi)>_{D}");
//sin_z->Draw("AP");
sR->cd(4);
sin_pt->GetXaxis()->SetTitleSize(0.05);
sin_pt->GetYaxis()->SetTitleSize(0.05);
sin_pt->SetMarkerSize(0.5);
sin_pt->SetMarkerStyle(21);
sin_pt->GetXaxis()->SetRangeUser(0.2,0.9);
sin_pt->GetXaxis()->SetTitle("z " );
sin_pt->GetYaxis()->SetTitle("<sin(#Phi)>_{Sn} / <sin(#Phi)>_{D}");
sin_pt->Draw("AP");

    sR->SaveAs("testsRatio.pdf");
    sR->SaveAs("testsRatio.root");
*/
    return 0;

}

