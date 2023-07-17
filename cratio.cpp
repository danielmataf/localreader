// Include necessary ROOT headers
//Run and compile with 
//g++ -o cratio cratio.cpp `root-config --cflags --libs`
// ./cratio

////////////////////////////////////////////////////////////////
//Code separated from the first one. 
// June 2023 
// Need to test " plot cratio.pdf" and "testcratio.pdf" recovered from recoveryposCLAS.cpp
//  cratio = <cphi>_Sn/ <cphi>_D        //cratio(Q2,v,z,pt)
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
using namespace std;



int main() {
    TCanvas* cR=new TCanvas("ratio","ratio"); //create new canvas
    cR->Divide(2,2);
    TGraphErrors *cos_Q = new TGraphErrors();
    TGraphErrors *cos_v = new TGraphErrors();
    TGraphErrors *cos_z = new TGraphErrors();
    TGraphErrors *cos_pt = new TGraphErrors();
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
    TH2F* h2DQ_D = dynamic_cast<TH2F*>(fileD->Get("cosQ"));  // Retrieve the TH2F histogram
    TH2F* h2Dv_D = dynamic_cast<TH2F*>(fileD->Get("cosv"));  
    TH2F* h2Dz_D = dynamic_cast<TH2F*>(fileD->Get("cosz"));
    TH2F* h2Dp_D = dynamic_cast<TH2F*>(fileD->Get("cosp"));
    TH2F* h2DQ_Sn = dynamic_cast<TH2F*>(fileSn->Get("cosQ"));
    TH2F* h2Dv_Sn = dynamic_cast<TH2F*>(fileSn->Get("cosv"));
    TH2F* h2Dz_Sn = dynamic_cast<TH2F*>(fileSn->Get("cosz"));
    TH2F* h2Dp_Sn = dynamic_cast<TH2F*>(fileSn->Get("cosp"));

    int numBinsQ = h2DQ_D->GetNbinsY();
    int numBinsc = h2DQ_D->GetNbinsX();

    for (int binQ = 1; binQ <= numBinsQ; binQ++) {
        double sumc = 0.0;
        double sumWeights = 0.0;
        double sumc_Sn = 0.0;
        double sumWeights_Sn = 0.0;

        for (int binc = 1; binc <= numBinsc; binc++) {
            double c = h2DQ_D->GetXaxis()->GetBinCenter(binc);
            double weight = h2DQ_D->GetBinContent(binc, binQ);
            double c_Sn = h2DQ_Sn->GetXaxis()->GetBinCenter(binc);
            double weight_Sn = h2DQ_Sn->GetBinContent(binc, binQ);
            sumc += c * weight;
            sumWeights += weight;
            sumc_Sn += c_Sn * weight_Sn;
            sumWeights_Sn += weight_Sn;

                

        }
        double averagec = sumc / sumWeights;
        double averagec_Sn = sumc_Sn / sumWeights_Sn;
    }



    int numBinsv = h2Dv_D->GetNbinsY();

    for (int binv = 1; binv <= numBinsv; binv++) {
        double sumc = 0.0;
        double sumWeightsv = 0.0;
        double sumc_Sn = 0.0;
        double sumWeightsv_Sn = 0.0;

        for (int binc = 1; binc <= numBinsc; binc++) {
            double c = h2Dv_D->GetXaxis()->GetBinCenter(binc);
            double weightv = h2Dv_D->GetBinContent(binc, binv);
            double c_Sn = h2Dv_Sn->GetXaxis()->GetBinCenter(binv);
            double weightv_Sn = h2Dv_Sn->GetBinContent(binc, binv);
            sumc += c * weightv;
            sumWeightsv += weightv;
            sumc_Sn += c_Sn * weightv_Sn;
            sumWeightsv_Sn += weightv_Sn;

                

        }
        double averagecv = sumc / sumWeightsv;
        double averagecv_Sn = sumc_Sn / sumWeightsv_Sn;
    }


    int numBinsz = h2Dz_D->GetNbinsY();

    for (int binz = 1; binz <= numBinsz; binz++) {
        double sumc = 0.0;
        double sumWeightsz = 0.0;
        double sumc_Sn = 0.0;
        double sumWeightsz_Sn = 0.0;

        for (int binc = 1; binc <= numBinsc; binc++) {
            double c = h2Dz_D->GetXaxis()->GetBinCenter(binc);
            double weightz = h2Dz_D->GetBinContent(binc, binz);
            double c_Sn = h2Dz_Sn->GetXaxis()->GetBinCenter(binz);
            double weightz_Sn = h2Dz_Sn->GetBinContent(binc, binz);
            sumc += c * weightz;
            sumWeightsz += weightz;
            sumc_Sn += c_Sn * weightz_Sn;
            sumWeightsz_Sn += weightz_Sn;

                

        }
        double averagecz = sumc / sumWeightsz;
        double averagecz_Sn = sumc_Sn / sumWeightsz_Sn;
        cout<<averagecz_Sn<< " / " <<averagecz<<" = "<<averagecz_Sn / averagecz<< endl;
    }

    int numBinsp = h2Dp_D->GetNbinsY();

    for (int binp = 1; binp <= numBinsp; binp++) {
        double sumc = 0.0;
        double sumWeightsp = 0.0;
        double sumc_Sn = 0.0;
        double sumWeightsp_Sn = 0.0;

        for (int binc = 1; binc <= numBinsc; binc++) {
            double c = h2Dp_D->GetXaxis()->GetBinCenter(binc);
            double weightp = h2Dp_D->GetBinContent(binc, binp);
            double c_Sn = h2Dp_Sn->GetXaxis()->GetBinCenter(binp);
            double weightp_Sn = h2Dp_Sn->GetBinContent(binc, binp);
            sumc += c * weightp;
            sumWeightsp += weightp;
            sumc_Sn += c_Sn * weightp_Sn;
            sumWeightsp_Sn += weightp_Sn;

                

        }
        double averagecp = sumc / sumWeightsp;
        double averagecp_Sn = sumc_Sn / sumWeightsp_Sn;
        //cout<<averagecp_Sn<< " / " <<averagecp<<" = "<<averagecp_Sn / averagecp<< endl;
    }

    // Close the first ROOT file
    fileD->Close();
    fileSn->Close();
    

    /*
cR->cd(1);
cos_Q->GetXaxis()->SetTitleSize(0.05);
cos_Q->GetYaxis()->SetTitleSize(0.05);
cos_Q->SetMarkerSize(0.5);
cos_Q->SetMarkerStyle(21);
cos_Q->GetXaxis()->SetRangeUser(1,6.5);
cos_Q->GetXaxis()->SetTitle("Q^{2} " );
cos_Q->GetYaxis()->SetTitle("<cos(#Phi)>_{Sn} / <cos(#Phi)>_{De}");
cos_Q->Draw("AP");
cR->cd(2);
cos_v->GetXaxis()->SetTitleSize(0.05);
cos_v->GetYaxis()->SetTitleSize(0.05);
cos_v->GetXaxis()->SetRangeUser(2,9);
cos_v->SetMarkerSize(0.5);
cos_v->SetMarkerStyle(21);
cos_v->GetXaxis()->SetTitle("#nu " );
cos_v->GetYaxis()->SetTitle("<cos(#Phi)>_{Sn} / <cos(#Phi)>_{De}");
cos_v->Draw("AP");
cR->cd(3);
cos_z->GetXaxis()->SetTitleSize(0.05);
cos_z->GetYaxis()->SetTitleSize(0.05);
cos_z->SetMarkerSize(0.5);
cos_z->SetMarkerStyle(21);
cos_z->GetXaxis()->SetTitle("p_{t}^{2} " );
cos_z->GetYaxis()->SetTitle("<cos(#Phi)>_{Sn} / <cos(#Phi)>_{De}");
cos_z->Draw("AP");
cR->cd(4);
cos_pt->GetXaxis()->SetTitleSize(0.05);
cos_pt->GetYaxis()->SetTitleSize(0.05);
cos_pt->SetMarkerSize(0.5);
cos_pt->SetMarkerStyle(21);
cos_pt->GetXaxis()->SetRangeUser(0.2,0.9);
cos_pt->GetXaxis()->SetTitle("z " );
cos_pt->GetYaxis()->SetTitle("<cos(#Phi)>_{Sn} / <cos(#Phi)>_{De}");
cos_pt->Draw("AP");
    cR->SaveAs("cratio.pdf");
    cR->SaveAs("cratio.root");
*/
    return 0;
}

