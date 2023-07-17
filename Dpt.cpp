// Include necessary ROOT headers
//Run and compile with 
//g++ -o Dpt Dpt.cpp `root-config --cflags --libs`
// ./Dpt


////////////////////////////////////////////////////////////////
//Code separated from the first one. 
//It contains only a function that determines (maybe) the transvers m broadening Dpt from a given histogram found in a root file 
// June 2023
//TBD: ERRORS!!!
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







int main() {
    h_Dpt hist_Dpt;
    TCanvas* cDpt=new TCanvas("Dptplot","Dptplot"); //create new canvas
    TCanvas* ctest=new TCanvas("test","test"); //create new canvas
    cDpt->Divide(2,2);
    TGraphErrors *Dpt_Q = new TGraphErrors();   //creation of emmpty graphs to be filled at the end 
    TGraphErrors *Dpt_v = new TGraphErrors();
    TGraphErrors *Dpt_z = new TGraphErrors();
    // Open the first ROOT file
    TFile* fileD = new TFile("../files2read/output_D.root", "READ");
    TFile* fileSn = new TFile("../files2read/output_Sn.root", "READ");
    // Retrieve the histograms from the ROOT file
    //TH1F* h_Q1test = dynamic_cast<TH1F*>(file1->Get("Q2"));
    
    TH1F* h_QtestD = dynamic_cast<TH1F*>(fileD->Get("Q2"));
    TH1F* h_QtestSn = dynamic_cast<TH1F*>(fileSn->Get("Q2"));
    TH1F* h_vtestD = dynamic_cast<TH1F*>(fileD->Get("gamnu"));
    TH1F* h_vtestSn = dynamic_cast<TH1F*>(fileSn->Get("gamnu"));
    TH1F* h_ztestD = dynamic_cast<TH1F*>(fileD->Get("z"));
    TH1F* h_ztestSn = dynamic_cast<TH1F*>(fileSn->Get("z"));
    TH1F* h_pttestD = dynamic_cast<TH1F*>(fileD->Get("P_t"));
    TH1F* h_pttestSn = dynamic_cast<TH1F*>(fileSn->Get("P_t"));
    TH2F* h2DQ_D = dynamic_cast<TH2F*>(fileD->Get("pt2Q2gen"));  // Retrieve the TH2F histogram
    TH2F* h2Dv_D = dynamic_cast<TH2F*>(fileD->Get("pt2v2gen"));  
    TH2F* h2Dz_D = dynamic_cast<TH2F*>(fileD->Get("pt2z2gen"));
    TH2F* h2DQ_Sn = dynamic_cast<TH2F*>(fileSn->Get("pt2Q2gen"));
    TH2F* h2Dv_Sn = dynamic_cast<TH2F*>(fileSn->Get("pt2v2gen"));
    TH2F* h2Dz_Sn = dynamic_cast<TH2F*>(fileSn->Get("pt2z2gen"));

    //TH2F* h_2DtestD  = dynamic_cast<TH2F*>(fileD->Get("pt2Q2gen"));
    //TH2F* h_2DtestSn = dynamic_cast<TH2F*>(fileSn->Get("pt2Q2gen"));
    //maybe not even need to dynamic cast, maybe hjust enough with calling the class
    //-> Yes, dynamic cast is needed, to make the difference between Sn and De
    //-> an adjustment has to be done to sumW2 histograms in order to properly recover them making the distinction on them. 
  
    
    int numBinsQ = h2DQ_Sn->GetNbinsY();
    int numBinsP = h2DQ_Sn->GetNbinsX();
/*
    for (int binQ = 1; binQ <= numBinsQ; binQ++) {
        double sumP = 0.0;
        double sumWeights = 0.0;
        double sumP_Sn = 0.0;
        double sumWeights_Sn = 0.0;

        for (int binP = 1; binP <= numBinsP; binP++) {
            double p = h2DQ_D->GetXaxis()->GetBinCenter(binP);
            double weight = h2DQ_D->GetBinContent(binP, binQ);
            double p_Sn = h2DQ_Sn->GetXaxis()->GetBinCenter(binP);
            double weight_Sn = h2DQ_Sn->GetBinContent(binP, binQ);
            sumP += p * weight;
            sumWeights += weight;
            sumP_Sn += p_Sn * weight_Sn;
            sumWeights_Sn += weight_Sn;

                

        }
        double averageP = sumP / sumWeights;
        double averageP_Sn = sumP_Sn / sumWeights_Sn;
        //cout<<averageP_Sn<< " - " <<averageP<<" = "<<averageP_Sn - averageP<< endl;
    }



    int numBinsv = h2Dv_D->GetNbinsY();

    for (int binv = 1; binv <= numBinsv; binv++) {
        double sumP = 0.0;
        double sumWeightsv = 0.0;
        double sumP_Sn = 0.0;
        double sumWeightsv_Sn = 0.0;

        for (int binP = 1; binP <= numBinsP; binP++) {
            double p = h2Dv_D->GetXaxis()->GetBinCenter(binP);
            double weightv = h2Dv_D->GetBinContent(binP, binv);
            double p_Sn = h2Dv_Sn->GetXaxis()->GetBinCenter(binv);
            double weightv_Sn = h2Dv_Sn->GetBinContent(binP, binv);
            sumP += p * weightv;
            sumWeightsv += weightv;
            sumP_Sn += p_Sn * weightv_Sn;
            sumWeightsv_Sn += weightv_Sn;

                

        }
        double averagePv = sumP / sumWeightsv;
        double averagePv_Sn = sumP_Sn / sumWeightsv_Sn;
    }


    int numBinsz = h2Dz_D->GetNbinsY();

    for (int binz = 1; binz <= numBinsz; binz++) {
        double sumP = 0.0;
        double sumWeightsz = 0.0;
        double sumP_Sn = 0.0;
        double sumWeightsz_Sn = 0.0;

        for (int binP = 1; binP <= numBinsP; binP++) {
            double p = h2Dz_D->GetXaxis()->GetBinCenter(binP);
            double weightz = h2Dz_D->GetBinContent(binP, binz);
            double p_Sn = h2Dz_Sn->GetXaxis()->GetBinCenter(binz);
            double weightz_Sn = h2Dz_Sn->GetBinContent(binP, binz);
            sumP += p * weightz;
            sumWeightsz += weightz;
            sumP_Sn += p_Sn * weightz_Sn;
            sumWeightsz_Sn += weightz_Sn;

                

        }
        double averagePz = sumP / sumWeightsz;
        double averagePz_Sn = sumP_Sn / sumWeightsz_Sn;
        cout<<averagePz_Sn<< " - " <<averagePz<<" = "<<averagePz_Sn - averagePz<< endl;
    }

    // Close the first ROOT file
    fileD->Close();
    fileSn->Close();
    cDpt->cd(1);
    //R_Q->SetTitle(title);
    //R_Q->GetXaxis()->SetTitle("Q{2} (GeV^{2}) " );
    //ctest->GetYaxis()->SetTitle("R^{#Pi+}_{Sn}");
    hist_Dpt.h_pQ_De->Draw("COLZ");
    cDpt->cd(2);
    hist_Dpt.h_pv_De->Draw("COLZ");
    cDpt->cd(3);
    //h_QtestD->Draw("COLZ");
    cDpt->cd(4);
    //h_2DtestD->Draw("COLZ");
    cDpt->SaveAs("dpttest.pdf");
    cDpt->SaveAs("dpttest.root");
*/
    return 0;
}

