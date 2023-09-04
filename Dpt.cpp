// Include necessary ROOT headers
//Run and compile with 
//g++ -o Dpt Dpt.cpp `root-config --cflags --libs`
// ./Dpt


////////////////////////////////////////////////////////////////
//Code separated from the first one. 
//function  determines the transvers m broadening Dpt from a given histogram found in a root file 
// august 2023
//
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
#include <cmath>
#include "histograms.hpp"
using namespace std;


void calculateWAvg(TH1F* hD, TH1F* h_weightD, TH1F* h_sqD,TH1F* hSn, TH1F* h_weightSn,TH1F* h_sqSn, TGraphErrors* g_result) {
    int numBins = hD->GetNbinsX();

    for (int bin = 1; bin <= numBins; bin++) {
        double N_D = hD->GetBinContent(bin);
        double P_D = h_weightD->GetBinContent(bin);
        double Psq_D = h_sqD->GetBinContent(bin);
        double N_Sn = hSn->GetBinContent(bin);
        double P_Sn = h_weightSn->GetBinContent(bin);
        double Psq_Sn = h_sqSn->GetBinContent(bin);
        double x = hD->GetXaxis()->GetBinCenter(bin);
        
        double wavg_D  = (N_D > 0) ? P_D/N_D : 0.0;
        double wavg_Sn = (N_Sn > 0) ? P_Sn/N_Sn : 0.0;
        double dpt_point = wavg_Sn - wavg_D ; 
        double variance_D = Psq_D/N_D - pow(wavg_D,2);     
        double variance_Sn=   Psq_Sn/N_Sn - pow(wavg_Sn,2) ;
        double Err_D  = sqrt(variance_D/N_D);
        double Err_Sn = sqrt(variance_Sn/N_Sn);
        double unc = sqrt( pow( Err_D,2) + pow(Err_Sn ,2) );

        g_result->SetPoint(bin - 1, x, dpt_point);
        g_result->SetPointError(bin - 1, 0, unc);
    }
}




int main(){
h_Dpt hist_Dpt;

    TCanvas* cDpt=new TCanvas("Dptplot","Dptplot"); //create new canvas
    TCanvas* cDpt3=new TCanvas("Dptplot3","Dptplot3"); //create new canvas
    TCanvas* ctest=new TCanvas("test","test"); //create new canvas
    

    TGraphErrors *Dpt_Qter = new TGraphErrors();
    TGraphErrors *Dpt_vter = new TGraphErrors();
    TGraphErrors *Dpt_zter = new TGraphErrors();
    cDpt->Divide(2,2);

    TFile* fileD = new TFile("../files2read/REoutput_D.root", "READ");
    TFile* fileSn = new TFile("../files2read/REoutput_Sn.root", "READ");

    
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

    TH1F* h_QweightD = dynamic_cast<TH1F*>(fileD->Get("weightedQp"));
    TH1F* h_QweightSn = dynamic_cast<TH1F*>(fileSn->Get("weightedQp"));
    TH1F* h_vweightD = dynamic_cast<TH1F*>(fileD->Get("weightedvp"));
    TH1F* h_vweightSn = dynamic_cast<TH1F*>(fileSn->Get("weightedvp"));
    TH1F* h_zweightD = dynamic_cast<TH1F*>(fileD->Get("weightedzp"));
    TH1F* h_zweightSn = dynamic_cast<TH1F*>(fileSn->Get("weightedzp"));

    TH1F* h_QsqD = dynamic_cast<TH1F*>(fileSn->Get("VarianceQp"));
    TH1F* h_QsqSn = dynamic_cast<TH1F*>(fileSn->Get("VarianceQp"));
    TH1F* h_vsqD = dynamic_cast<TH1F*>(fileD->Get("Variancevp"));
    TH1F* h_vsqSn = dynamic_cast<TH1F*>(fileSn->Get("Variancevp"));
    TH1F* h_zsqD = dynamic_cast<TH1F*>(fileD->Get("Variancezp"));
    TH1F* h_zsqSn = dynamic_cast<TH1F*>(fileSn->Get("Variancezp"));

    //TH2F* h_2DtestD  = dynamic_cast<TH2F*>(fileD->Get("pt2Q2gen"));
    //TH2F* h_2DtestSn = dynamic_cast<TH2F*>(fileSn->Get("pt2Q2gen"));
    //maybe not even need to dynamic cast, maybe hjust enough with calling the class
    //-> Yes, dynamic cast is needed, to make the difference between Sn and De
    //-> an adjustment has to be done to sumW2 histograms in order to properly recover them making the distinction on them. 
  
    
    int numBinsQ = h2DQ_Sn->GetNbinsY();
    int numBinsP = h2DQ_Sn->GetNbinsX();
    int numBinsv = h2Dv_D->GetNbinsY();
    int numBinsz = h2Dz_D->GetNbinsY();
    int nubinsQ = h_QtestD->GetNbinsX();
    int nubinsv = h_vtestD->GetNbinsX();
    int nubinsz = h_ztestD->GetNbinsX();

    calculateWAvg(h_QtestD, h_QweightD, h_QsqD, h_QtestSn, h_QweightSn, h_QsqSn, Dpt_Qter);
    calculateWAvg(h_vtestD, h_vweightD, h_vsqD, h_vtestSn, h_vweightSn, h_vsqSn, Dpt_vter);
    calculateWAvg(h_ztestD, h_zweightD, h_zsqD, h_ztestSn, h_zweightSn, h_zsqSn, Dpt_zter);

    // Close the first ROOT file
    fileD->Close();
    fileSn->Close();
    cDpt->cd(1);
    //Dpt_Q->SetTitle(title);
    Dpt_Qter->GetXaxis()->SetTitle("Q{2} (GeV^{2}) " );
    Dpt_Qter->GetYaxis()->SetTitle("P_{t}^{2}diff) ");
    Dpt_Qter->Draw("APE");
    cDpt->cd(2);
    //Dpt_v->SetTitle(title);
    Dpt_vter->GetXaxis()->SetTitle("v " );
    Dpt_vter->GetYaxis()->SetTitle("P_{t}^{2}diff) ");
    Dpt_vter->Draw("APE");
    cDpt->cd(3);
    //Dpt_z->SetTitle(title);
    Dpt_zter->GetXaxis()->SetTitle("z  " );
    Dpt_zter->GetYaxis()->SetTitle("P_{t}^{2}diff) ");
    Dpt_zter->Draw("APE");
    //h_QtestD->Draw("COLZ");
    cDpt->cd(4);
    //Dpt_Qter->GetXaxis()->SetTitle("null  " );
    //Dpt_Qter->GetYaxis()->SetTitle("null ");
    Dpt_Qter->Draw("APE");
    cDpt->SaveAs("dpttest.pdf");
    cDpt->SaveAs("dpttest.root");
}
