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

void calculateSRat(TH1F* hD, TH1F* h_weightD, TH1F* h_sqD,TH1F* hSn, TH1F* h_weightSn,TH1F* h_sqSn, TGraphErrors* g_result) {
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
        double srat_point = (wavg_D > 0) ? wavg_Sn/wavg_D : 0.0; 
        double variance_D = Psq_D/N_D - pow(wavg_D,2);     
        double variance_Sn=   Psq_Sn/N_Sn - pow(wavg_Sn,2) ;
        double Err_D  = sqrt(variance_D/N_D);
        double Err_Sn = sqrt(variance_Sn/N_Sn);
        double unc = (wavg_Sn/wavg_D) *sqrt( pow( Err_D/wavg_D,2) + pow(Err_Sn/wavg_Sn ,2) );
        if (unc>5){
            unc = 0.0;
            srat_point = 0.0;
        }

        g_result->SetPoint(bin - 1, x, srat_point);
        g_result->SetPointError(bin - 1, 0, unc);
            //w_sQ_Sn->Fill(Q2,sin(phih));

    }
}
int main() {
    h_sratio hist_sratio;
    TCanvas* sR=new TCanvas("ratio","ratio"); //sReate new canvas
    sR->Divide(2,2);
    TGraphErrors *sin_Q = new TGraphErrors();
    TGraphErrors *sin_v = new TGraphErrors();
    TGraphErrors *sin_z = new TGraphErrors();
    TGraphErrors *sin_pt = new TGraphErrors();
    TGraphErrors *sin_Qbis= new TGraphErrors();
    TGraphErrors *sin_vbis= new TGraphErrors();
    TGraphErrors *sin_zbis= new TGraphErrors();
    TGraphErrors *sin_pbis = new TGraphErrors();
    TGraphErrors *sin_Qter = new TGraphErrors();
    TGraphErrors *sin_vter = new TGraphErrors();
    TGraphErrors *sin_zter = new TGraphErrors();
    TGraphErrors *sin_pter = new TGraphErrors();


    // Open the first ROOT file
    //TFile* file1 = new TFile("build/output1.root", "READ");
    TFile* fileD = new TFile("../files2read/REoutput_D.root", "READ");
    TFile* fileSn = new TFile("../files2read/REoutput_Sn.root", "READ");

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

    TH1F* h_QweightD = dynamic_cast<TH1F*>(fileD->Get("weightedQsin"));
    TH1F* h_QweightSn = dynamic_cast<TH1F*>(fileSn->Get("weightedQsin"));
    TH1F* h_vweightD = dynamic_cast<TH1F*>(fileD->Get("weightedvsin"));
    TH1F* h_vweightSn = dynamic_cast<TH1F*>(fileSn->Get("weightedvsin"));
    TH1F* h_zweightD = dynamic_cast<TH1F*>(fileD->Get("weightedzsin"));
    TH1F* h_zweightSn = dynamic_cast<TH1F*>(fileSn->Get("weightedzsin"));
    TH1F* h_ptweightD = dynamic_cast<TH1F*>(fileD->Get("weightedpsin"));
    TH1F* h_ptweightSn = dynamic_cast<TH1F*>(fileSn->Get("weightedpsin"));

    TH1F* h_QsqD = dynamic_cast<TH1F*>(fileSn->Get("VarianceQsin"));
    TH1F* h_QsqSn = dynamic_cast<TH1F*>(fileSn->Get("VarianceQsin"));
    TH1F* h_vsqD = dynamic_cast<TH1F*>(fileD->Get("Variancevsin"));
    TH1F* h_vsqSn = dynamic_cast<TH1F*>(fileSn->Get("Variancevsin"));
    TH1F* h_zsqD = dynamic_cast<TH1F*>(fileD->Get("Variancezsin"));
    TH1F* h_zsqSn = dynamic_cast<TH1F*>(fileSn->Get("Variancezsin"));
    TH1F* h_ptsqD = dynamic_cast<TH1F*>(fileD->Get("Variancepsin"));
    TH1F* h_ptsqSn = dynamic_cast<TH1F*>(fileSn->Get("Variancepsin"));


    int numBinsQ = h2DQ_D->GetNbinsY();
    int numbinsS = h2DQ_D->GetNbinsX();
    int numBinsv = h2Dv_D->GetNbinsY();
    int numBinsp = h2Dp_D->GetNbinsY();
    int numBinsz = h2Dz_D->GetNbinsY();
    int nubinsQ = h_QtestD->GetNbinsX(); 
    int nubinsv = h_vtestD->GetNbinsX(); 
    int nubinsz = h_ztestD->GetNbinsX(); 
    int nubinsp = h_pttestD->GetNbinsX(); 

    calculateSRat( h_QtestD, h_QweightD, h_QsqD,h_QtestSn, h_QweightSn,h_QsqD, sin_Q);
    calculateSRat( h_vtestD, h_vweightD, h_vsqD,h_vtestSn, h_vweightSn,h_vsqD, sin_v);
    calculateSRat( h_ztestD, h_zweightD, h_zsqD,h_ztestSn, h_zweightSn,h_zsqD, sin_z);
    calculateSRat( h_pttestD, h_ptweightD, h_ptsqD,h_pttestSn, h_ptweightSn,h_ptsqD, sin_pt);

/*


    /////////fix weight test/////////////////
    /////////////////////////////////////////
    for (int binQ = 1; binQ <= nubinsQ; binQ++) {
        double N_D = h_QtestD->GetBinContent(binQ);
        double sin_D = h_QweightD->GetBinContent(binQ);
        double N_Sn = h_QtestSn->GetBinContent(binQ);        
        double sin_Sn = h_QweightSn->GetBinContent(binQ);  
        double x_Q_Sn = h_QtestSn->GetXaxis()->GetBinCenter(binQ);
        double Psq_D = h_QsqD->GetBinContent(binQ);
        double Psq_Sn = h_QsqSn->GetBinContent(binQ);


        //<sin (phi_D)> = sin_D/N_D
        double wavg_D  = (N_D > 0) ? sin_D/N_D : 0.0;
        double wavg_Sn = (N_Sn > 0) ? sin_Sn/N_Sn : 0.0;
        double sratio_point = (wavg_D > 0) ? wavg_Sn/wavg_D : 0.0; 
        //cout<<wavg_Sn<<"/" <<wavg_D <<" = " <<cratio_point<<"   cratiopoint"<< endl;

        double varianceQ_D = abs(Psq_D/N_D - pow(wavg_D,2));     
        double varianceQ_Sn= abs(Psq_Sn/N_Sn - pow(wavg_Sn,2)) ;
        double ErrQ_D  = sqrt(varianceQ_D);
        double ErrQ_Sn = sqrt(varianceQ_Sn);
        double unc = (wavg_Sn/wavg_D) *sqrt( pow( ErrQ_D/wavg_D,2) + pow(ErrQ_Sn/wavg_Sn ,2) );
        sin_Qter->SetPoint(binQ-1, x_Q_Sn, sratio_point );
        if (unc>5){
            unc = 3;
        }
        sin_Qter->SetPointError(binQ-1, 0, unc);
        //cout<<Errz_D<<"Errrrr....."<<endl;

    }

    for (int binv = 1; binv <= nubinsv; binv++) {
        double N_D = h_vtestD->GetBinContent(binv);
        double sin_D = h_vweightD->GetBinContent(binv);
        double N_Sn = h_vtestSn->GetBinContent(binv);        
        double sin_Sn = h_vweightSn->GetBinContent(binv);  
        double x_Q_Sn = h_vtestSn->GetXaxis()->GetBinCenter(binv);
        double Psq_D = h_vsqD->GetBinContent(binv); 
        double Psq_Sn = h_vsqSn->GetBinContent(binv);

        //<sin (phi_D)> = sin_D/N_D
        double wavg_D  = (N_D > 0) ? sin_D/N_D : 0.0;
        double wavg_Sn = (N_Sn > 0) ? sin_Sn/N_Sn : 0.0;
        double sratio_vpoint = (wavg_D > 0) ? wavg_Sn/wavg_D : 0.0; 
        //cout<<wavg_Sn<<"/" <<wavg_D <<" = " <<cratio_vpoint<<"   cratiopoint(v)"<< endl;
        double variancev_D = abs(Psq_D/N_D - pow(wavg_D,2));     
        double variancev_Sn= abs(Psq_Sn/N_Sn - pow(wavg_Sn,2)) ;
        double Errv_D  = sqrt(variancev_D);
        double Errv_Sn = sqrt(variancev_Sn);
        double unc = (wavg_Sn/wavg_D) *sqrt( pow( Errv_D/wavg_D,2) + pow(Errv_Sn/wavg_Sn ,2) );
        sin_vter->SetPoint(binv-1, x_Q_Sn, sratio_vpoint );
        sin_vter->SetPointError(binv-1, 0, unc);
        cout<<sin_Sn<<" sin_Sn .........-" <<endl;
    }

    

    for (int binz = 1; binz <= nubinsz; binz++) {
        double N_D = h_ztestD->GetBinContent(binz);
        double sin_D = h_zweightD->GetBinContent(binz);
        double N_Sn = h_ztestSn->GetBinContent(binz);        
        double sin_Sn = h_zweightSn->GetBinContent(binz);  
        double x_Q_Sn = h_ztestSn->GetXaxis()->GetBinCenter(binz);
        double Psq_D = h_vsqD->GetBinContent(binz);         
        double Psq_Sn = h_vsqSn->GetBinContent(binz);       

        //<sin (phi_D)> = sin_D/N_D
        double wavg_D  = (N_D > 0) ? sin_D/N_D : 0.0;
        double wavg_Sn = (N_Sn > 0) ? sin_Sn/N_Sn : 0.0;
        double sratio_zpoint = (wavg_D > 0) ? wavg_Sn/wavg_D : 0.0; 
        //cout<<wavg_Sn<<"/" <<wavg_D <<" = " <<cratio_zpoint<<"   cratiopoint(z)"<< endl;

        double variancez_D = abs(Psq_D/N_D - pow(wavg_D,2));     
        double variancez_Sn= abs(Psq_Sn/N_Sn - pow(wavg_Sn,2)) ;
        double Errz_D  = sqrt(variancez_D);
        double Errz_Sn = sqrt(variancez_Sn);
        double unc = (wavg_Sn/wavg_D) *sqrt( pow( Errz_D/wavg_D,2) + pow(Errz_Sn/wavg_Sn ,2) );
        sin_zter->SetPoint(binz-1, x_Q_Sn, sratio_zpoint );
        sin_zter->SetPointError(binz-1, 0, unc);


    }

    for (int binp = 1; binp <= nubinsp; binp++) {
        double N_D = h_pttestD->GetBinContent(binp);
        double sin_D = h_ptweightD->GetBinContent(binp);
        double N_Sn = h_pttestSn->GetBinContent(binp);        
        double sin_Sn = h_ptweightSn->GetBinContent(binp);  
        double x_Q_Sn = h_pttestSn->GetXaxis()->GetBinCenter(binp);
        double Psq_D = h_vsqD->GetBinContent(binp);         
        double Psq_Sn = h_vsqSn->GetBinContent(binp);       

        //<sin (phi_D)> = sin_D/N_D
        double wavg_D  = (N_D > 0) ? sin_D/N_D : 0.0;
        double wavg_Sn = (N_Sn > 0) ? sin_Sn/N_Sn : 0.0;
        double sratio_ppoint = (wavg_D > 0) ? wavg_Sn/wavg_D : 0.0; 
        //cout<<wavg_Sn<<"/" <<wavg_D <<" = " <<cratio_ppoint<<"   cratiopoint(p)"<< endl;
        double variancep_D = abs(Psq_D/N_D - pow(wavg_D,2));     
        double variancep_Sn= abs(Psq_Sn/N_Sn - pow(wavg_Sn,2));
        double Errp_D  = sqrt(variancep_D);
        double Errp_Sn = sqrt(variancep_Sn);
        double unc = (wavg_Sn/wavg_D) *sqrt( pow( Errp_D/wavg_D,2) + pow(Errp_Sn/wavg_Sn ,2) );
        sin_pter->SetPoint(binp-1, x_Q_Sn, sratio_ppoint );
        sin_pter->SetPointError(binp-1, 0, 0);

    }







    /////////fix weight test (end)///////////
    /////////////////////////////////////////

*/





    // Close the first ROOT file
    fileD->Close();
    fileSn->Close();
    

    
sR->cd(1);
sin_Q->GetXaxis()->SetTitleSize(0.05);
sin_Q->GetYaxis()->SetTitleSize(0.05);
sin_Q->SetMarkerSize(0.5);
sin_Q->SetMarkerStyle(21);
sin_Q->GetXaxis()->SetRangeUser(1,6.5);
sin_Q->GetXaxis()->SetTitle("Q^{2} " );
sin_Q->GetYaxis()->SetTitle("<sin(#Phi)>_{Sn} / <sin(#Phi)>_{D}");
sin_Q->Draw("APE");

sR->cd(2);
sin_v->GetXaxis()->SetTitleSize(0.05);
sin_v->GetYaxis()->SetTitleSize(0.05);
sin_v->GetXaxis()->SetRangeUser(2,9);
sin_v->SetMarkerSize(0.5);
sin_v->SetMarkerStyle(21);
sin_v->GetXaxis()->SetTitle("#nu " );
sin_v->GetYaxis()->SetTitle("<sin(#Phi)>_{Sn} / <sin(#Phi)>_{D}");
sin_v->Draw("APE");

sR->cd(3);
sin_z->GetXaxis()->SetTitleSize(0.05);
sin_z->GetYaxis()->SetTitleSize(0.05);
sin_z->SetMarkerSize(0.5);
sin_z->SetMarkerStyle(21);
sin_z->GetXaxis()->SetTitle("p_{t}^{2} " );
sin_z->GetYaxis()->SetTitle("<sin(#Phi)>_{Sn} / <sin(#Phi)>_{D}");
sin_z->Draw("APE");
sR->cd(4);
sin_pt->GetXaxis()->SetTitleSize(0.05);
sin_pt->GetYaxis()->SetTitleSize(0.05);
sin_pt->SetMarkerSize(0.5);
sin_pt->SetMarkerStyle(21);
sin_pt->GetXaxis()->SetRangeUser(0.2,0.9);
sin_pt->GetXaxis()->SetTitle("z " );
sin_pt->GetYaxis()->SetTitle("<sin(#Phi)>_{Sn} / <sin(#Phi)>_{D}");
sin_pt->Draw("APE");
    sR->SaveAs("testsRatio.pdf");
    sR->SaveAs("testsRatio.root");

    return 0;

}




