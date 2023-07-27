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
    TFile* fileD = new TFile("../files2read/output_D10binJul.root", "READ");
    TFile* fileSn = new TFile("../files2read/output_Sn10binJul.root", "READ");

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



    for (int binQ = 1; binQ <= numBinsQ; binQ++) {
        double sums = 0.0;
        double sumWeights = 0.0;
        double sums_Sn = 0.0;
        double sumWeights_Sn = 0.0;

        for (int binc = 1; binc <= numbinsS; binc++) {
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




    for (int binv = 1; binv <= numBinsv; binv++) {
        double sums = 0.0;
        double sumWeightsv = 0.0;
        double sums_Sn = 0.0;
        double sumWeightsv_Sn = 0.0;

        for (int binc = 1; binc <= numbinsS; binc++) {
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



    for (int binz = 1; binz <= numBinsz; binz++) {
        double sums = 0.0;
        double sumWeightsz = 0.0;
        double sums_Sn = 0.0;
        double sumWeightsz_Sn = 0.0;

        for (int binc = 1; binc <= numbinsS; binc++) {
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


    for (int binp = 1; binp <= numBinsp; binp++) {
        double sums = 0.0;
        double sumWeightsp = 0.0;
        double sums_Sn = 0.0;
        double sumWeightsp_Sn = 0.0;

        for (int binc = 1; binc <= numbinsS; binc++) {
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







for (int binQ = 1; binQ <= numBinsQ; binQ++) {
        double sinsumS = 0.0;
        double sumWeightsQ = 0.0;
        double sinsumS_Sn = 0.0;
        double sumWeightsQ_Sn = 0.0;
        double sinsumS_squared = 0.0;
        double sumWeightsQ_squared = 0.0;
        double sinsumS_Sn_squared = 0.0;
        double sumWeightsQ_Sn_squared = 0.0;
        double x_Q_Sn = h_QtestSn->GetXaxis()->GetBinCenter(binQ);
        for (int binS = 1; binS <= numbinsS; binS++) {
            double sinphi = h2DQ_D->GetXaxis()->GetBinCenter(binS);
            double weightQ = h2DQ_D->GetBinContent(binS, binQ);
            double sinphi_Sn = h2DQ_Sn->GetXaxis()->GetBinCenter(binQ);
            double weightQ_Sn = h2DQ_Sn->GetBinContent(binS, binQ);
            sinsumS += sinphi * weightQ;
            sumWeightsQ += weightQ;
            sinsumS_Sn += sinphi_Sn * weightQ_Sn;
            sumWeightsQ_Sn += weightQ_Sn;
            sinsumS_squared += sinphi * sinphi * weightQ;
            sumWeightsQ_squared += weightQ * weightQ;
            sinsumS_Sn_squared += sinphi_Sn * sinphi_Sn * weightQ_Sn;
            sumWeightsQ_Sn_squared += weightQ_Sn * weightQ_Sn;
        }
        double averagePC = sinsumS / sumWeightsQ;
        double averagePC_Sn = sinsumS_Sn / sumWeightsQ_Sn;
        double varPC = (sinsumS_squared / sumWeightsQ) - (averagePC * averagePC / (sumWeightsQ));
        double varPC_Sn = (sinsumS_Sn_squared / sumWeightsQ_Sn) - (averagePC_Sn * averagePC_Sn / (sumWeightsQ_Sn) );
  //      lines 2up and 2down added for error TBC  
        double uncertaintyPC = sqrt(varPC /sumWeightsQ);
        double uncertaintyPC_Sn = sqrt(varPC_Sn /sumWeightsQ_Sn);
        if (sumWeightsQ != 0 && sumWeightsQ_Sn != 0){
            double sin_point = averagePC_Sn / averagePC  ;
            double unc= (averagePC_Sn/averagePC) *sqrt( pow((sqrt(uncertaintyPC)/averagePC),2) + pow((sqrt(uncertaintyPC_Sn)/averagePC_Sn),2)    );
            sin_Qbis->SetPoint(binQ-1, x_Q_Sn, sin_point );
            sin_Qbis->SetPointError(binQ-1, unc );

        }

    }

    //int numbinsS = h2Dv_Sn->GetNbinsX();        //C is for sin; for sin(phi)

    for (int binv = 1; binv <= numBinsv; binv++) {
        double sinsumS = 0.0;
        double sumWeightsv = 0.0;
        double sinsumS_Sn = 0.0;
        double sumWeightsv_Sn = 0.0;
        double sinsumS_squared = 0.0;
        double sumWeightsv_squared = 0.0;
        double sinsumS_Sn_squared = 0.0;
        double sumWeightsv_Sn_squared = 0.0;
        double x_v_Sn = h_vtestSn->GetXaxis()->GetBinCenter(binv);
        for (int binS = 1; binS <= numbinsS; binS++) {
            double sinphi = h2Dv_D->GetXaxis()->GetBinCenter(binS);
            double weightv = h2Dv_D->GetBinContent(binS, binv);
            double sinphi_Sn = h2Dv_Sn->GetXaxis()->GetBinCenter(binv);
            double weightv_Sn = h2Dv_Sn->GetBinContent(binS, binv);
            sinsumS += sinphi * weightv;
            sumWeightsv += weightv;
            sinsumS_Sn += sinphi_Sn * weightv_Sn;
            sumWeightsv_Sn += weightv_Sn;
            sinsumS_squared += sinphi * sinphi * weightv;
            sumWeightsv_squared += weightv * weightv;
            sinsumS_Sn_squared += sinphi_Sn * sinphi_Sn * weightv_Sn;
            sumWeightsv_Sn_squared += weightv_Sn * weightv_Sn;
        }
        double averagePC = sinsumS / sumWeightsv;
        double averagePC_Sn = sinsumS_Sn / sumWeightsv_Sn;
        double varPC = (sinsumS_squared / sumWeightsv) - (averagePC * averagePC / (sumWeightsv));
        double varPC_Sn = (sinsumS_Sn_squared / sumWeightsv_Sn) - (averagePC_Sn * averagePC_Sn / (sumWeightsv_Sn) );
  //      lines 2up and 2down added for error TBC  
        double uncertaintyPC = sqrt(varPC /sumWeightsv);
        double uncertaintyPC_Sn = sqrt(varPC_Sn /sumWeightsv_Sn);
        if (sumWeightsv != 0 && sumWeightsv_Sn != 0){
            double sin_point = averagePC_Sn / averagePC  ;
            double unc= (averagePC_Sn/averagePC) *sqrt( pow((sqrt(uncertaintyPC)/averagePC),2) + pow((sqrt(uncertaintyPC_Sn)/averagePC_Sn),2)    );
            sin_vbis->SetPoint(binv-1, x_v_Sn, sin_point );
            sin_vbis->SetPointError(binv-1, unc );

        }

    }


     for (int binz = 1; binz <= numBinsz; binz++) {
        double sinsumS = 0.0;
        double sumWeightsz = 0.0;
        double sinsumS_Sn = 0.0;
        double sumWeightsz_Sn = 0.0;
        double sinsumS_squared = 0.0;
        double sumWeightsz_squared = 0.0;
        double sinsumS_Sn_squared = 0.0;
        double sumWeightsz_Sn_squared = 0.0;
        double x_z_Sn = h_ztestSn->GetXaxis()->GetBinCenter(binz);
        for (int binS = 1; binS <= numbinsS; binS++) {
            double sinphi = h2Dz_D->GetXaxis()->GetBinCenter(binS);
            double weightz = h2Dz_D->GetBinContent(binS, binz);
            double sinphi_Sn = h2Dz_Sn->GetXaxis()->GetBinCenter(binz);
            double weightz_Sn = h2Dz_Sn->GetBinContent(binS, binz);
            sinsumS += sinphi * weightz;
            sumWeightsz += weightz;
            sinsumS_Sn += sinphi_Sn * weightz_Sn;
            sumWeightsz_Sn += weightz_Sn;
            sinsumS_squared += sinphi * sinphi * weightz;
            sumWeightsz_squared += weightz * weightz;
            sinsumS_Sn_squared += sinphi_Sn * sinphi_Sn * weightz_Sn;
            sumWeightsz_Sn_squared += weightz_Sn * weightz_Sn;
        }
        double averagePC = sinsumS / sumWeightsz;
        double averagePC_Sn = sinsumS_Sn / sumWeightsz_Sn;
        double varPC = (sinsumS_squared / sumWeightsz) - (averagePC * averagePC / (sumWeightsz));
        double varPC_Sn = (sinsumS_Sn_squared / sumWeightsz_Sn) - (averagePC_Sn * averagePC_Sn / (sumWeightsz_Sn) );
  //      lines 2up and 2down added for error TBC  
        double uncertaintyPC = sqrt(varPC /sumWeightsz);
        double uncertaintyPC_Sn = sqrt(varPC_Sn /sumWeightsz_Sn);
        if (sumWeightsz != 0 && sumWeightsz_Sn != 0){
            double sin_point = averagePC_Sn / averagePC  ;
            double unc= (averagePC_Sn/averagePC) *sqrt( pow((sqrt(uncertaintyPC)/averagePC),2) + pow((sqrt(uncertaintyPC_Sn)/averagePC_Sn),2)    );
            sin_zbis->SetPoint(binz-1, x_z_Sn, sin_point );
            sin_zbis->SetPointError(binz-1, unc );

        }

    }


 for (int binp = 1; binp <= numBinsp; binp++) {
        double sinsumS = 0.0;
        double sumWeightsp = 0.0;
        double sinsumS_Sn = 0.0;
        double sumWeightsp_Sn = 0.0;
        double sinsumS_squared = 0.0;
        double sumWeightsp_squared = 0.0;
        double sinsumS_Sn_squared = 0.0;
        double sumWeightsp_Sn_squared = 0.0;
        double x_p_Sn = h_pttestSn->GetXaxis()->GetBinCenter(binp);
        for (int binS = 1; binS <= numbinsS; binS++) {
            double sinphi = h2Dp_D->GetXaxis()->GetBinCenter(binS);
            double weightp = h2Dp_D->GetBinContent(binS, binp);
            double sinphi_Sn = h2Dp_Sn->GetXaxis()->GetBinCenter(binp);
            double weightp_Sn = h2Dp_Sn->GetBinContent(binS, binp);
            sinsumS += sinphi * weightp;
            sumWeightsp += weightp;
            sinsumS_Sn += sinphi_Sn * weightp_Sn;
            sumWeightsp_Sn += weightp_Sn;
            sinsumS_squared += sinphi * sinphi * weightp;
            sumWeightsp_squared += weightp * weightp;
            sinsumS_Sn_squared += sinphi_Sn * sinphi_Sn * weightp_Sn;
            sumWeightsp_Sn_squared += weightp_Sn * weightp_Sn;
        }
        double averagePC = sinsumS / sumWeightsp;
        double averagePC_Sn = sinsumS_Sn / sumWeightsp_Sn;
        double varPC = (sinsumS_squared / sumWeightsp) - (averagePC * averagePC / (sumWeightsp));
        double varPC_Sn = (sinsumS_Sn_squared / sumWeightsp_Sn) - (averagePC_Sn * averagePC_Sn / (sumWeightsp_Sn) );
  //      lines 2up and 2down added for error TBC  
        double uncertaintyPC = sqrt(varPC /sumWeightsp);
        double uncertaintyPC_Sn = sqrt(varPC_Sn /sumWeightsp_Sn);
        if (sumWeightsp != 0 && sumWeightsp_Sn != 0){
            double sin_point = averagePC_Sn / averagePC  ;
            double unc= (sin_point) *sqrt( pow((sqrt(uncertaintyPC)/averagePC),2) + pow((sqrt(uncertaintyPC_Sn)/averagePC_Sn),2)    );
            sin_pbis->SetPoint(binp-1, x_p_Sn, sin_point );
            sin_pbis->SetPointError(binp-1, unc );
            cout<<unc<<"..."<<endl;
        }

    }








    // Close the first ROOT file
    fileD->Close();
    fileSn->Close();
    
sR->cd(1);
sin_Qter->GetXaxis()->SetTitleSize(0.05);
sin_Qter->GetYaxis()->SetTitleSize(0.05);
sin_Qter->SetMarkerSize(0.5);
sin_Qter->SetMarkerStyle(21);
sin_Qter->GetXaxis()->SetRangeUser(1,6.5);
sin_Qter->GetXaxis()->SetTitle("Q^{2} " );
sin_Qter->GetYaxis()->SetTitle("<sin(#Phi)>_{Sn} / <sin(#Phi)>_{D}");
sin_Qter->Draw("APE");

sR->cd(2);
sin_vter->GetXaxis()->SetTitleSize(0.05);
sin_vter->GetYaxis()->SetTitleSize(0.05);
sin_vter->GetXaxis()->SetRangeUser(2,9);
sin_vter->SetMarkerSize(0.5);
sin_vter->SetMarkerStyle(21);
sin_vter->GetXaxis()->SetTitle("#nu " );
sin_vter->GetYaxis()->SetTitle("<sin(#Phi)>_{Sn} / <sin(#Phi)>_{D}");
sin_vter->Draw("APE");

sR->cd(3);
sin_zter->GetXaxis()->SetTitleSize(0.05);
sin_zter->GetYaxis()->SetTitleSize(0.05);
sin_zter->SetMarkerSize(0.5);
sin_zter->SetMarkerStyle(21);
sin_zter->GetXaxis()->SetTitle("p_{t}^{2} " );
sin_zter->GetYaxis()->SetTitle("<sin(#Phi)>_{Sn} / <sin(#Phi)>_{D}");
sin_zter->Draw("APE");
sR->cd(4);
sin_pter->GetXaxis()->SetTitleSize(0.05);
sin_pter->GetYaxis()->SetTitleSize(0.05);
sin_pter->SetMarkerSize(0.5);
sin_pter->SetMarkerStyle(21);
sin_pter->GetXaxis()->SetRangeUser(0.2,0.9);
sin_pter->GetXaxis()->SetTitle("z " );
sin_pter->GetYaxis()->SetTitle("<sin(#Phi)>_{Sn} / <sin(#Phi)>_{D}");
sin_pter->Draw("APE");

    sR->SaveAs("testsRatio.pdf");
    sR->SaveAs("testsRatio.root");

    return 0;

}




