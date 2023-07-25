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
    TCanvas* cDpt3=new TCanvas("Dptplot3","Dptplot3"); //create new canvas
    TCanvas* ctest=new TCanvas("test","test"); //create new canvas
    cDpt->Divide(2,2);
    TGraphErrors *Dpt_Q = new TGraphErrors();   //creation of emmpty graphs to be filled at the end 
    TGraphErrors *Dpt_Qbis = new TGraphErrors();
    TGraphErrors *Dpt_Qsec = new TGraphErrors();
    TGraphErrors *Dpt_v = new TGraphErrors();
    TGraphErrors *Dpt_vbis = new TGraphErrors();
    TGraphErrors *Dpt_vsec = new TGraphErrors();    

    TGraphErrors *Dpt_z = new TGraphErrors();
    TGraphErrors *Dpt_zbis = new TGraphErrors();
    TGraphErrors *Dpt_zsec = new TGraphErrors();

    TGraphErrors *Dpt_Qter = new TGraphErrors();

    // Open the first ROOT file
    TFile* fileD = new TFile("../files2read/testOutput_D.root", "READ");
    TFile* fileSn = new TFile("../files2read/testOutput_Sn.root", "READ");
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

    //TH2F* h_2DtestD  = dynamic_cast<TH2F*>(fileD->Get("pt2Q2gen"));
    //TH2F* h_2DtestSn = dynamic_cast<TH2F*>(fileSn->Get("pt2Q2gen"));
    //maybe not even need to dynamic cast, maybe hjust enough with calling the class
    //-> Yes, dynamic cast is needed, to make the difference between Sn and De
    //-> an adjustment has to be done to sumW2 histograms in order to properly recover them making the distinction on them. 
  
    
    int numBinsQ = h2DQ_Sn->GetNbinsY();
    int numBinsP = h2DQ_Sn->GetNbinsX();
    int numBinsv = h2Dv_D->GetNbinsY();
    int numBinsz = h2Dz_D->GetNbinsY();








//////test weighted average fix //////////

    for (int binQ = 1; binQ <= nubinsQ; binQ++) {
        //pt is already squared we call it P = p_t * p_t as stocked in initial histo in main.cpp
        double N_D = h_QtestD->GetBinContent(binQ);
        double P_D = h_QweightD->GetBinContent(binQ);
        double N_Sn = h_QtestSn->GetBinContent(binQ);        
        double P_Sn = h_QweightSn->GetBinContent(binQ);  
        double x_Q_Sn = h_QtestSn->GetXaxis()->GetBinCenter(binQ);


        //<cos (phi_D)> = cos_D/N_D
        double wavg_D  = (N_D > 0) ? P_D/N_D : 0.0;
        double wavg_Sn = (N_Sn > 0) ? P_Sn/N_Sn : 0.0;
        double cratio_point = wavg_Sn - wavg_D ; 
        cout<<wavg_Sn<<"/" <<wavg_D <<" = " <<cratio_point<<"   cratiopoint"<< endl;
        Dpt_Qter->SetPoint(binQ-1, x_Q_Sn, cratio_point );

    }





//////test weighted average fix //////////



    for (int binQ = 1; binQ <= numBinsQ; binQ++) {
        double sumP = 0.0;
        double countsP = 0.0;
        double countsP_Sn = 0.0;
        double sumWeights = 0.0;
        double sumP_Sn = 0.0;
        double sumWeights_Sn = 0.0;
        double x_q_Sn = h_QtestSn->GetXaxis()->GetBinCenter(binQ);

        for (int binP = 1; binP <= numBinsP; binP++) {
            double p = h2DQ_D->GetXaxis()->GetBinCenter(binP);
            double weight = h2DQ_D->GetBinContent(binP, binQ);
            double p_Sn = h2DQ_Sn->GetXaxis()->GetBinCenter(binP);
            double weight_Sn = h2DQ_Sn->GetBinContent(binP, binQ);
            countsP += p; 
            countsP_Sn += p_Sn; 
            sumP += p * weight;
            sumWeights += weight;
            sumP_Sn += p_Sn * weight_Sn;
            sumWeights_Sn += weight_Sn;
        }
        
        double averageP = sumP / sumWeights;
        //double unc_sumP = sumP * sqrt(  pow( (sqrt(countsP)/countsP)   ,2) + pow((sqrt(sumWeights)/sumWeights)  ,2)      ) // uncertainty up
        //double unc_w = sqrt(sumWeights); //unc down
        double averageP_Sn = sumP_Sn / sumWeights_Sn;
        double unc_z = averageP * (  sqrt( pow(sqrt(countsP_Sn)/countsP_Sn  ,2) + pow(sqrt(sumWeights)/sumWeights  ,2) ) + pow( sqrt(sumWeights)/sumWeights ,2)    ) ;
        double unc_zSn = averageP_Sn * (  sqrt( pow(sqrt(countsP)/countsP  ,2) + pow(sqrt(sumWeights_Sn)/sumWeights_Sn  ,2) ) + pow( sqrt(sumWeights_Sn)/sumWeights_Sn ,2)    ); 
        
        if (sumWeights != 0 && sumWeights_Sn != 0){
                double pt_point = averageP_Sn - averageP  ;
                double unc_tot = unc_zSn + unc_z ; 
                Dpt_Qsec->SetPoint(binQ-1, x_q_Sn, pt_point );
                Dpt_Qsec->SetPointError(binQ-1, 0,unc_tot );


        }
        //cout<<averageP_Sn<< " - " <<averageP<<" = "<<averageP_Sn - averageP<< endl;
    }




    for (int binQ = 1; binQ <= numBinsQ; binQ++) {
        double sumP = 0.0;
        double sumWeightsQ = 0.0;
        double sumP_Sn = 0.0;
        double sumWeightsQ_Sn = 0.0;
        double sumP_squared = 0.0;
        double sumWeightsQ_squared = 0.0;
        double sumP_Sn_squared = 0.0;
        double sumWeightsQ_Sn_squared = 0.0;
        double x_Q_Sn = h_QtestSn->GetXaxis()->GetBinCenter(binQ);

        for (int binP = 1; binP <= numBinsP; binP++) {
            double p = h2DQ_D->GetXaxis()->GetBinCenter(binP);
            double weightQ = h2DQ_D->GetBinContent(binP, binQ);
            double p_Sn = h2DQ_Sn->GetXaxis()->GetBinCenter(binQ);
            double weightQ_Sn = h2DQ_Sn->GetBinContent(binP, binQ);
            sumP += p * weightQ;
            sumWeightsQ += weightQ;
            sumP_Sn += p_Sn * weightQ_Sn;
            sumWeightsQ_Sn += weightQ_Sn;
            sumP_squared += p * p * weightQ;
            sumWeightsQ_squared += weightQ * weightQ;
            sumP_Sn_squared += p_Sn * p_Sn * weightQ_Sn;
            sumWeightsQ_Sn_squared += weightQ_Sn * weightQ_Sn;
        }

        double averagePQ = sumP / sumWeightsQ;
        double averagePQ_Sn = sumP_Sn / sumWeightsQ_Sn;
        //uncomment below maybe errors to be rechecked (corrected) TBC
        //double varPQ = (sumP_squared / sumWeightsQ) - (averagePQ * averagePQ);
        //double varPQ_Sn = (sumP_Sn_squared / sumWeightsQ_Sn) - (averagePQ_Sn * averagePQ_Sn);
        //double uncertaintyPQ = sqrt(varPQ / (sumWeightsQ - 1));
        //double uncertaintyPQ_Sn = sqrt(varPQ_Sn / (sumWeightsQ_Sn - 1));

        double varPQ = (sumP_squared / sumWeightsQ) - (averagePQ * averagePQ / (sumWeightsQ));
        double varPQ_Sn = (sumP_Sn_squared / sumWeightsQ_Sn) - (averagePQ_Sn * averagePQ_Sn / (sumWeightsQ_Sn) );
  //      lines 2up and 2down added for error TBC  
        double uncertaintyPQ = sqrt(varPQ /sumWeightsQ);
        double uncertaintyPQ_Sn = sqrt(varPQ_Sn /sumWeightsQ_Sn);



        //double averagePz = sumP / sumWeightsz;
        //double averagePz_Sn = sumP_Sn / sumWeightsz_Sn;
        if (sumWeightsQ != 0 && sumWeightsQ_Sn != 0){
            double pt_point = averagePQ_Sn - averagePQ  ;
            cout     << sumWeightsQ_Sn<<endl; 
            //cout     << averagePQ_Sn << " ± " << uncertaintyPQ_Sn << " - " << averagePQ << " ± " << uncertaintyPQ << " = " << averagePQ_Sn - averagePQ << " ± " << sqrt(uncertaintyPQ_Sn * uncertaintyPQ_Sn + uncertaintyPQ * uncertaintyPQ) << endl;
            double unc= sqrt(uncertaintyPQ_Sn * uncertaintyPQ_Sn + uncertaintyPQ * uncertaintyPQ);
            //cout<<averagePQ_Sn<< " - " <<averagePQ<<" = "<<averagePQ_Sn - averagePQ<< endl;
            Dpt_Qbis->SetPoint(binQ-1, x_Q_Sn, pt_point );
            Dpt_Qbis->SetPointError(binQ-1, 0,unc );

        }
    }


    





    for (int binv = 1; binv <= numBinsv; binv++) {
        double sumP = 0.0;
        double sumWeightsv = 0.0;
        double sumP_Sn = 0.0;
        double sumWeightsv_Sn = 0.0;
        double x_v_Sn = h_vtestSn->GetXaxis()->GetBinCenter(binv);


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
        if (sumWeightsv != 0 && sumWeightsv_Sn != 0){
                double pt_point = averagePv_Sn - averagePv  ;
                Dpt_v->SetPoint(binv-1, x_v_Sn, pt_point );

        }
    }

    for (int binv = 1; binv <= numBinsv; binv++) {
        double sumP = 0.0;
        double sumWeightsv = 0.0;
        double sumP_Sn = 0.0;
        double sumWeightsv_Sn = 0.0;
        double sumP_squared = 0.0;
        double sumWeightsv_squared = 0.0;
        double sumP_Sn_squared = 0.0;
        double sumWeightsv_Sn_squared = 0.0;
        double x_v_Sn = h_vtestSn->GetXaxis()->GetBinCenter(binv);

        for (int binP = 1; binP <= numBinsP; binP++) {
            double p = h2Dv_D->GetXaxis()->GetBinCenter(binP);
            double weightv = h2Dv_D->GetBinContent(binP, binv);
            double p_Sn = h2Dv_Sn->GetXaxis()->GetBinCenter(binv);
            double weightv_Sn = h2Dv_Sn->GetBinContent(binP, binv);
            sumP += p * weightv;
            sumWeightsv += weightv;
            sumP_Sn += p_Sn * weightv_Sn;
            sumWeightsv_Sn += weightv_Sn;
            sumP_squared += p * p * weightv;
            sumWeightsv_squared += weightv * weightv;
            sumP_Sn_squared += p_Sn * p_Sn * weightv_Sn;
            sumWeightsv_Sn_squared += weightv_Sn * weightv_Sn;
        }

        double averagePv = sumP / sumWeightsv;
        double averagePv_Sn = sumP_Sn / sumWeightsv_Sn;
        //uncomment this part
        //double varPv = (sumP_squared / sumWeightsv) - (averagePv * averagePv);
        //double varPv_Sn = (sumP_Sn_squared / sumWeightsv_Sn) - (averagePv_Sn * averagePv_Sn);
        //double uncertaintyPv = sqrt(varPv / (sumWeightsv - 1));
        //double uncertaintyPv_Sn = sqrt(varPv_Sn / (sumWeightsv_Sn - 1));

        double varPv = (sumP_squared / sumWeightsv) - (averagePv * averagePv / (sumWeightsv));
        double varPv_Sn = (sumP_Sn_squared / sumWeightsv_Sn) - (averagePv_Sn * averagePv_Sn / (sumWeightsv_Sn) );
  //      lines 2up and 2down added for error TBC  
        double uncertaintyPv = sqrt(varPv /sumWeightsv);
        double uncertaintyPv_Sn = sqrt(varPv_Sn /sumWeightsv_Sn);

        //double averagePz = sumP / sumWeightsz;
        //double averagePz_Sn = sumP_Sn / sumWeightsz_Sn;
        if (sumWeightsv != 0 && sumWeightsv_Sn != 0){
            double pt_point = averagePv_Sn - averagePv  ;
            //cout     << averagePv_Sn << " ± " << uncertaintyPv_Sn << " - " << averagePv << " ± " << uncertaintyPv << " = " << averagePv_Sn - averagePv << " ± " << sqrt(uncertaintyPv_Sn * uncertaintyPv_Sn + uncertaintyPv * uncertaintyPv) << endl;
            double unc= sqrt(uncertaintyPv_Sn * uncertaintyPv_Sn + uncertaintyPv * uncertaintyPv);
            //cout<<averagePv_Sn<< " - " <<averagePv<<" = "<<averagePv_Sn - averagePv<< endl;
            Dpt_vbis->SetPoint(binv-1, x_v_Sn, pt_point );
            Dpt_vbis->SetPointError(binv-1, 0,unc );

        }
    }











    for (int binz = 1; binz <= numBinsz; binz++) {
        double sumP = 0.0;
        double sumWeightsz = 0.0;
        double sumP_Sn = 0.0;
        double sumWeightsz_Sn = 0.0;
        double sumP_squared = 0.0;
        double sumWeightsz_squared = 0.0;
        double sumP_Sn_squared = 0.0;
        double sumWeightsz_Sn_squared = 0.0;
        double x_z_Sn = h_ztestSn->GetXaxis()->GetBinCenter(binz);

        for (int binP = 1; binP <= numBinsP; binP++) {
            double p = h2Dz_D->GetXaxis()->GetBinCenter(binP);
            double weightz = h2Dz_D->GetBinContent(binP, binz);
            double p_Sn = h2Dz_Sn->GetXaxis()->GetBinCenter(binz);
            double weightz_Sn = h2Dz_Sn->GetBinContent(binP, binz);
            sumP += p * weightz;
            sumWeightsz += weightz;
            sumP_Sn += p_Sn * weightz_Sn;
            sumWeightsz_Sn += weightz_Sn;
            sumP_squared += p * p * weightz;
            sumWeightsz_squared += weightz * weightz;
            sumP_Sn_squared += p_Sn * p_Sn * weightz_Sn;
            sumWeightsz_Sn_squared += weightz_Sn * weightz_Sn;
        }

        double averagePz = sumP / sumWeightsz;
        double averagePz_Sn = sumP_Sn / sumWeightsz_Sn;

        //double varPz = (sumP_squared / sumWeightsz) - (averagePz * averagePz);
        //double varPz_Sn = (sumP_Sn_squared / sumWeightsz_Sn) - (averagePz_Sn * averagePz_Sn);
//
        //double uncertaintyPz = sqrt(varPz / (sumWeightsz-1 ));
        //double uncertaintyPz_Sn = sqrt(varPz_Sn / (sumWeightsz_Sn -1 ));
 
        double varPz = (sumP_squared / sumWeightsz) - (averagePz * averagePz / (sumWeightsz));
        double varPz_Sn = (sumP_Sn_squared / sumWeightsz_Sn) - (averagePz_Sn * averagePz_Sn / (sumWeightsz_Sn) );
  //      
        double uncertaintyPz = sqrt(varPz /sumWeightsz);
        double uncertaintyPz_Sn = sqrt(varPz_Sn /sumWeightsz_Sn);
        
        //double uncertaintyPz = sqrt(varPz );
        //double uncertaintyPz_Sn = sqrt(varPz_Sn );

        //double averagePz = sumP / sumWeightsz;
        //double averagePz_Sn = sumP_Sn / sumWeightsz_Sn;
        if (sumWeightsz != 0 && sumWeightsz_Sn != 0){
            double pt_point = averagePz_Sn - averagePz  ;
            //cout     << averagePz_Sn << " ± " << uncertaintyPz_Sn << " - " << averagePz << " ± " << uncertaintyPz << " = " << averagePz_Sn - averagePz << " ± " << sqrt(uncertaintyPz_Sn * uncertaintyPz_Sn + uncertaintyPz * uncertaintyPz) << endl;
            double unc= sqrt(uncertaintyPz_Sn * uncertaintyPz_Sn + uncertaintyPz * uncertaintyPz);
            //cout<<averagePz_Sn<< " - " <<averagePz<<" = "<<averagePz_Sn - averagePz<< endl;
            Dpt_zbis->SetPoint(binz-1, x_z_Sn, pt_point );
            Dpt_zbis->SetPointError(binz-1, 0,unc );

        }
    }



    for (int binz = 1; binz <= numBinsz; binz++) {
        double sumP = 0.0;
        double sumWeightsz = 0.0;
        double sumP_Sn = 0.0;
        double sumWeightsz_Sn = 0.0;
        double x_z_Sn = h_ztestSn->GetXaxis()->GetBinCenter(binz);

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
        if (sumWeightsz != 0 && sumWeightsz_Sn != 0){
                double pt_point = averagePz_Sn - averagePz  ;
                //cout<<averagePz_Sn<< " - " <<averagePz<<" = "<<averagePz_Sn - averagePz<< endl;
                Dpt_z->SetPoint(binz-1, x_z_Sn, pt_point );


        }


    }

    // Close the first ROOT file
    fileD->Close();
    fileSn->Close();
    cDpt->cd(1);
    //Dpt_Q->SetTitle(title);
    Dpt_Qbis->GetXaxis()->SetTitle("Q{2} (GeV^{2}) " );
    Dpt_Qbis->GetYaxis()->SetTitle("P_{t}^{2}diff) ");
    Dpt_Qbis->Draw("APE");
    cDpt->cd(2);
    //Dpt_v->SetTitle(title);
    Dpt_vbis->GetXaxis()->SetTitle("v " );
    Dpt_vbis->GetYaxis()->SetTitle("P_{t}^{2}diff) ");
    Dpt_vbis->Draw("APE");
    cDpt->cd(3);
    //Dpt_z->SetTitle(title);
    Dpt_zbis->GetXaxis()->SetTitle("Q  " );
    Dpt_zbis->GetYaxis()->SetTitle("P_{t}^{2}diff) ");
    Dpt_zbis->Draw("APE");
    //h_QtestD->Draw("COLZ");
    cDpt->cd(4);
    Dpt_zbis->GetXaxis()->SetTitle("z  " );
    Dpt_zbis->GetYaxis()->SetTitle("P_{t}^{2}diff) ");
    Dpt_zbis->Draw("APE");
    cDpt->SaveAs("dpttest.pdf");
    cDpt->SaveAs("dpttest.root");


    cDpt3->cd(1);
    //Dpt_Q->SetTitle(title);
    Dpt_Qsec->GetXaxis()->SetTitle("Q{2} (GeV^{2}) " );
    Dpt_Qsec->GetYaxis()->SetTitle("P_{t}^{2}diff) ");
    Dpt_Qsec->Draw("APE");
    cDpt3->cd(2);
    //Dpt_v->SetTitle(title);
    Dpt_vsec->GetXaxis()->SetTitle("v " );
    Dpt_vsec->GetYaxis()->SetTitle("P_{t}^{2}diff) ");
    Dpt_vsec->Draw("APE");
    cDpt3->cd(3);
    //Dpt_z->SetTitle(title);
    Dpt_zsec->GetXaxis()->SetTitle("z  " );
    Dpt_zsec->GetYaxis()->SetTitle("P_{t}^{2}diff) ");
    Dpt_zsec->Draw("APE");
    //h_QtestD->Draw("COLZ");
    cDpt3->cd(4);
    Dpt_zsec->GetXaxis()->SetTitle("z  " );
    Dpt_zsec->GetYaxis()->SetTitle("P_{t}^{2}diff) ");
    Dpt_zsec->Draw("APE");
    cDpt3->SaveAs("dptsec.pdf");
    cDpt3->SaveAs("dptsec.root");

    return 0;
}

