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
    TGraphErrors *cos_Qbis= new TGraphErrors();
    TGraphErrors *cos_vbis= new TGraphErrors();
    TGraphErrors *cos_zbis= new TGraphErrors();
    TGraphErrors *cos_pbis = new TGraphErrors();
    TGraphErrors *cos_Qter = new TGraphErrors();
    // Open the first ROOT file
    //TFile* file1 = new TFile("build/output1.root", "READ");
    //TFile* fileD = new TFile("../files2read/testOutput_D.root", "READ");
    //TFile* fileSn = new TFile("../files2read/testOutput_Sn.root", "READ");
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
    
    TH1F* h_QweightD = dynamic_cast<TH1F*>(fileD->Get("weightedQcos"));
       
    TH1F* h_QweightSn = dynamic_cast<TH1F*>(fileSn->Get("weightedQcos"));
    TH1F* h_vweightD = dynamic_cast<TH1F*>(fileD->Get("weightedvcos"));
    TH1F* h_vweightSn = dynamic_cast<TH1F*>(fileSn->Get("weightedvcos"));
    TH1F* h_zweightD = dynamic_cast<TH1F*>(fileD->Get("weightedzcos"));
    TH1F* h_zweightSn = dynamic_cast<TH1F*>(fileSn->Get("weightedzcos"));
    TH1F* h_ptweightD = dynamic_cast<TH1F*>(fileD->Get("weightedP_tcos"));
    TH1F* h_ptweightSn = dynamic_cast<TH1F*>(fileSn->Get("weightedP_tcos"));

//  ward to see if the histogram is existent
/*
    if (!h_QweightD) {
        std::cerr << "Error: Histogram not found in the file." << std::endl;
        fileD->Close();
        return 0;
    }
*/
//


    int nubinsQ = h_QtestD->GetNbinsX(); 
    int numBinsQ = h2DQ_D->GetNbinsY();
    int numBinsc = h2DQ_D->GetNbinsX();
    int numBinsv = h2Dv_D->GetNbinsY();
    int numBinsz = h2Dz_D->GetNbinsY();
    int numBinsp = h2Dp_D->GetNbinsY();
    int numBinsC = h2DQ_Sn->GetNbinsX();        //C is for cos; for cos(phi)


    ///////////weighted TH1F test/////////////////
    for (int binQ = 1; binQ <= nubinsQ; binQ++) {
        double N_D = h_QtestD->GetBinContent(binQ);
        double cos_D = h_QweightD->GetBinContent(binQ);
        double N_Sn = h_QtestSn->GetBinContent(binQ);        
        double cos_Sn = h_QweightSn->GetBinContent(binQ);  
        double x_Q_Sn = h_QtestSn->GetXaxis()->GetBinCenter(binQ);


        //<cos (phi_D)> = cos_D/N_D
        double wavg_D  = (N_D > 0) ? cos_D/N_D : 0.0;
        double wavg_Sn = (N_Sn > 0) ? cos_Sn/N_Sn : 0.0;
        double cratio_point = (wavg_D > 0) ? wavg_Sn/wavg_D : 0.0; 
        cout<<wavg_Sn<<"/" <<wavg_D <<" = " <<cratio_point<<"   cratiopoint"<< endl;
        cos_Qter->SetPoint(binQ-1, x_Q_Sn, cratio_point );





    }

    ///////////weighted TH1F test/////////////////




















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





    //int numBinsQ = h2DQ_Sn->GetNbinsY();

    for (int binQ = 1; binQ <= numBinsQ; binQ++) {
        double sumC = 0.0;
        double sumWeightsQ = 0.0;
        double sumC_Sn = 0.0;
        double sumWeightsQ_Sn = 0.0;
        double sumC_squared = 0.0;
        double sumWeightsQ_squared = 0.0;
        double sumC_Sn_squared = 0.0;
        double sumWeightsQ_Sn_squared = 0.0;
        double x_Q_Sn = h_QtestSn->GetXaxis()->GetBinCenter(binQ);
        for (int binC = 1; binC <= numBinsC; binC++) {
            double cosphi = h2DQ_D->GetXaxis()->GetBinCenter(binC);
            double weightQ = h2DQ_D->GetBinContent(binC, binQ);
            double cosphi_Sn = h2DQ_Sn->GetXaxis()->GetBinCenter(binQ);
            double weightQ_Sn = h2DQ_Sn->GetBinContent(binC, binQ);
            sumC += cosphi * weightQ;
            sumWeightsQ += weightQ;
            sumC_Sn += cosphi_Sn * weightQ_Sn;
            sumWeightsQ_Sn += weightQ_Sn;
            sumC_squared += cosphi * cosphi * weightQ;
            sumWeightsQ_squared += weightQ * weightQ;
            sumC_Sn_squared += cosphi_Sn * cosphi_Sn * weightQ_Sn;
            sumWeightsQ_Sn_squared += weightQ_Sn * weightQ_Sn;
        }
        double averagePC = sumC / sumWeightsQ;
        double averagePC_Sn = sumC_Sn / sumWeightsQ_Sn;
        double varPC = (sumC_squared / sumWeightsQ) - (averagePC * averagePC / (sumWeightsQ));
        double varPC_Sn = (sumC_Sn_squared / sumWeightsQ_Sn) - (averagePC_Sn * averagePC_Sn / (sumWeightsQ_Sn) );
  //      lines 2up and 2down added for error TBC  
        double uncertaintyPC = sqrt(varPC /sumWeightsQ);
        double uncertaintyPC_Sn = sqrt(varPC_Sn /sumWeightsQ_Sn);
        if (sumWeightsQ != 0 && sumWeightsQ_Sn != 0){
            double cos_point = averagePC_Sn / averagePC  ;
            double unc= (averagePC_Sn/averagePC) *sqrt( pow((sqrt(uncertaintyPC)/averagePC),2) + pow((sqrt(uncertaintyPC_Sn)/averagePC_Sn),2)    );
            cos_Qbis->SetPoint(binQ-1, x_Q_Sn, cos_point );
            cos_Qbis->SetPointError(binQ-1, 0,0.5 );

        }

    }

    //int numBinsC = h2Dv_Sn->GetNbinsX();        //C is for cos; for cos(phi)

    for (int binv = 1; binv <= numBinsv; binv++) {
        double sumC = 0.0;
        double sumWeightsv = 0.0;
        double sumC_Sn = 0.0;
        double sumWeightsv_Sn = 0.0;
        double sumC_squared = 0.0;
        double sumWeightsv_squared = 0.0;
        double sumC_Sn_squared = 0.0;
        double sumWeightsv_Sn_squared = 0.0;
        double x_v_Sn = h_vtestSn->GetXaxis()->GetBinCenter(binv);
        for (int binC = 1; binC <= numBinsC; binC++) {
            double cosphi = h2Dv_D->GetXaxis()->GetBinCenter(binC);
            double weightv = h2Dv_D->GetBinContent(binC, binv);
            double cosphi_Sn = h2Dv_Sn->GetXaxis()->GetBinCenter(binv);
            double weightv_Sn = h2Dv_Sn->GetBinContent(binC, binv);
            sumC += cosphi * weightv;
            sumWeightsv += weightv;
            sumC_Sn += cosphi_Sn * weightv_Sn;
            sumWeightsv_Sn += weightv_Sn;
            sumC_squared += cosphi * cosphi * weightv;
            sumWeightsv_squared += weightv * weightv;
            sumC_Sn_squared += cosphi_Sn * cosphi_Sn * weightv_Sn;
            sumWeightsv_Sn_squared += weightv_Sn * weightv_Sn;
        }
        double averagePC = sumC / sumWeightsv;
        double averagePC_Sn = sumC_Sn / sumWeightsv_Sn;
        double varPC = (sumC_squared / sumWeightsv) - (averagePC * averagePC / (sumWeightsv));
        double varPC_Sn = (sumC_Sn_squared / sumWeightsv_Sn) - (averagePC_Sn * averagePC_Sn / (sumWeightsv_Sn) );
  //      lines 2up and 2down added for error TBC  
        double uncertaintyPC = sqrt(varPC /sumWeightsv);
        double uncertaintyPC_Sn = sqrt(varPC_Sn /sumWeightsv_Sn);
        if (sumWeightsv != 0 && sumWeightsv_Sn != 0){
            double cos_point = averagePC_Sn / averagePC  ;
            double unc= (averagePC_Sn/averagePC) *sqrt( pow((sqrt(uncertaintyPC)/averagePC),2) + pow((sqrt(uncertaintyPC_Sn)/averagePC_Sn),2)    );
            cos_vbis->SetPoint(binv-1, x_v_Sn, cos_point );
            cos_vbis->SetPointError(binv-1, 0,0.5 );

        }

    }


     for (int binz = 1; binz <= numBinsz; binz++) {
        double sumC = 0.0;
        double sumWeightsz = 0.0;
        double sumC_Sn = 0.0;
        double sumWeightsz_Sn = 0.0;
        double sumC_squared = 0.0;
        double sumWeightsz_squared = 0.0;
        double sumC_Sn_squared = 0.0;
        double sumWeightsz_Sn_squared = 0.0;
        double x_z_Sn = h_ztestSn->GetXaxis()->GetBinCenter(binz);
        for (int binC = 1; binC <= numBinsC; binC++) {
            double cosphi = h2Dz_D->GetXaxis()->GetBinCenter(binC);
            double weightz = h2Dz_D->GetBinContent(binC, binz);
            double cosphi_Sn = h2Dz_Sn->GetXaxis()->GetBinCenter(binz);
            double weightz_Sn = h2Dz_Sn->GetBinContent(binC, binz);
            sumC += cosphi * weightz;
            sumWeightsz += weightz;
            sumC_Sn += cosphi_Sn * weightz_Sn;
            sumWeightsz_Sn += weightz_Sn;
            sumC_squared += cosphi * cosphi * weightz;
            sumWeightsz_squared += weightz * weightz;
            sumC_Sn_squared += cosphi_Sn * cosphi_Sn * weightz_Sn;
            sumWeightsz_Sn_squared += weightz_Sn * weightz_Sn;
        }
        double averagePC = sumC / sumWeightsz;
        double averagePC_Sn = sumC_Sn / sumWeightsz_Sn;
        double varPC = (sumC_squared / sumWeightsz) - (averagePC * averagePC / (sumWeightsz));
        double varPC_Sn = (sumC_Sn_squared / sumWeightsz_Sn) - (averagePC_Sn * averagePC_Sn / (sumWeightsz_Sn) );
  //      lines 2up and 2down added for error TBC  
        double uncertaintyPC = sqrt(varPC /sumWeightsz);
        double uncertaintyPC_Sn = sqrt(varPC_Sn /sumWeightsz_Sn);
        if (sumWeightsz != 0 && sumWeightsz_Sn != 0){
            double cos_point = averagePC_Sn / averagePC  ;
            double unc= (averagePC_Sn/averagePC) *sqrt( pow((sqrt(uncertaintyPC)/averagePC),2) + pow((sqrt(uncertaintyPC_Sn)/averagePC_Sn),2)    );
            cos_zbis->SetPoint(binz-1, x_z_Sn, cos_point );
            cos_zbis->SetPointError(binz-1, 0,0.5 );

        }

    }


 for (int binp = 1; binp <= numBinsp; binp++) {
        double sumC = 0.0;
        double sumWeightsp = 0.0;
        double sumC_Sn = 0.0;
        double sumWeightsp_Sn = 0.0;
        double sumC_squared = 0.0;
        double sumWeightsp_squared = 0.0;
        double sumC_Sn_squared = 0.0;
        double sumWeightsp_Sn_squared = 0.0;
        double x_p_Sn = h_pttestSn->GetXaxis()->GetBinCenter(binp);
        for (int binC = 1; binC <= numBinsC; binC++) {
            double cosphi = h2Dp_D->GetXaxis()->GetBinCenter(binC);
            double weightp = h2Dp_D->GetBinContent(binC, binp);
            double cosphi_Sn = h2Dp_Sn->GetXaxis()->GetBinCenter(binp);
            double weightp_Sn = h2Dp_Sn->GetBinContent(binC, binp);
            sumC += cosphi * weightp;
            sumWeightsp += weightp;
            sumC_Sn += cosphi_Sn * weightp_Sn;
            sumWeightsp_Sn += weightp_Sn;
            sumC_squared += cosphi * cosphi * weightp;
            sumWeightsp_squared += weightp * weightp;
            sumC_Sn_squared += cosphi_Sn * cosphi_Sn * weightp_Sn;
            sumWeightsp_Sn_squared += weightp_Sn * weightp_Sn;
        }
        double averagePC = sumC / sumWeightsp;
        double averagePC_Sn = sumC_Sn / sumWeightsp_Sn;
        double varPC = (sumC_squared / sumWeightsp) - (averagePC * averagePC / (sumWeightsp));
        double varPC_Sn = (sumC_Sn_squared / sumWeightsp_Sn) - (averagePC_Sn * averagePC_Sn / (sumWeightsp_Sn) );
  //      lines 2up and 2down added for error TBC  
        double uncertaintyPC = sqrt(varPC /sumWeightsp);
        double uncertaintyPC_Sn = sqrt(varPC_Sn /sumWeightsp_Sn);
        if (sumWeightsp != 0 && sumWeightsp_Sn != 0){
            double cos_point = averagePC_Sn / averagePC  ;
            double unc= (cos_point) *sqrt( pow((sqrt(uncertaintyPC)/averagePC),2) + pow((sqrt(uncertaintyPC_Sn)/averagePC_Sn),2)    );
            cos_pbis->SetPoint(binp-1, x_p_Sn, cos_point );
            cos_pbis->SetPointError(binp-1, 0,0.5 );

        }

    }





////////////////trying cratio w/ clascoll methods (?)////////////////////////////

//----cratio @ v try2---//		//"WORKING" in CLAS 


//h2DQ_D ; h2DQ_Sn
/*
double* avgc3 = new double[numBinsQ];
double* errc3 = new double[numBinsQ];
double* avgc4 = new double[nbinsxc4];
double* errc4 = new double[nbinsxc4];

for (int xbinc3 = 1; xbinc3 <= numBinsc; xbinc3++) {
  double sumc3 = 0.0;
  double sumsqc3 = 0.0;
  int countc3 = 0;
  double sumc4 = 0.0;
  double sumsqc4 = 0.0;
  int countc4 = 0;
  for (int ybinc3 = 1; ybinc3 <= numBinsQ; ybinc3++) {
    double contentc3 = h2DQ_Sn->GetBinContent(xbinc3, ybinc3);
    double yvaluec3 = h2DQ_Sn->GetYaxis()->GetBinCenter(ybinc3);
    sumc3 += contentc3 * yvaluec3; //calc sum val*weight
    sumsqc3 += contentc3 * yvaluec3 * yvaluec3; //calc val* weight²
    countc3 += contentc3;
    double contentc4 = h2DQ_D->GetBinContent(xbinc3, ybinc3);
    double yvaluec4 = h2DQ_D->GetYaxis()->GetBinCenter(ybinc3);
    sumc4 += contentc4 * yvaluec4;
    sumsqc4 += contentc4 * yvaluec4 * yvaluec4;
    countc4 += contentc4;
  }
  avgc3[xbinc3-1] = (countc3 > 0) ? sumc3 / countc3 : 0.0;
  double variancec3 = (countc3 > 1) ? (sumsqc3 - sumc3*sumc3/countc3) / (countc3 - 1) : 0.0;
  errc3[xbinc3-1] = (countc3 > 1) ? sqrt(variancec3/countc3) : 0.0;
  avgc4[xbinc4-1] = (countc4 > 0) ? sumc4 / countc4 : 0.0;
  double variancec4 = (countc4 > 1) ? (sumsqc4 - sumc4*sumc4/countc4) / (countc4 - 1) : 0.0;
  errc4[xbinc4-1] = (countc4 > 1) ? sqrt(variancec4/countc4) : 0.0;
}

for (int xbinc4 = 1; xbinc4 <= nbinsxc4; xbinc4++) {
  
  for (int ybinc4 = 1; ybinc4 <= h_cv_De->GetNbinsY(); ybinc4++) {
    
  }
  
}

double* avgcv = new double[nbinsxc4];
double* errcv = new double[nbinsxc4];
for (int i_avgcv = 1; i_avgcv <= nbinsxc4-1; i_avgcv++) {
    double x_cv = hist_v_Sn->GetXaxis()->GetBinCenter(i_avgcv);
    avgcv[i_avgcv]  = avgc3[i_avgcv] / avgc4[i_avgcv];
    double dc3 = 1/(avgc3[i_avgcv]);	
    double dc4 = (avgc3[i_avgcv]) / (avgc4[i_avgcv]*avgc4[i_avgcv]);	//partial deriv in (2nd var) denominatr --- needs to be squared when calculating errcQ
    errcv[i_avgcv] = sqrt( (dc3*dc3)  *(errc3[i_avgcv]*errc3[i_avgcv])  + (dc4*dc4) *(errc4[i_avgcv]*errc4[i_avgcv]) );
    if (x_cv<8.5 && x_cv>2.5){
    	plotcratio_v->SetPoint(i_avgcv-1,x_cv,avgcv[i_avgcv]);
	plotcratio_v->SetPointError(i_avgcv-1, 0, errcv[i_avgcv]);
    }
}


/////////////////////////////////////////////////////////////////////////////////
*/

    // Close the first ROOT file
    fileD->Close();
    fileSn->Close();
    

    
cR->cd(1);
cos_Qbis->GetXaxis()->SetTitleSize(0.05);
cos_Qbis->GetYaxis()->SetTitleSize(0.05);
cos_Qbis->SetMarkerSize(0.5);
cos_Qbis->SetMarkerStyle(21);
cos_Qbis->GetXaxis()->SetRangeUser(1,6.5);
cos_Qbis->GetXaxis()->SetTitle("Q^{2} " );
cos_Qbis->GetYaxis()->SetTitle("<cos(#Phi)>_{Sn} / <cos(#Phi)>_{De}");
cos_Qbis->Draw("AP");
cR->cd(2);
cos_vbis->GetXaxis()->SetTitleSize(0.05);
cos_vbis->GetYaxis()->SetTitleSize(0.05);
cos_vbis->GetXaxis()->SetRangeUser(2,9);
cos_vbis->SetMarkerSize(0.5);
cos_vbis->SetMarkerStyle(21);
cos_vbis->GetXaxis()->SetTitle("#nu " );
cos_vbis->GetYaxis()->SetTitle("<cos(#Phi)>_{Sn} / <cos(#Phi)>_{De}");
cos_vbis->Draw("AP");
cR->cd(3);
cos_zbis->GetXaxis()->SetTitleSize(0.05);
cos_zbis->GetYaxis()->SetTitleSize(0.05);
cos_zbis->SetMarkerSize(0.5);
cos_zbis->SetMarkerStyle(21);
cos_zbis->GetXaxis()->SetTitle("p_{t}^{2} " );
cos_zbis->GetYaxis()->SetTitle("<cos(#Phi)>_{Sn} / <cos(#Phi)>_{De}");
cos_zbis->Draw("AP");
cR->cd(4);
cos_Qter->GetXaxis()->SetTitleSize(0.05);
cos_Qter->GetYaxis()->SetTitleSize(0.05);
cos_Qter->SetMarkerSize(0.5);
cos_Qter->SetMarkerStyle(21);
cos_Qter->GetXaxis()->SetRangeUser(0.2,0.9);
cos_Qter->GetXaxis()->SetTitle("z " );
cos_Qter->GetYaxis()->SetTitle("<cos(#Phi)>_{Sn} / <cos(#Phi)>_{De}");
cos_Qter->Draw("AP");
    cR->SaveAs("cratio.pdf");
    cR->SaveAs("cratio.root");

    return 0;
}

