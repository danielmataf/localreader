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

void calculateCRat(TH1F* hD, TH1F* h_weightD, TH1F* h_sqD,TH1F* hSn, TH1F* h_weightSn,TH1F* h_sqSn, TGraphErrors* g_result) {
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
        double crat_point = (wavg_D > 0) ? wavg_Sn/wavg_D : 0.0; 
        double variance_D = Psq_D/N_D - pow(wavg_D,2);     
        double variance_Sn=   Psq_Sn/N_Sn - pow(wavg_Sn,2) ;
        cout<<"varD= "<< variance_D<<";   varSn= "<< variance_Sn<<endl;
        double Err_D  = sqrt(variance_D/N_D);
        double Err_Sn = sqrt(variance_Sn/N_Sn);
        //double unc = (wavg_Sn/wavg_D) *sqrt( pow( Err_D/wavg_D,2) + pow(Err_Sn/wavg_Sn ,2) );
        //if (unc>5){
        //    unc = 0.0;
        //    crat_point = 0.0;
        //}

        g_result->SetPoint(bin - 1, x, crat_point);
        g_result->SetPointError(bin - 1, 0, Err_D);
    }
    cout<<"............"<<endl;
}

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
    TGraphErrors *cos_vter = new TGraphErrors();
    TGraphErrors *cos_zter = new TGraphErrors();
    TGraphErrors *cos_pter = new TGraphErrors();

    // Open the first ROOT file
    //TFile* file1 = new TFile("build/output1.root", "READ");
    //TFile* fileD = new TFile("../files2read/testOutput_D.root", "READ");
    //TFile* fileSn = new TFile("../files2read/testOutput_Sn.root", "READ");
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
    TH1F* h_ptweightD = dynamic_cast<TH1F*>(fileD->Get("weightedpcos"));
    TH1F* h_ptweightSn = dynamic_cast<TH1F*>(fileSn->Get("weightedpcos"));

    TH1F* h_QsqD = dynamic_cast<TH1F*>(fileSn->Get("VarianceQcos"));
    TH1F* h_QsqSn = dynamic_cast<TH1F*>(fileSn->Get("VarianceQcos"));
    TH1F* h_vsqD = dynamic_cast<TH1F*>(fileD->Get("Variancevcos"));
    TH1F* h_vsqSn = dynamic_cast<TH1F*>(fileSn->Get("Variancevcos"));
    TH1F* h_zsqD = dynamic_cast<TH1F*>(fileD->Get("Variancezcos"));
    TH1F* h_zsqSn = dynamic_cast<TH1F*>(fileSn->Get("Variancezcos"));
    TH1F* h_ptsqD = dynamic_cast<TH1F*>(fileD->Get("Variancepcos"));
    TH1F* h_ptsqSn = dynamic_cast<TH1F*>(fileSn->Get("Variancepcos"));

//  ward to see if the histogram is existent
/*
    if (!h_ptweightD) {
        std::cerr << "Error: Histogram not found in the file." << std::endl;
        fileD->Close();
        return 0;
    }
*/
//
    calculateCRat( h_QtestD, h_QweightD, h_QsqD,h_QtestSn, h_QweightSn,h_QsqD, cos_Q);
    calculateCRat( h_vtestD, h_vweightD, h_vsqD,h_vtestSn, h_vweightSn,h_vsqD, cos_v);
    calculateCRat( h_ztestD, h_zweightD, h_zsqD,h_ztestSn, h_zweightSn,h_zsqD, cos_z);
    calculateCRat( h_pttestD, h_ptweightD, h_ptsqD,h_pttestSn, h_ptweightSn,h_ptsqD,cos_pt);



/*
    int nubinsQ = h_QtestD->GetNbinsX(); 
    int nubinsv = h_vtestD->GetNbinsX(); 
    int nubinsz = h_ztestD->GetNbinsX(); 
    int nubinsp = h_pttestD->GetNbinsX(); 
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
        double Psq_D = h_QsqD->GetBinContent(binQ);
        double Psq_Sn = h_QsqSn->GetBinContent(binQ);


        //<cos (phi_D)> = cos_D/N_D
        double wavg_D  = (N_D > 0) ? cos_D/N_D : 0.0;
        double wavg_Sn = (N_Sn > 0) ? cos_Sn/N_Sn : 0.0;
        double cratio_point = (wavg_D > 0) ? wavg_Sn/wavg_D : 0.0; 
        //cout<<wavg_Sn<<"/" <<wavg_D <<" = " <<cratio_point<<"   cratiopoint"<< endl;

        double varianceQ_D = Psq_D/N_D - pow(wavg_D,2);     
        double varianceQ_Sn= Psq_Sn/N_Sn - pow(wavg_Sn,2) ;
        double ErrQ_D  = sqrt(varianceQ_D);
        double ErrQ_Sn = sqrt(varianceQ_Sn);
        double unc = (wavg_Sn/wavg_D) *sqrt( pow( ErrQ_D/wavg_D,2) + pow(ErrQ_Sn/wavg_Sn ,2) );
        cos_Qter->SetPoint(binQ -1, x_Q_Sn, cratio_point );
        cos_Qter->SetPointError(binQ-1, 0, unc);
        //cout<<Errz_D<<"Errrrr....."<<endl;

    }

    for (int binv = 1; binv <= nubinsv; binv++) {
        double N_D = h_vtestD->GetBinContent(binv);
        double cos_D = h_vweightD->GetBinContent(binv);
        double N_Sn = h_vtestSn->GetBinContent(binv);        
        double cos_Sn = h_vweightSn->GetBinContent(binv);  
        double x_Q_Sn = h_vtestSn->GetXaxis()->GetBinCenter(binv);
        double Psq_D = h_vsqD->GetBinContent(binv); 
        double Psq_Sn = h_vsqSn->GetBinContent(binv);

        //<cos (phi_D)> = cos_D/N_D
        double wavg_D  = (N_D > 0) ? cos_D/N_D : 0.0;
        double wavg_Sn = (N_Sn > 0) ? cos_Sn/N_Sn : 0.0;
        double cratio_vpoint = (wavg_D > 0) ? wavg_Sn/wavg_D : 0.0; 
        //cout<<wavg_Sn<<"/" <<wavg_D <<" = " <<cratio_vpoint<<"   cratiopoint(v)"<< endl;
        double variancev_D = Psq_D/N_D - pow(wavg_D,2);     
        double variancev_Sn= Psq_Sn/N_Sn - pow(wavg_Sn,2) ;
        double Errv_D  = sqrt(variancev_D);
        double Errv_Sn = sqrt(variancev_Sn);
        double unc = (wavg_Sn/wavg_D) *sqrt( pow( Errv_D/wavg_D,2) + pow(Errv_Sn/wavg_Sn ,2) );
        cos_vter->SetPoint(binv-1, x_Q_Sn, cratio_vpoint );
        cos_vter->SetPointError(binv-1, 0, unc);
        cout<<unc<<" unc .........-" <<endl;
    }

    

    for (int binz = 1; binz <= nubinsz; binz++) {
        double N_D = h_ztestD->GetBinContent(binz);
        double cos_D = h_zweightD->GetBinContent(binz);
        double N_Sn = h_ztestSn->GetBinContent(binz);        
        double cos_Sn = h_zweightSn->GetBinContent(binz);  
        double x_Q_Sn = h_ztestSn->GetXaxis()->GetBinCenter(binz);
        double Psq_D = h_vsqD->GetBinContent(binz);         
        double Psq_Sn = h_vsqSn->GetBinContent(binz);       

        //<cos (phi_D)> = cos_D/N_D
        double wavg_D  = (N_D > 0) ? cos_D/N_D : 0.0;
        double wavg_Sn = (N_Sn > 0) ? cos_Sn/N_Sn : 0.0;
        double cratio_zpoint = (wavg_D > 0) ? wavg_Sn/wavg_D : 0.0; 
        //cout<<wavg_Sn<<"/" <<wavg_D <<" = " <<cratio_zpoint<<"   cratiopoint(z)"<< endl;

        double variancez_D = Psq_D/N_D - pow(wavg_D,2);     
        double variancez_Sn= Psq_Sn/N_Sn - pow(wavg_Sn,2) ;
        double Errz_D  = sqrt(variancez_D);
        double Errz_Sn = sqrt(variancez_Sn);
        double unc = (wavg_Sn/wavg_D) *sqrt( pow( Errz_D/wavg_D,2) + pow(Errz_Sn/wavg_Sn ,2) );
        cos_zter->SetPoint(binz-1, x_Q_Sn, cratio_zpoint );
        cos_zter->SetPointError(binz-1, 0, unc);


    }

    for (int binp = 1; binp <= nubinsp; binp++) {
        double N_D = h_pttestD->GetBinContent(binp);
        double cos_D = h_ptweightD->GetBinContent(binp);
        double N_Sn = h_pttestSn->GetBinContent(binp);        
        double cos_Sn = h_ptweightSn->GetBinContent(binp);  
        double x_Q_Sn = h_pttestSn->GetXaxis()->GetBinCenter(binp);
        double Psq_D = h_vsqD->GetBinContent(binp);         
        double Psq_Sn = h_vsqSn->GetBinContent(binp);       

        //<cos (phi_D)> = cos_D/N_D
        double wavg_D  = (N_D > 0) ? cos_D/N_D : 0.0;
        double wavg_Sn = (N_Sn > 0) ? cos_Sn/N_Sn : 0.0;
        double cratio_ppoint = (wavg_D > 0) ? wavg_Sn/wavg_D : 0.0; 
        //cout<<wavg_Sn<<"/" <<wavg_D <<" = " <<cratio_ppoint<<"   cratiopoint(p)"<< endl;
        double variancep_D = Psq_D/N_D - pow(wavg_D,2);     
        double variancep_Sn= Psq_Sn/N_Sn - pow(wavg_Sn,2) ;
        double Errp_D  = sqrt(variancep_D);
        double Errp_Sn = sqrt(variancep_Sn);
        double unc = (wavg_Sn/wavg_D) *sqrt( pow( Errp_D/wavg_D,2) + pow(Errp_Sn/wavg_Sn ,2) );
        cos_pter->SetPoint(binp-1, x_Q_Sn, cratio_ppoint );
        cos_pter->SetPointError(binp-1, 0, 0);

    }
    ///////////weighted TH1F test/////////////////
    //////////////////////////////////////////////


*/

    // Close the first ROOT file
    fileD->Close();
    fileSn->Close();
    

    


cR->cd(1);
cos_Q->GetXaxis()->SetTitleSize(0.05);
cos_Q->GetYaxis()->SetTitleSize(0.05);
cos_Q->SetMarkerSize(0.5);
cos_Q->SetMarkerStyle(21);
cos_Q->GetXaxis()->SetRangeUser(1,6.5);
cos_Q->GetXaxis()->SetTitle("Q^{2} " );
cos_Q->GetYaxis()->SetTitle("<cos(#Phi)>_{Sn} / <cos(#Phi)>_{De}");
cos_Q->Draw("APE");
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
cos_z->GetXaxis()->SetTitle("z");
cos_z->GetYaxis()->SetTitle("<cos(#Phi)>_{Sn} / <cos(#Phi)>_{De}");
cos_z->Draw("AP");
cR->cd(4);
cos_pt->GetXaxis()->SetTitleSize(0.05);
cos_pt->GetYaxis()->SetTitleSize(0.05);
cos_pt->SetMarkerSize(0.5);
cos_pt->SetMarkerStyle(21);
cos_pt->GetXaxis()->SetRangeUser(0.2,0.9);
cos_pt->GetXaxis()->SetTitle("p_{t}^{2} " );
cos_pt->GetYaxis()->SetTitle("<cos(#Phi)>_{Sn} / <cos(#Phi)>_{De}");
cos_pt->Draw("AP");
    //cR->SaveAs("cosD.pdf");
    //cR->SaveAs("cosD.root");

    return 0;
}

