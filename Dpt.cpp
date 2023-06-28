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
    TGraphErrors *Dpt_Q = new TGraphErrors();
    TGraphErrors *Dpt_v = new TGraphErrors();
    TGraphErrors *Dpt_z = new TGraphErrors();
    // Open the first ROOT file
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
    //maybe not even need to dynamic cast, maybe hjust enough with calling the class
    //-> Yes, dynamic cast is needed, to make the difference between Sn and De
    //-> an adjustment has to be done to sumW2 histograms in order to properly recover them making the distinction on them. 
    hist_Dpt.w_ptQ_Sn->Sumw2();
    hist_Dpt.w_ptQ_De->Sumw2();
    hist_Dpt.w_ptv_Sn->Sumw2();
    hist_Dpt.w_ptv_De->Sumw2();
    hist_Dpt.w_ptz_Sn->Sumw2();
    hist_Dpt.w_ptz_De->Sumw2();
    //hist_Dpt.h_QtestD->Sumw2();
    for (Int_t i = 1; i <= 10; ++i)
    {
        Double_t Q2_D = h_QtestD->GetBinCenter(i);
        Double_t Q2_Sn = h_QtestSn->GetBinCenter(i);
        //change name to the following TBD :
        Double_t v_D = h_vtestD->GetBinCenter(i);
        Double_t v_Sn = h_vtestSn->GetBinCenter(i);
        Double_t z_D = h_ztestD->GetBinCenter(i);
        Double_t z_Sn = h_ztestSn->GetBinCenter(i);
        Double_t P_t_D = h_pttestD->GetBinCenter(i);
        Double_t P_t_Sn = h_pttestSn->GetBinContent(i);
        
        hist_Dpt.w_ptQ_Sn->Fill(Q2_Sn, P_t_Sn * P_t_Sn);
        hist_Dpt.w_ptQ_De->Fill(Q2_D, P_t_D * P_t_D);
        hist_Dpt.w_ptv_Sn->Fill(v_Sn, P_t_Sn * P_t_Sn);
        hist_Dpt.w_ptv_De->Fill(v_D, P_t_D * P_t_D);
        hist_Dpt.w_ptz_Sn->Fill(z_Sn, P_t_Sn * P_t_Sn);
        hist_Dpt.w_ptz_De->Fill(z_D, P_t_D * P_t_D);
        cout<< Q2_D<<";"  <<Q2_Sn<< "; "<<v_D <<";"<<v_Sn <<"; "<<z_D<<";"<<z_Sn<<endl;
        cout<<P_t_D<<";"<<P_t_Sn<<endl;
        cout<<" ---------------- "<<endl; 
    }



    //TH1F* h_pttestD = dynamic_cast<TH1F*>(fileD->Get("P_t"));
    //TH1F* h_pttestSn = dynamic_cast<TH1F*>(fileSn->Get("P_t"));

    // Access properties of h_Q1
    //std::cout>>numBins1>>std::endl;  
    // ...
    //hist_Dpt.w_ptQ_Sn;

    cout << hist_Dpt.w_ptQ_Sn->GetXaxis()->GetBinCenter(5) << endl;
    //this outbut is just displaying the binning, not actual content the sumW hists were not filled with Sn/D distinction. Empty. no entries
    
    
    for (int i_q=1; i_q<=hist_Dpt.w_ptQ_Sn->GetNbinsX(); i_q++) 
    {
        double x_q_w = hist_Dpt.w_ptQ_Sn->GetXaxis()->GetBinCenter(i_q);		
        double y_q_wSn = hist_Dpt.w_ptQ_Sn->GetBinContent(i_q);				
        double x_q = hist_Dpt.w_ptQ_De->GetXaxis()->GetBinCenter(i_q);			
        double y_qSn = h_pttestSn->GetBinContent(i_q) ;
    	double y_q_wDe = hist_Dpt.w_ptQ_De->GetBinContent(i_q);
	    double y_qDe = h_pttestD->GetBinContent(i_q) ;
	    double meandptQSn =(y_qSn > 0) ?  y_q_wSn / y_qSn : 0.0;
	    double meandptQDe = ( y_qDe > 0) ?  y_q_wDe / y_qDe : 0.0;
	    double Rq_err=  (1/y_q_Sn) + (1/y_q_De) + (1/y_q_Sn_e) + (1/y_q_De_e);
      	cout<<Rq_err <<" = meanDPTQ_ERR"<<endl;
 
   	    //WE HAD A CNDITION ON X VALUE HERE... USELESS NOWbc zere not retrieving xvalues --- just make sure you do the cuts properly on main.cpp
        Dpt_Q->SetPoint(i_q-1, x_q_w,meandptQSn - meandptQDe);
    }
    for (int i_v=1; i_v<=hist_Dpt.w_ptv_Sn->GetNbinsX(); i_v++) 
    {
        double x_v_w = hist_Dpt.w_ptv_Sn->GetXaxis()->GetBinCenter(i_v);		
        double y_v_wSn = hist_Dpt.w_ptv_Sn->GetBinContent(i_v);				
        double x_v = hist_Dpt.w_ptv_De->GetXaxis()->GetBinCenter(i_v);			
        double y_vSn = h_pttestSn->GetBinContent(i_v) ;
    	double y_v_wDe = hist_Dpt.w_ptv_De->GetBinContent(i_v);
	    double y_vDe = h_pttestD->GetBinContent(i_v) ;
	    double meandptvSn =(y_vSn > 0) ?  y_v_wSn / y_vSn : 0.0;
	    double meandptvDe = ( y_vDe > 0) ?  y_v_wDe / y_vDe : 0.0;
	//cout<<meandptQSn - meandptQDe <<"meanDPTQ"<<endl;
	//double Rq_err=  (1/y_v_Sn) + (1/y_v_De) + (1/y_v_Sn_e) + (1/y_v_De_e); 
   	    //WE HAD A CNDITION ON X VALUE HERE... USELESS NOWbc zere not retrieving xvalues --- just make sure you do the cuts properly on main.cpp
        Dpt_v->SetPoint(i_v-1, x_v_w,meandptvSn - meandptvDe);
    }
    for (int i_z=1; i_z<=hist_Dpt.w_ptv_Sn->GetNbinsX(); i_z++) 
    {
        double x_z_w = hist_Dpt.w_ptz_Sn->GetXaxis()->GetBinCenter(i_z);		
        double y_z_wSn = hist_Dpt.w_ptz_Sn->GetBinContent(i_z);				
        double x_z = hist_Dpt.w_ptz_De->GetXaxis()->GetBinCenter(i_z);			
        double y_zSn = h_pttestSn->GetBinContent(i_z) ;
    	double y_z_wDe = hist_Dpt.w_ptz_De->GetBinContent(i_z);
	    double y_zDe = h_pttestD->GetBinContent(i_z) ;
	    double meandptzSn =(y_zSn > 0) ?  y_z_wSn / y_zSn : 0.0;
	    double meandptzDe = ( y_zDe > 0) ?  y_z_wDe / y_zDe : 0.0;
	//cout<<meandptQSn - meandptQDe <<"meanDPTQ"<<endl;
	//double Rq_err=  (1/y_z_Sn) + (1/y_z_De) + (1/y_z_Sn_e) + (1/y_z_De_e); 
   	    //WE HAD A CNDITION ON X VALUE HERE... USELESS NOWbc zere not retrieving xvalues --- just make sure you do the cuts properly on main.cpp
        Dpt_z->SetPoint(i_z-1, x_z_w,meandptzSn - meandptzDe);
    }
    

    // Close the first ROOT file
    fileD->Close();
    fileSn->Close();
    cDpt->cd(1);
    //R_Q->SetTitle(title);
    //R_Q->GetXaxis()->SetTitle("Q{2} (GeV^{2}) " );
    //ctest->GetYaxis()->SetTitle("R^{#Pi+}_{Sn}");
    Dpt_Q->Draw("AP");
    cDpt->cd(2);
    Dpt_v->Draw("AP");
    cDpt->cd(3);
    Dpt_z->Draw("AP");

    cDpt->SaveAs("dpttest.pdf");
    cDpt->SaveAs("dpttest.root");

    return 0;
}

