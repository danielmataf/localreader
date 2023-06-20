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




  for (Int_t i = 1; i <= 10; ++i)
    {
        double  Q2_D = h_QtestD->GetBinCenter(i);
        double  Q2_Sn = h_QtestSn->GetBinCenter(i);
        double  v_D = h_vtestD->GetBinCenter(i);
        double  v_Sn = h_vtestSn->GetBinCenter(i);
        double  z_D = h_ztestD->GetBinCenter(i);
        double  z_Sn = h_ztestSn->GetBinCenter(i);
        double  P_t_D = h_pttestD->GetBinCenter(i);
        double  P_t_Sn = h_pttestSn->GetBinContent(i);
        double  phih_D = h_phiD->GetBinContent(i);
        double  phih_Sn = h_phiSn->GetBinContent(i);
        
        hist_sratio.w_sQ_Sn->Fill(Q2_Sn, sin(phih_Sn));
        hist_sratio.w_sQ_De->Fill(Q2_D, sin(phih_D));
        hist_sratio.w_sv_Sn->Fill(v_Sn, sin(phih_Sn));
        hist_sratio.w_sv_De->Fill(v_D, sin(phih_D));
        hist_sratio.w_sz_Sn->Fill(z_Sn, sin(phih_Sn));
        hist_sratio.w_sz_De->Fill(z_D, sin(phih_D));
        hist_sratio.w_spt_De->Fill(P_t_Sn, sin(phih_Sn));
        hist_sratio.w_spt_De->Fill(P_t_D, sin(phih_D));
    }




    // Access properties of h_Q1
    //std::cout>>numBins1>>std::endl;  
 
for (int i_q=1; i_q<=hist_sratio.w_sQ_Sn->GetNbinsX(); i_q++) {
        double x_q_w1 = hist_sratio.w_sQ_Sn->GetXaxis()->GetBinCenter(i_q);		
        double y_q_wSn1 = hist_sratio.w_sQ_Sn->GetBinContent(i_q);				
        double x_q1 = h_QtestSn->GetXaxis()->GetBinCenter(i_q);			
        double y_qSn1 = h_QtestSn->GetBinContent(i_q) ;
	    double y_q_wDe1 = hist_sratio.w_sQ_De->GetBinContent(i_q);
	    double y_qDe1 = h_QtestD->GetBinContent(i_q) ;
	    double meansQSn1 = y_q_wSn1 / y_qSn1;
	    double meansQDe1 = y_q_wDe1 / y_qDe1; 
	//cout<<meandptQSn - meandptQDe <<"meanDPTQ"<<endl;
	//double Rq_err=  (1/y_q_Sn) + (1/y_q_De) + (1/y_q_Sn_e) + (1/y_q_De_e); 
		if(x_q_w1>1.5){ 
        	    sin_Q->SetPoint(i_q-1, x_q_w1,meansQSn1/ meansQDe1);
    	}
	//}
}

for (int i_v=1; i_v<=hist_sratio.w_sQ_Sn->GetNbinsX(); i_v++) {
        double x_v_w1 = hist_sratio.w_sQ_Sn->GetXaxis()->GetBinCenter(i_v);		
        double y_v_wSn1 = hist_sratio.w_sQ_Sn->GetBinContent(i_v);				
        double x_v1 = h_vtestSn->GetXaxis()->GetBinCenter(i_v);			
        double y_vSn1 = h_vtestSn->GetBinContent(i_v) ;
	    double y_v_wDe1 = hist_sratio.w_sQ_De->GetBinContent(i_v);
	    double y_vDe1 = h_vtestD->GetBinContent(i_v) ;
	    double meansvSn1 = y_v_wSn1 / y_vSn1;
	    double meansvDe1 = y_v_wDe1 / y_vDe1; 
	//cout<<meandptQSn - meandptQDe <<"meanDPTQ"<<endl;
	//double Rq_err=  (1/y_v_Sn) + (1/y_v_De) + (1/y_v_Sn_e) + (1/y_v_De_e); 
		//if(x_v_w<8.5 && x_v_w>2.5){ 
        	    sin_v->SetPoint(i_v-1, x_v_w1,meansvSn1/ meansvDe1);
    	//}
	//}
}

for (int i_z=1; i_z<=hist_sratio.w_sQ_Sn->GetNbinsX(); i_z++) {
        double x_z_w1 = hist_sratio.w_sQ_Sn->GetXaxis()->GetBinCenter(i_z);		
        double y_z_wSn1 = hist_sratio.w_sQ_Sn->GetBinContent(i_z);				
        double x_z1 = h_ztestSn->GetXaxis()->GetBinCenter(i_z);			
        double y_zSn1 = h_ztestSn->GetBinContent(i_z) ;
	    double y_z_wDe1 = hist_sratio.w_sQ_De->GetBinContent(i_z);
	    double y_zDe1 = h_ztestD->GetBinContent(i_z) ;
	    double meanszSn1 = y_z_wSn1 / y_zSn1;
	    double meanszDe1 = y_z_wDe1 / y_zDe1; 
	//cout<<meandptQSn - meandptQDe <<"meanDPTQ"<<endl;
	//double Rq_err=  (1/y_z_Sn) + (1/y_z_De) + (1/y_z_Sn_e) + (1/y_z_De_e); 
		//if(x_z_w<8.5 && x_z_w>2.5){ 
        	    sin_z->SetPoint(i_z-1, x_z_w1,meanszSn1/ meanszDe1);
    	//}
	//}
}

for (int i_pt=1; i_pt<=hist_sratio.w_sQ_Sn->GetNbinsX(); i_pt++) {
        double x_pt_w1 = hist_sratio.w_sQ_Sn->GetXaxis()->GetBinCenter(i_pt);		
        double y_pt_wSn1 = hist_sratio.w_sQ_Sn->GetBinContent(i_pt);				
        double x_pt1 = h_pttestSn->GetXaxis()->GetBinCenter(i_pt);			
        double y_ptSn1 = h_pttestSn->GetBinContent(i_pt) ;
	    double y_pt_wDe1 = hist_sratio.w_sQ_De->GetBinContent(i_pt);
	    double y_ptDe1 = h_pttestD->GetBinContent(i_pt) ;
	    double meansptSn1 = y_pt_wSn1 / y_ptSn1;
	    double meansptDe1 = y_pt_wDe1 / y_ptDe1; 
	//cout<<meandptQSn - meandptQDe <<"meanDPTQ"<<endl;
	//double Rq_err=  (1/y_pt_Sn) + (1/y_pt_De) + (1/y_pt_Sn_e) + (1/y_pt_De_e); 
		//if(x_pt_w<8.5 && x_pt_w>2.5){ 
        	    sin_pt->SetPoint(i_pt-1, x_pt_w1,meansptSn1/ meansptDe1);
    	//}
	//}
}






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

    return 0;
}

