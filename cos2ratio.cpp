// Include necessary ROOT headers
//Run and compile with 
//g++ -o cos2ratio cos2ratio.cpp `root-config --cflags --libs`
// ./cos2ratio

////////////////////////////////////////////////////////////////
//Code separated from the first one. 
// June 2023 
// Need to test " plot c2Ratio.pdf" and "testc2Ratio.pdf" recovered from recoveryposCLAS.cpp
//  c2Ratio = <cphi>_Sn/ <cphi>_D        //c2Ratio(Q2,v,z,pt)
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


            //w_c2Q_Sn->Fill(Q2,cos2(phih));


int main() {
    h_c2ratio hist_c2Ratio;
    TCanvas* c2R=new TCanvas("ratio","ratio"); //c2Reate new canvas
    c2R->Divide(2,2);
    TGraphErrors *cos2_Q = new TGraphErrors();
    TGraphErrors *cos2_v = new TGraphErrors();
    TGraphErrors *cos2_z = new TGraphErrors();
    TGraphErrors *cos2_pt = new TGraphErrors();
    TGraphErrors *cos2_Qbis = new TGraphErrors();
    TGraphErrors *cos2_vbis = new TGraphErrors();
    TGraphErrors *cos2_zbis = new TGraphErrors();
    TGraphErrors *cos2_pbis  = new TGraphErrors();
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
    TH2F* h2DQ_D = dynamic_cast<TH2F*>(fileD->Get("cos2Q"));  // Retrieve the TH2F histogram
    TH2F* h2Dv_D = dynamic_cast<TH2F*>(fileD->Get("cos2v"));  
    TH2F* h2Dz_D = dynamic_cast<TH2F*>(fileD->Get("cos2z"));
    TH2F* h2Dp_D = dynamic_cast<TH2F*>(fileD->Get("cos2p"));
    TH2F* h2DQ_Sn = dynamic_cast<TH2F*>(fileSn->Get("cos2Q"));
    TH2F* h2Dv_Sn = dynamic_cast<TH2F*>(fileSn->Get("cos2v"));
    TH2F* h2Dz_Sn = dynamic_cast<TH2F*>(fileSn->Get("cos2z"));
    TH2F* h2Dp_Sn = dynamic_cast<TH2F*>(fileSn->Get("cos2p"));
    //pt is already squared in the histogrammmm





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
        
        hist_c2Ratio.w_c2Q_Sn->Fill(Q2_Sn, cos(2* phih_Sn));
        hist_c2Ratio.w_c2Q_De->Fill(Q2_D, cos(2* phih_D));
        hist_c2Ratio.w_c2v_Sn->Fill(v_Sn, cos(2* phih_Sn));
        hist_c2Ratio.w_c2v_De->Fill(v_D, cos(2* phih_D));
        hist_c2Ratio.w_c2z_Sn->Fill(z_Sn, cos(2* phih_Sn));
        hist_c2Ratio.w_c2z_De->Fill(z_D, cos(2* phih_D));
        hist_c2Ratio.w_c2pt_De->Fill(P_t_Sn, cos(2* phih_Sn));
        hist_c2Ratio.w_c2pt_De->Fill(P_t_D, cos(2* phih_D));
    }




    // Access properties of h_Q1
    //std::cout>>numBins1>>std::endl;  
 
for (int i_q=1; i_q<=hist_c2Ratio.w_c2Q_Sn->GetNbinsX(); i_q++) {
        double x_q_w1 = hist_c2Ratio.w_c2Q_Sn->GetXaxis()->GetBinCenter(i_q);		
        double y_q_wSn1 = hist_c2Ratio.w_c2Q_Sn->GetBinContent(i_q);				
        double x_q1 = h_QtestSn->GetXaxis()->GetBinCenter(i_q);			
        double y_qSn1 = h_QtestSn->GetBinContent(i_q) ;
	    double y_q_wDe1 = hist_c2Ratio.w_c2Q_De->GetBinContent(i_q);
	    double y_qDe1 = h_QtestD->GetBinContent(i_q) ;
	    double meancQSn1 = y_q_wSn1 / y_qSn1;
	    double meancQDe1 = y_q_wDe1 / y_qDe1; 
	//cout<<meandptQSn - meandptQDe <<"meanDPTQ"<<endl;
	//double Rq_err=  (1/y_q_Sn) + (1/y_q_De) + (1/y_q_Sn_e) + (1/y_q_De_e); 
		if(x_q_w1>1.5){ 
        	    cos2_Q->SetPoint(i_q-1, x_q_w1,meancQSn1/ meancQDe1);
    	}
	//}
}

for (int i_v=1; i_v<=hist_c2Ratio.w_c2Q_Sn->GetNbinsX(); i_v++) {
        double x_v_w1 = hist_c2Ratio.w_c2Q_Sn->GetXaxis()->GetBinCenter(i_v);		
        double y_v_wSn1 = hist_c2Ratio.w_c2Q_Sn->GetBinContent(i_v);				
        double x_v1 = h_vtestSn->GetXaxis()->GetBinCenter(i_v);			
        double y_vSn1 = h_vtestSn->GetBinContent(i_v) ;
	    double y_v_wDe1 = hist_c2Ratio.w_c2Q_De->GetBinContent(i_v);
	    double y_vDe1 = h_vtestD->GetBinContent(i_v) ;
	    double meancvSn1 = y_v_wSn1 / y_vSn1;
	    double meancvDe1 = y_v_wDe1 / y_vDe1; 
	//cout<<meandptQSn - meandptQDe <<"meanDPTQ"<<endl;
	//double Rq_err=  (1/y_v_Sn) + (1/y_v_De) + (1/y_v_Sn_e) + (1/y_v_De_e); 
		//if(x_v_w<8.5 && x_v_w>2.5){ 
        	    cos2_v->SetPoint(i_v-1, x_v_w1,meancvSn1/ meancvDe1);
    	//}
	//}
}

for (int i_z=1; i_z<=hist_c2Ratio.w_c2Q_Sn->GetNbinsX(); i_z++) {
        double x_z_w1 = hist_c2Ratio.w_c2Q_Sn->GetXaxis()->GetBinCenter(i_z);		
        double y_z_wSn1 = hist_c2Ratio.w_c2Q_Sn->GetBinContent(i_z);				
        double x_z1 = h_ztestSn->GetXaxis()->GetBinCenter(i_z);			
        double y_zSn1 = h_ztestSn->GetBinContent(i_z) ;
	    double y_z_wDe1 = hist_c2Ratio.w_c2Q_De->GetBinContent(i_z);
	    double y_zDe1 = h_ztestD->GetBinContent(i_z) ;
	    double meanczSn1 = y_z_wSn1 / y_zSn1;
	    double meanczDe1 = y_z_wDe1 / y_zDe1; 
	//cout<<meandptQSn - meandptQDe <<"meanDPTQ"<<endl;
	//double Rq_err=  (1/y_z_Sn) + (1/y_z_De) + (1/y_z_Sn_e) + (1/y_z_De_e); 
		//if(x_z_w<8.5 && x_z_w>2.5){ 
        	    cos2_z->SetPoint(i_z-1, x_z_w1,meanczSn1/ meanczDe1);
    	//}
	//}
}









for (int i_pt=1; i_pt<=hist_c2Ratio.w_c2Q_Sn->GetNbinsX(); i_pt++) {
        double x_pt_w1 = hist_c2Ratio.w_c2Q_Sn->GetXaxis()->GetBinCenter(i_pt);		
        double y_pt_wSn1 = hist_c2Ratio.w_c2Q_Sn->GetBinContent(i_pt);				
        double x_pt1 = h_pttestSn->GetXaxis()->GetBinCenter(i_pt);			
        double y_ptSn1 = h_pttestSn->GetBinContent(i_pt) ;
	    double y_pt_wDe1 = hist_c2Ratio.w_c2Q_De->GetBinContent(i_pt);
	    double y_ptDe1 = h_pttestD->GetBinContent(i_pt) ;
	    double meancptSn1 = y_pt_wSn1 / y_ptSn1;
	    double meancptDe1 = y_pt_wDe1 / y_ptDe1; 
	//cout<<meandptQSn - meandptQDe <<"meanDPTQ"<<endl;
	//double Rq_err=  (1/y_pt_Sn) + (1/y_pt_De) + (1/y_pt_Sn_e) + (1/y_pt_De_e); 
		//if(x_pt_w<8.5 && x_pt_w>2.5){ 
        	    cos2_pt->SetPoint(i_pt-1, x_pt_w1,meancptSn1/ meancptDe1);
    	//}
	//}
}


    // Close the first ROOT file
    fileD->Close();
    fileSn->Close();
    


    c2R->cd(1);
cos2_Qbis->GetXaxis()->SetTitleSize(0.05);
cos2_Qbis->GetYaxis()->SetTitleSize(0.05);
cos2_Qbis->SetMarkerSize(0.5);
cos2_Qbis->SetMarkerStyle(21);
cos2_Qbis->GetXaxis()->SetRangeUser(1,6.5);
cos2_Qbis->GetXaxis()->SetTitle("Q^{2} " );
cos2_Qbis->GetYaxis()->SetTitle("<cos(  #phih)>_{Sn} / <cos(  #phih)>_{D}");
cos2_Qbis->Draw("APE");



    c2R->SaveAs("cRatiotest.pdf");
    c2R->SaveAs("cRatiotest.root");

c2R->cd(1);
cos2_Q->GetXaxis()->SetTitleSize(0.05);
cos2_Q->GetYaxis()->SetTitleSize(0.05);
cos2_Q->SetMarkerSize(0.5);
cos2_Q->SetMarkerStyle(21);
cos2_Q->GetXaxis()->SetRangeUser(1,6.5);
cos2_Q->GetXaxis()->SetTitle("Q^{2} " );
cos2_Q->GetYaxis()->SetTitle("<cos( 2 #phih)>_{Sn} / <cos( 2 #phih)>_{D}");
cos2_Q->Draw("AP");

c2R->cd(2);
cos2_v->GetXaxis()->SetTitleSize(0.05);
cos2_v->GetYaxis()->SetTitleSize(0.05);
cos2_v->GetXaxis()->SetRangeUser(2,9);
cos2_v->SetMarkerSize(0.5);
cos2_v->SetMarkerStyle(21);
cos2_v->GetXaxis()->SetTitle("#nu " );
cos2_v->GetYaxis()->SetTitle("<cos(2 #phih)>_{Sn} / <cos(2 #phih)>_{D}");
//cos_v->Draw("AP");

c2R->cd(3);
cos2_z->GetXaxis()->SetTitleSize(0.05);
cos2_z->GetYaxis()->SetTitleSize(0.05);
cos2_z->SetMarkerSize(0.5);
cos2_z->SetMarkerStyle(21);
cos2_z->GetXaxis()->SetTitle("p_{t}^{2} " );
cos2_z->GetYaxis()->SetTitle("<cos(2 #phih)>_{Sn} / <cos(2 #phih)>_{D}");
//cos_z->Draw("AP");
c2R->cd(4);
cos2_pt->GetXaxis()->SetTitleSize(0.05);
cos2_pt->GetYaxis()->SetTitleSize(0.05);
cos2_pt->SetMarkerSize(0.5);
cos2_pt->SetMarkerStyle(21);
cos2_pt->GetXaxis()->SetRangeUser(0.2,0.9);
cos2_pt->GetXaxis()->SetTitle("z " );
cos2_pt->GetYaxis()->SetTitle("<cos(2 #phih)>_{Sn} / <cos(2 #phih)>_{D}");
cos2_pt->Draw("AP");

    c2R->SaveAs("testc2Ratio.pdf");
    c2R->SaveAs("testc2Ratio.root");



    return 0;
}

