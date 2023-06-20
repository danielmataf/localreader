// Include necessary ROOT headers
//Run and compile with 
//g++ -o multipratio MultipRatio.cpp `root-config --cflags --libs`
// ./multipratio




////////////////////////////////////////////////////////////////
//Code separated from the first one. 
//It contains only a function that determines (maybe) the multipratio from a given histogram found in a root file 
// May 2023
// June 2023 
//  a condition has been added to plot thins according to kinematic cuts. Do not always consider it. also it is not the same cut for every variable (depends on the cut on that variable.)
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
using namespace std;



int main() {
    TCanvas* cR=new TCanvas("ratio","ratio"); //create new canvas
    cR->Divide(2,2);
    TGraphErrors *R_Q = new TGraphErrors();
    TGraphErrors *R_v = new TGraphErrors();
    TGraphErrors *R_z = new TGraphErrors();
    TGraphErrors *R_pt = new TGraphErrors();
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

    // Access properties of h_Q1
    //std::cout>>numBins1>>std::endl;  
    // ...
    cout << h_pttestD->GetNbinsX() << endl;
    cout << h_pttestSn->GetNbinsX() << endl;
    for (int i_q=1; i_q<=h_QtestD->GetNbinsX(); i_q++) {
        double x_q_Sn = h_QtestSn->GetXaxis()->GetBinCenter(i_q);		//no need maybe ? it will be always the same (I think /!\)
        double y_q_Sn = h_QtestSn->GetBinContent(i_q);				//hist_Q_Sn is the equivalent to  N^{A=Sn}_{h=pi}
        double x_q_De = h_QtestD->GetXaxis()->GetBinCenter(i_q);			//no need maybe ?
        double y_q_De = h_QtestD->GetBinContent(i_q) ;				//hist_Q_De is the equivalent to  N^{De}_{h=pi}
        double y_q_Sn_e = h_QtestSn->GetBinContent(i_q) ;			//hist_Q_Sn_e is the equivalent to  N^{A=Sn}_{h=pi}
        double y_q_De_e = h_QtestD->GetBinContent(i_q) ;
        double interm1 = y_q_Sn/y_q_Sn_e;
        double interm2 = y_q_De/y_q_De_e;
        double interm3 = interm1/interm2;
        double y_R_q = interm1;
        double Rq_err=  (1/y_q_Sn) + (1/y_q_De) + (1/y_q_Sn_e) + (1/y_q_De_e); 
		if(x_q_Sn>1.5){ 
        	    R_Q->SetPoint(i_q-1, x_q_Sn, interm3 );
		    R_Q->SetPointError(i_q-1, 0, Rq_err);
    		}
    }

    for (int i_v=1; i_v<=h_vtestD->GetNbinsX(); i_v++) {
        double x_v_Sn = h_vtestSn->GetXaxis()->GetBinCenter(i_v);		//no need maybe ? it will be always the same (I think /!\)
        double y_v_Sn = h_vtestSn->GetBinContent(i_v);				//hist_Q_Sn is the equivalent to  N^{A=Sn}_{h=pi}
        double x_v_De = h_vtestD->GetXaxis()->GetBinCenter(i_v);			//no need maybe ?
        double y_v_De = h_vtestD->GetBinContent(i_v) ;				//hist_Q_De is the equivalent to  N^{De}_{h=pi}
        double y_v_Sn_e = h_vtestSn->GetBinContent(i_v) ;			//hist_Q_Sn_e is the equivalent to  N^{A=Sn}_{h=pi}
        double y_v_De_e = h_vtestD->GetBinContent(i_v) ;
        double interm1 = y_v_Sn/y_v_Sn_e;
        double interm2 = y_v_De/y_v_De_e;
        double interm3 = interm1/interm2;
        double y_R_v = interm1;
        double Rv_err=  (1/y_v_Sn) + (1/y_v_De) + (1/y_v_Sn_e) + (1/y_v_De_e); 
		if(x_v_Sn>1.5){ 
        	    R_v->SetPoint(i_v-1, x_v_Sn, interm3 );
		    R_v->SetPointError(i_v-1, 0, Rv_err);
    		}
    }

    for (int i_z=1; i_z<=h_ztestD->GetNbinsX(); i_z++) {
        double x_z_Sn = h_ztestSn->GetXaxis()->GetBinCenter(i_z);		//no need maybe ? it will be always the same (I think /!\)
        double y_z_Sn = h_ztestSn->GetBinContent(i_z);				//hist_Q_Sn is the equivalent to  N^{A=Sn}_{h=pi}
        double x_z_De = h_ztestD->GetXaxis()->GetBinCenter(i_z);			//no need maybe ?
        double y_z_De = h_ztestD->GetBinContent(i_z) ;				//hist_Q_De is the equivalent to  N^{De}_{h=pi}
        double y_z_Sn_e = h_ztestSn->GetBinContent(i_z) ;			//hist_Q_Sn_e is the equivalent to  N^{A=Sn}_{h=pi}
        double y_z_De_e = h_ztestD->GetBinContent(i_z) ;
        double interm1 = y_z_Sn/y_z_Sn_e;
        double interm2 = y_z_De/y_z_De_e;
        double interm3 = interm1/interm2;
        double y_R_z = interm1;
        double Rz_err=  (1/y_z_Sn) + (1/y_z_De) + (1/y_z_Sn_e) + (1/y_z_De_e); 
		//if(x_z_Sn>1.5){                                               //not the sale condition
        	    R_z->SetPoint(i_z-1, x_z_Sn, interm3 );
		    R_z->SetPointError(i_z-1, 0, Rz_err);
    	//	}
    }
    for (int i_pt=1; i_pt<=h_pttestD->GetNbinsX(); i_pt++) {
        double x_pt_Sn = h_pttestSn->GetXaxis()->GetBinCenter(i_pt);		//no need maybe ? it will be always the same (I think /!\)
        double y_pt_Sn = h_pttestSn->GetBinContent(i_pt);				//hist_Q_Sn is the equivalent to  N^{A=Sn}_{h=pi}
        double x_pt_De = h_pttestD->GetXaxis()->GetBinCenter(i_pt);			//no need maybe ?
        double y_pt_De = h_pttestD->GetBinContent(i_pt) ;				//hist_Q_De is the equivalent to  N^{De}_{h=pi}
        double y_pt_Sn_e = h_pttestSn->GetBinContent(i_pt) ;			//hist_Q_Sn_e is the equivalent to  N^{A=Sn}_{h=pi}
        double y_pt_De_e = h_pttestD->GetBinContent(i_pt) ;
        double interm1 = y_pt_Sn/y_pt_Sn_e;
        double interm2 = y_pt_De/y_pt_De_e;
        double interm3 = interm1/interm2;
        double y_R_pt = interm1;
        double Rpt_err=  (1/y_pt_Sn) + (1/y_pt_De) + (1/y_pt_Sn_e) + (1/y_pt_De_e); 
		//if(x_pt_Sn>1.5){ 
        	    R_pt->SetPoint(i_pt-1, x_pt_Sn, interm3 );
		    R_pt->SetPointError(i_pt-1, 0, Rpt_err);
    		//}
    }

    // Close the first ROOT file
    fileD->Close();
    fileSn->Close();
    cR->cd(1);
    //R_Q->SetTitle(title);
    R_Q->GetXaxis()->SetTitle("Q{2} (GeV^{2}) " );
    R_Q->GetYaxis()->SetTitle("R^{#Pi+}_{Sn}");
    R_Q->Draw("AP");
    //cR->SaveAs("ratio.pdf");
    cR->cd(2);
    R_v->GetXaxis()->SetTitle("#nu " );
    R_v->GetYaxis()->SetTitle("R^{#Pi+}_{Sn}");
    R_v->Draw("AP");
    cR->cd(3);
    R_z->GetXaxis()->SetTitle("z " );
    R_z->GetYaxis()->SetTitle("R^{#Pi+}_{Sn}");
    R_z->Draw("AP");
    cR->cd(4);
    R_pt->GetXaxis()->SetTitle("pt " );
    R_pt->GetYaxis()->SetTitle("R^{#Pi+}_{Sn}");
    R_pt->Draw("AP");
    cR->SaveAs("ratio.pdf");
    

    return 0;
}

