// Include necessary ROOT headers
//Run and compile with 
//g++ -o multipratio MultipRatio.cpp `root-config --cflags --libs`
// ./multipratio


////////////////////////////////////////////////////////////////
//Code separated from the first one. 
//It contains only a function that determines (maybe) the multipratio from a given histogram found in a root file 
// May 2023
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
    TCanvas* cDpt=new TCanvas("Dptplot","Dptplot"); //create new canvas
    TGraphErrors *testDptQ = new TGraphErrors();
    // Open the first ROOT file
    TFile* file1 = new TFile("build/output1.root", "READ");

    // Retrieve the first histogram from the ROOT file
    //TH1F* h_Q1test = dynamic_cast<TH1F*>(file1->Get("Q2"));

    // Access properties of h_Q1
    //std::cout>>numBins1>>std::endl;  
    // ...
    cout << h_Q1test->GetNbinsX() << endl;
    /*
    for (int i_q=1; i_q<=w_ptQ_Sn->GetNbinsX(); i_q++) {
        double x_q_w = w_ptQ_Sn->GetXaxis()->GetBinCenter(i_q);		
        double y_q_wSn = w_ptQ_Sn->GetBinContent(i_q);				
        double x_q = hist_Q_Sn->GetXaxis()->GetBinCenter(i_q);			
        double y_qSn = hist_P_t_Sn->GetBinContent(i_q) ;
	double y_q_wDe = w_ptQ_De->GetBinContent(i_q);
	double y_qDe = hist_P_t->GetBinContent(i_q) ;
	double meandptQSn =(y_qSn > 0) ?  y_q_wSn / y_qSn : 0.0;
	double meandptQDe = ( y_qDe > 0) ?  y_q_wDe / y_qDe : 0.0;
	//cout<<meandptQSn - meandptQDe <<"meanDPTQ"<<endl;
	//double Rq_err=  (1/y_q_Sn) + (1/y_q_De) + (1/y_q_Sn_e) + (1/y_q_De_e); 
		if(x_q_w>1.5){ 
        	    testDptQ->SetPoint(i_q-1, x_q_w,meandptQSn - meandptQDe);
    		}
	//}
    }
    */

    // Close the first ROOT file
    file1->Close();
    cR->cd(1);
    //R_Q->SetTitle(title);
    R_Q->GetXaxis()->SetTitle("Q{2} (GeV^{2}) " );
    R_Q->GetYaxis()->SetTitle("R^{#Pi+}_{Sn}");
    R_Q->Draw("AP");
    cR->SaveAs("ratio.pdf");

    return 0;
}

