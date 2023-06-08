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
    TCanvas* cR=new TCanvas("ratio","ratio"); //create new canvas
    TGraphErrors *R_Q = new TGraphErrors();
    // Open the first ROOT file
    TFile* file1 = new TFile("build/output1.root", "READ");

    // Retrieve the first histogram from the ROOT file
    TH1F* h_Q1test = dynamic_cast<TH1F*>(file1->Get("Q2"));

    // Access properties of h_Q1
    //std::cout>>numBins1>>std::endl;  
    // ...
    cout << h_Q1test->GetNbinsX() << endl;
    for (int i_q=1; i_q<=h_Q1test->GetNbinsX(); i_q++) {
        double x_q_Sn = h_Q1test->GetXaxis()->GetBinCenter(i_q);		//no need maybe ? it will be always the same (I think /!\)
        double y_q_Sn = h_Q1test->GetBinContent(i_q);				//hist_Q_Sn is the equivalent to  N^{A=Sn}_{h=pi}
        double x_q_De = h_Q1test->GetXaxis()->GetBinCenter(i_q);			//no need maybe ?
        double y_q_De = h_Q1test->GetBinContent(i_q) ;				//hist_Q_De is the equivalent to  N^{De}_{h=pi}
        double y_q_Sn_e = h_Q1test->GetBinContent(i_q) ;			//hist_Q_Sn_e is the equivalent to  N^{A=Sn}_{h=pi}
        double y_q_De_e = h_Q1test->GetBinContent(i_q) ;
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

