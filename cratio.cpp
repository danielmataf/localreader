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
 
	//----Dpt @ Q try2---//		//WORKING!
int nbinsx = h_ptQ_Sn->GetNbinsX();
double* avg = new double[nbinsx];
double* err = new double[nbinsx];
for (int xbin = 1; xbin <= nbinsx; xbin++) {
  double sum = 0.0;
  double sumsq = 0.0;
  int count = 0;
  for (int ybin = 1; ybin <= h_ptQ_Sn->GetNbinsY(); ybin++) {
    double content = h_ptQ_Sn->GetBinContent(xbin, ybin);
    double yvalue = h_ptQ_Sn->GetYaxis()->GetBinCenter(ybin);
    sum += content * yvalue;
    sumsq += content * yvalue * yvalue;	//somme des carrés, weighted
    count += content;			// weighted average divise (la somme des valeurs*poids) par la somme de poids, et non par la somme d'events 
  }
  avg[xbin-1] = (count > 0) ? sum / count : 0.0;
  double variance = (count > 1) ? (sumsq - sum*sum/count) / (count - 1) : 0.0;
  err[xbin-1] = (count > 1) ? sqrt(variance/count) : 0.0;
  cout<<"variance= "<<variance<<endl;
}

int nbinsx2 = h_ptQ_De->GetNbinsX();
double* avg2 = new double[nbinsx2];
double* err2 = new double[nbinsx2];
for (int xbin2 = 1; xbin2 <= nbinsx2; xbin2++) {
  double sum2 = 0.0;
  double sumsq2 = 0.0;
  int count2 = 0;
  for (int ybin2 = 1; ybin2 <= h_ptQ_De->GetNbinsY(); ybin2++) {
    double content2 = h_ptQ_De->GetBinContent(xbin2, ybin2);
    double yvalue2 = h_ptQ_De->GetYaxis()->GetBinCenter(ybin2);
    sum2 += content2 * yvalue2;
    sumsq2 += content2 * yvalue2 * yvalue2;
    count2 += content2;
  }
  avg2[xbin2-1] = (count2 > 0) ? sum2 / count2 : 0.0;
  double variance2 = (count2 > 1) ? (sumsq2 - sum2*sum2/count2) / (count2 - 1) : 0.0;
  err2[xbin2-1] = (count2 > 1) ? sqrt(variance2/count2) : 0.0;
}

double* avgDptQ = new double[nbinsx2];
double* errDptQ = new double[nbinsx2];
for (int i_avgQ = 1; i_avgQ <= nbinsx2-1; i_avgQ++) {
    double x_Dpt = hist_Q_Sn->GetXaxis()->GetBinCenter(i_avgQ);
    avgDptQ[i_avgQ]  = avg[i_avgQ] - avg2[i_avgQ];
    errDptQ[i_avgQ] = sqrt(err[i_avgQ]*err[i_avgQ] + err2[i_avgQ]*err2[i_avgQ]);
    //cout<<"err1 = "<<err[i_avgQ-1]<< "and err 2= "<<err2[i_avgQ-1]<<endl;
    if(x_Dpt>1.5){ 
	plotDpt_Q->SetPoint(i_avgQ-1,x_Dpt,avgDptQ[i_avgQ-1]);
	plotDpt_Q->SetPointError(i_avgQ-1, 0, errDptQ[i_avgQ-1]);
	//cout<<errDptQ[i_avgQ]<<endl;
	//plotDpt_Q->SetPoint(i_avgQ,x_Dpt,avgDptQ[i_avgQ]);
	//plotDpt_Q->SetPointError(i_avgQ, 0, errDptQ[i_avgQ]);
    }
}



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
cos_Q->Draw("AP");
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
cos_z->GetXaxis()->SetTitle("p_{t}^{2} " );
cos_z->GetYaxis()->SetTitle("<cos(#Phi)>_{Sn} / <cos(#Phi)>_{De}");
cos_z->Draw("AP");
cR->cd(4);
cos_pt->GetXaxis()->SetTitleSize(0.05);
cos_pt->GetYaxis()->SetTitleSize(0.05);
cos_pt->SetMarkerSize(0.5);
cos_pt->SetMarkerStyle(21);
cos_pt->GetXaxis()->SetRangeUser(0.2,0.9);
cos_pt->GetXaxis()->SetTitle("z " );
cos_pt->GetYaxis()->SetTitle("<cos(#Phi)>_{Sn} / <cos(#Phi)>_{De}");
cos_pt->Draw("AP");
    cR->SaveAs("cratio.pdf");
    cR->SaveAs("cratio.root");

    return 0;
}

