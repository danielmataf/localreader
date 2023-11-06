// Include necessary ROOT headers
//Run and compile with 
//g++ -o multipratio tryRatio.cpp `root-config --cflags --libs`
// ./multipratio




////////////////////////////////////////////////////////////////
//Reusing old code for multipratio ... 
//It contains only a function that determines (maybe) the multipratio from a given histogram found in a root file 
// 
//  Oct 2023 
//Only plot R as a fct of pt2 for CLAS presentation.. 
//CleanCode is supposed to be applied to a set of data (cooked and MC) for LD2 and Sn
//Process MC data first with a given # evts, get the # of evts and process cooked with this #
// CleanCode will produce a root file everytime you run it. depending on the curts applied and the files in input
//Only change filenames and cuts  
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



void calculateMRat(int option, TH1F* hD, TH1F* h_onlyeD, int elecD, TH1F* hSn, TH1F* h_onlyeSn, int elecSn,  TGraphErrors* g_result, TGraphErrors* g_resulta) {
    int numBins = hD->GetNbinsX();

    for (int bin = 1; bin <= numBins; bin++) {
        double x_Sn = hSn->GetXaxis()->GetBinCenter(bin);		
        double y_Sn = hSn->GetBinContent(bin);				
        double x_D = hD->GetXaxis()->GetBinCenter(bin);			
        double y_D = hD->GetBinContent(bin) ;				
        double y_Sn_e = h_onlyeSn->GetBinContent(bin) ;			
        double y_D_e = h_onlyeD->GetBinContent(bin) ;
        double interm1;
        double interm2;
        double interm3;
        double Rq_err;
        if (option == 1){
            interm1 = (y_Sn_e > 0) ? y_Sn/y_Sn_e : 0.0;
            interm2 = (y_D_e > 0) ? y_D/y_D_e : 0.0;
            interm3 = (interm2 > 0) ? interm1/interm2 : 0.0;
            Rq_err= interm3 * sqrt(1/y_Sn + 1/y_D + 1/y_Sn_e + 1/y_D_e);
        }
        if (option==2) {
            interm1 = (elecSn > 0) ? y_Sn/elecSn : 0.0;
            interm2 = (elecD > 0) ? y_D/elecD : 0.0;
            interm3 = (interm2 > 0) ? interm1/interm2 : 0.0;
            cout<<interm3<<"->"<<endl;
            if (y_Sn!=0 && y_D!=0 && y_Sn_e != 0 && y_D_e != 0){
                Rq_err= interm3 * sqrt(1/y_Sn + 1/y_D + 1/y_Sn_e + 1/y_D_e);
                cout<<Rq_err<<"->"<<endl;
            }
        }
        g_result->SetPoint(bin-1, x_Sn, interm3 );
		g_result->SetPointError(bin-1, 0, Rq_err);
    }
}






int main() {
    TCanvas* cR=new TCanvas("ratio","ratio"); //create new canvas
    cR->Divide(2,2);
    TGraphErrors *R_pt = new TGraphErrors();
    TGraphErrors *R_pta = new TGraphErrors();
    TGraphErrors *ReRpt = new TGraphErrors();
    // Open ROOT files
    //change here to path of MC or cooked
    //  Choped has all histograms pre & pos 
    //  Use here Monitoring, it has only histograms pos cuts [and the electron counter]
    //TFile* fileD = new TFile("../files2read/REoutput_D.root", "READ");
    TFile* fileD = new TFile("/home/matamoros/Desktop/reboot/CleanCode/build/chop_ULTsimLD2.root","READ");
    //TFile* fileSn = new TFile("../files2read/REoutput_Sn.root", "READ");
    TFile* fileSn = new TFile("/home/matamoros/Desktop/reboot/CleanCode/build/chop_ULTsimSn.root","READ");
    TFile* fileDbis = new TFile("/home/matamoros/Desktop/reboot/CleanCode/build/output_ULTsimLD2.root","READ");
    TFile* fileSnbis = new TFile("/home/matamoros/Desktop/reboot/CleanCode/build/output_ULTsimSn.root","READ");

    //ERROR HANDLER - file open
    if (!fileD->IsOpen()) {
        std::cerr << "Error: Cant open fileLD2" << std::endl;
        return 0 ;
    }
    if (!fileSn->IsOpen()) {
        std::cerr << "Error: Cant open fileSn" << std::endl;
        return 0 ;
    }
    if (!fileDbis->IsOpen()) {
        std::cerr << "Error: Cant open fileLD2bis" << std::endl;
        return 0 ;
    }
    if (!fileSnbis->IsOpen()) {
        std::cerr << "Error: Cant open fileSnbis" << std::endl;
        return 0 ;
    }



    // Retrieve the first histogram from the ROOT file
    //  USING MONITORING histNAMES 
    //TH1F* h_pttestD = dynamic_cast<TH1F*>(fileD->Get("pt2p"));
    TH1F* h_pttestD = dynamic_cast<TH1F*>(fileD->Get("pt2R"));
    TH1F* h_pttestSn = dynamic_cast<TH1F*>(fileSn->Get("pt2R"));
    TH1F* h_QtestD = dynamic_cast<TH1F*>(fileD->Get("Q2R"));
    TH1F* h_QtestSn = dynamic_cast<TH1F*>(fileSn->Get("Q2R"));
    TH1F* h_nutestD = dynamic_cast<TH1F*>(fileD->Get("nuR"));
    TH1F* h_nutestSn = dynamic_cast<TH1F*>(fileSn->Get("nuR"));
    TH1F* h_ztestD = dynamic_cast<TH1F*>(fileD->Get("zR"));
    TH1F* h_ztestSn = dynamic_cast<TH1F*>(fileSn->Get("zR"));
    TH1F* h_QtestDpos = dynamic_cast<TH1F*>(fileDbis->Get("Q2truePOS"));
    TH1F* h_QtestSnpos = dynamic_cast<TH1F*>(fileSnbis->Get("Q2truePOS"));
    TH1F* h_nutestDpos = dynamic_cast<TH1F*>(fileDbis->Get("nutruePOS"));
    TH1F* h_nutestSnpos = dynamic_cast<TH1F*>(fileSnbis->Get("nutruePOS"));


    if (!h_pttestD || !h_pttestSn) {
            std::cerr << "Error: Cant retrieve histograms" << std::endl;
            return 0 ;
        }
    
    //TH1F* h_ptonlyeD = dynamic_cast<TH1F*>(fileDbis->Get("onlye_pt"));
    //TH1F* h_ptonlyeSn = dynamic_cast<TH1F*>(fileSnbis->Get("onlye_pt"));
    //if (!h_ptonlyeD || !h_ptonlyeSn) {
    //        std::cerr << "Error: Cant retrieve histogram singles" << std::endl;
    //        return 0 ;
    //    }

        
    TTree* treeD = dynamic_cast<TTree*>(  fileDbis->Get("treecounter"));
    TTree* treeSn = dynamic_cast<TTree*>(fileSnbis->Get("treecounter"));
        if (!treeD || !treeSn) {
            std::cerr << "Error: Cant retrieve trees" << std::endl;
            return 0 ;
        }
    

    
    
    



    int counter_e_D; 
    
    treeD->SetBranchAddress("counterel_R", &counter_e_D);
    treeD->GetEntry(0);
    int counter_e_Sn;
    treeSn->SetBranchAddress("counterel_R", &counter_e_Sn);
    treeSn->GetEntry(0);
    double intermD = counter_e_D;
    double intermSn= counter_e_Sn;
    
    calculateMRat( 2, h_pttestD, h_pttestD,  intermD, h_pttestD, h_pttestD,  intermSn, R_pt,R_pta);



    //Close the first ROOT file
    fileD->Close();
    fileSn->Close();
    cR->cd(1);
    R_pta->GetXaxis()->SetTitle("Q{2} (GeV^{2}) " );
    R_pta->GetYaxis()->SetTitle("R^{#Pi+}_{Sn}");
    R_pta->Draw("AP");
    ////cR->SaveAs("ratio.pdf");
    fileSn->Close();
    cR->cd(2);
    R_pt->GetXaxis()->SetTitle("Q{2} (GeV^{2}) " );
    R_pt->GetYaxis()->SetTitle("R^{#Pi+}_{Sn}");
    R_pt->Draw("AP");
    cR->SaveAs("ratiotest.pdf");
    

    return 0;
}


    /*
    for (int i_pt=1; i_pt<=h_pttestD->GetNbinsX(); i_pt++) {
        double x_pt_Sn = h_pttestSn->GetXaxis()->GetBinCenter(i_pt);		//no need maybe ? it will be always the same (I think /!\)
        float y_pt_Sn = h_pttestSn->GetBinContent(i_pt);				//hist_Q_Sn is the equivalent to  N^{A=Sn}_{h=pi}
        double x_pt_De = h_pttestD->GetXaxis()->GetBinCenter(i_pt);			//no need maybe ?
        float y_pt_De = h_pttestD->GetBinContent(i_pt) ;				//hist_Q_De is the equivalent to  N^{De}_{h=pi}
        double y_pt_Sn_e = intermSn; //h_ptonlyeD->GetBinContent(i_pt) ;			//hist_Q_Sn_e is the equivalent to  N^{A=Sn}_{h=pi}
        double y_pt_De_e = intermD; //h_ptonlyeD->GetBinContent(i_pt) ;
        double interm1 = y_pt_Sn/y_pt_Sn_e;
        double interm2 = y_pt_De/y_pt_De_e;
        double interm3 = interm1/interm2;
        double y_R_pt = interm1;
        double Rpt_err= sqrt( pow((sqrt(y_pt_Sn)/y_pt_Sn),2) + pow((sqrt(y_pt_De)/y_pt_De),2) + pow((sqrt(y_pt_Sn_e)/y_pt_Sn_e),2) + pow((sqrt(y_pt_De_e)/y_pt_De_e),2) );
        //double Rpt_err=  (1/y_pt_Sn) + (1/y_pt_De) + (1/y_pt_Sn_e) + (1/y_pt_De_e);     // false uncertainty 
        //cout<< "y1/y_e1 = "<< interm1<<endl;
        //cout<< "y2/y_e2 = "<< interm2<<endl;
        cout<< "y_pt_Sn = "<< y_pt_Sn<<endl;
        cout<< "y_pt_De = "<< y_pt_De <<endl;
        cout<< "y_pt_Sn_e = "<< y_pt_Sn_e <<endl;
        cout<< "y_pt_De_e = "<< y_pt_De_e <<endl;
        cout<< "R error= "<<Rpt_err<<endl;
        cout<< "...................."<<endl;
		if(interm2!=0 &&  x_pt_Sn<1.2 ){                //added condition on the pt value (x axis) bc pt histograms dont go further than these value  
        	    ReRpt->SetPoint(i_pt-1, x_pt_Sn, interm3 ); //
		        ReRpt->SetPointError(i_pt-1, 0,Rpt_err );
    	}
    }
*/
