#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPDF.h>
#include "Event.h" 
#include "Ratio.h"
#include "CutSet.h"

Ratio::Ratio(CutSet cutsD, CutSet cutsA) //: cutsD(cutsD), cutsSn(cutsSn) {
//Ratio::Ratio(CutSet a)
  {
    //cut1 = a;
    cutd = cutsD;
    cuta = cutsA;
    //class takes a set of cuts anyway 
    //in order to determine Ratio from the histos with these cuts 
}

//~ add deletes 


void Ratio::FillHistograms(const Event& event, const std::string target) {
                                                //using a flag for targets 
    //forget condition on taarget, that should directly be a cut !!!! TBD 
    if (target == "D" && cutd.PassCutsElectrons(event)==true) {


        
        //h_nu_D->Fill(event.Getnu());
        //set a counter that increases when electroncuts = passed; in order for it to be called when R is  computed in had variables (?)
        
        for (const Particle& hadron : event.GetHadrons()) {
            if (cutd.PassCutsHadrons(hadron)==true){
                
                //h_nu_D_had->Fill(event.Getnu());
                ////not using the if (==false) return statement 
                //h_z_D->Fill(hadron.Getz());
                //h_pt2_D->Fill(hadron.Getpt2());
                h_nu_z_pt2D->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
                
                
                //!!!!!!!!!!    
                //define a hist TH3 w/ nu, z  et pt2
                //!!!!!!!!!!    
            
            
            }
        }
    }
    else if (target == "Sn" && cuta.PassCutsElectrons(event)==true) {
        //here change the else if to just else in order to have a generic target 
        //h_nu_A->Fill(event.Getnu());
        for (const Particle& hadron : event.GetHadrons()) {
            if (cuta.PassCutsHadrons(hadron)==true){
                //h_nu_A_had->Fill(event.Getnu());
                //h_z_A->Fill(hadron.Getz());
                //h_pt2_A->Fill(hadron.Getpt2());
                h_nu_z_pt2A->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
            }
        }
    }
}






void Ratio::calcR(){
    //all histos are supposed to be filled.
    //no need to include them as argument of the function just recover themm by calling em 
    //maybe need of two different fillers for the h_nu histo (only electron)
    


    //create a general function that determines R for all variables (?) 
        //then call R with the specific variable as an argument  argument 
        //and use switch in order to have a generic histogram ? 

    //IN 3d HISTO x=nu; y= z; z=pt2;
    int counter_3D = 0;
    int numBinsX = h_nu_z_pt2D->GetNbinsX();    //same bins 4 Deut and A 
    int numBinsY = h_nu_z_pt2D->GetNbinsY(); 
    int numBinsZ = h_nu_z_pt2D->GetNbinsZ(); 
    //TGraphErrors *R_v = new TGraphErrors();
    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {  
        double nu_D  = h_nu_z_pt2D->GetBinContent(Xbin,0,0) ;	
        double nu_D_had = h_nu_z_pt2A->GetBinContent(Xbin,0,0) ;
        for (int Ybin = 1; Ybin <= numBinsY; Ybin++ ){
            double nu_D_y  =    h_nu_z_pt2D->GetBinContent(Xbin,Ybin,0) ;	
            double nu_D_had_y = h_nu_z_pt2A->GetBinContent(Xbin,Ybin,0) ;
            for (int Ybin = 1; Ybin <= numBinsY; Ybin++ ){

                counter_3D ++;  //1000 elements here
            }
        }
      
    //    double x_axis =h_nu_D->GetXaxis()->GetBinCenter(bin);
    //            //should be the same for all 		
    //    double nu_D  = h_nu_D->GetBinContent(bin) ;	
    //    double nu_D_had = h_nu_D_had->GetBinContent(bin) ;
    //    double z_D   = h_z_D->GetBinContent(bin) ;				
    //    double pt2_D = h_pt2_D->GetBinContent(bin) ;				
    //    double nu_A  = h_nu_A->GetBinContent(bin) ;				
    //    double nu_A_had = h_nu_D_had->GetBinContent(bin) ;
    //    double z_A   = h_z_A->GetBinContent(bin) ;				
    //    double pt2_A = h_pt2_A->GetBinContent(bin) ;	
    //    double Rq_err;
    //    double interm1 = (nu_A > 0) ? nu_A_had/nu_A : 0.0;
    //    double interm2 = (nu_D > 0) ? nu_D_had/nu_D : 0.0;
    //    double interm3 = (interm2 > 0) ? interm1/interm2 : 0.0;
    //    Rq_err= interm3 * sqrt(1/nu_A_had + 1/nu_D_had + 1/nu_A_had + 1/nu_D_had);
    //    R_v->SetPoint(bin-1, x_axis, interm3 );
	//	//	
    }   
}

//TH3 getbins bini binj etc 

void Ratio::PlotRatio(const std::string filename) {
    TCanvas Rcanv("Ratio canvas", "Ratio Plots");
    Rcanv.Divide(1, 2);
    Rcanv.cd(1);
    h_nu_z_pt2A->Draw("lego");
    Rcanv.cd(2);
    h_nu_z_pt2D->Draw("lego");
    Rcanv.Print((filename + ".pdf").c_str());

}

    //outputFile = new TFile("ratio_output.root", "RECREATE");
