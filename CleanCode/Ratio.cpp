#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
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
    if (target == "D" && cutd.PassCutsElectrons(event)==true) {
        h_nu_D->Fill(event.Getnu());
        for (const Particle& hadron : event.GetHadrons()) {
            if (cutd.PassCutsHadrons(hadron)==true){
                //not using the if (==false) return statement 
                h_z_D->Fill(hadron.Getz());
                h_pt2_D->Fill(hadron.Getpt2());
            }
        }
    }
    else if (target == "Sn" && cuta.PassCutsElectrons(event)==true) {
        h_nu_A->Fill(event.Getnu());
        for (const Particle& hadron : event.GetHadrons()) {
            if (cuta.PassCutsHadrons(hadron)==true){
                h_z_A->Fill(hadron.Getz());
                h_pt2_A->Fill(hadron.Getpt2());
            }
        }
    }
}


void Ratio::calcR(){
    //all histos are supposed to be filled.
    //no need to include them as argument of the function just recover themm by calling em 
    //maybe need of two different fillers for the h_nu histo (only electron)
    int numBins = h_nu_D->GetNbinsX();  //all bins should be the same 
                                        // DONT FORGET TO DEFINE THE BINNING OF KINVAR TO Rbinning = 10
    for (int bin = 1; bin <= numBins; bin++) {

        double nu_D  = h_nu_D->GetBinContent(bin) ;				
        double z_D   = h_z_D->GetBinContent(bin) ;				
        double pt2_D = h_pt2_D->GetBinContent(bin) ;				
        double nu_A  = h_nu_A->GetBinContent(bin) ;				
        double z_A   = h_z_A->GetBinContent(bin) ;				
        double pt2_A = h_pt2_A->GetBinContent(bin) ;				
    }   

}



//void Ratio::WriteHistos(const std::string filename) {
//    outputFile->cd();
//    h_nu_D->Write();
//    h_z_D->Write();
//    h_pt2_D->Write();
//    h_nu_Sn->Write();
//    h_z_Sn->Write();
//    h_pt2_Sn->Write();
//    outputFile->Write();
//    outputFile->Close();
//}

    //outputFile = new TFile("ratio_output.root", "RECREATE");
