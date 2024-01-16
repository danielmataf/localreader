#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <fstream>
#include <TPDF.h>
#include "Event.h" 
#include "Ratio.h"
#include "CutSet.h"

Ratio::Ratio(CutSet cutsD, CutSet cutsA): //: cutsD(cutsD), cutsSn(cutsSn) {
    ratMatrix(Rbin, std::vector<std::vector<double>>(Rbin, std::vector<double>(Rbin, 0.0))),
    errorMatrix(Rbin, std::vector<std::vector<double>>(Rbin, std::vector<double>(Rbin, 0.0))) {
//Ratio::Ratio(CutSet a)
  
    //cut1 = a;
    cutd = cutsD;
    cuta = cutsA;
    
    //class takes a set of cuts anyway 
    //in order to determine Ratio from the histos with these cuts 
}

//~ add deletes 


void Ratio::FillHistograms(const Event& event, const std::string target) {
        
    int targetType = event.GetTargetType();
        
    //add passcut el onlyu once ten cut on target. inverse to =false
    //then add  ALu o cut needed + pass cut hadrons 


                                                //using a flag for targets 
    //forget condition on taarget, that should directly be a cut !!!! TBD 
    if (target == "D" && cutd.PassCutsElectrons(event)==true) {
        counter_elLD2 ++;
        //set a counter that increases when electroncuts = passed; in order for it to be called when R is  computed in had variables (?)
        h_nuD->Fill(event.Getnu());
        for (const Particle& hadron : event.GetHadrons()) {
            if (cutd.PassCutsHadrons(hadron)==true){
                ////not using the if (==false) return statement 

                h_nu_z_pt2D->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
            }
        }
    }
    else if (target == "Sn" && cuta.PassCutsElectrons(event)==true) {
        counter_elSn++; //counter for only electrons for z and pt
        //here change the else if to just else in order to have a generic target 
        h_nuA->Fill(event.Getnu());
        for (const Particle& hadron : event.GetHadrons()) {
            if (cuta.PassCutsHadrons(hadron)==true){
                h_nu_z_pt2A->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
            }
        }
        //Add Here Cu && change cut for Sn here. 1st cut is passelectrons 
    }

    //Add CxC

}






void Ratio::calcR(){
    //histos are supposed to be filled.
    //no need to include them as argument of the function just recover themm by calling em 
    //maybe need of two different fillers for the h_nu histo (only electron)
    

    int graph_pointnb = 0;
    //IN 3d HISTO x=nu; y= z; z=pt2;
    int counter_3D = 0;
    int numBinsX = h_nu_z_pt2D->GetNbinsX();    //same bins 4 Deut and A 
    int numBinsY = h_nu_z_pt2D->GetNbinsY(); 
    int numBinsZ = h_nu_z_pt2D->GetNbinsZ(); 
    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {  
        double val_nuelD = h_nuD->GetBinContent(Xbin);   
        double val_nuelA = h_nuA->GetBinContent(Xbin); 
          
        for (int Ybin = 1; Ybin <= numBinsY; Ybin++ ){
            for (int Zbin = 1; Zbin <= numBinsZ; Zbin++ ){
                // loop with 125 values (5 per axis)
                double valD = h_nu_z_pt2D->GetBinContent(Xbin,Ybin,Zbin);   //seems to properly recover value 
                double valA = h_nu_z_pt2A->GetBinContent(Xbin,Ybin,Zbin);
                double interm1nu = (val_nuelA > 0) ? valA / val_nuelA : 0.0;
                double interm2nu = (val_nuelD > 0) ? valD / val_nuelD : 0.0;
                double ratvalue  = (interm2nu > 0) ?interm1nu/interm2nu : 0.0;
                double raterr = ratvalue * sqrt(1/valA + 1/valD + 1/val_nuelA + 1/val_nuelD);
                ratMatrix[Xbin - 1][Ybin - 1][Zbin - 1] = ratvalue;
                errorMatrix[Xbin - 1][Ybin - 1][Zbin - 1] = raterr;
            }
        }

    //    Rq_err= interm3 * sqrt(1/nu_A_had + 1/nu_D_had + 1/nu_A_had + 1/nu_D_had);
    //    R_v->SetPoint(bin-1, x_axis, interm3 );
	//	//	
    }   
}

//TH3 getbins bini binj etc 

void Ratio::writeMatrixToFile(const std::string& filename) {
    std::ofstream outputFile(filename);
    //IN 3d HISTO x=nu; y= z; z=pt2;
    for (int x = 0; x < Rbin; ++x) {
        double xaxisval = h_nu_z_pt2A->GetXaxis()->GetBinCenter(x + 1);
        outputFile << "nu = " << xaxisval << std::endl;
        for (int z = 0; z < Rbin; ++z) {
            for (int y = 0; y < Rbin; ++y) {
                double value = ratMatrix[x][y][z];
                double error = errorMatrix[x][y][z];
                //outputFile << "rat[" << x << "][" << y << "][" << z << "] +- err[" << x << "][" << y << "][" << z << "] = "
                //           << value << " +- " << error << "\t\t";
                outputFile << value << " +- " << error << "\t\t";
            
            }
            outputFile << std::endl;
        }
        outputFile << std::endl << std::endl;
    }
    outputFile.close();
}


void Ratio::multiplotR() {
    for (int x = 0; x < Rbin; ++x) {

        double nuValue = h_nu_z_pt2D->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "multiplotR_nu" + std::to_string(nuValue) + ".pdf";
        TCanvas canvas("c", "Multiplot R", 1200, 800);
        canvas.Divide(3, 2); 
        for (int z = 0; z < Rbin; ++z) {
            double pt2Value = h_nu_z_pt2D->GetZaxis()->GetBinCenter(z + 1);
            canvas.cd(z + 1);
            TGraphErrors *graph = new TGraphErrors();
            for (int y = 0; y < Rbin; ++y) {
                double zValue = h_nu_z_pt2D->GetYaxis()->GetBinCenter(y + 1);
                double value = ratMatrix[x][y][z];
                double error = errorMatrix[x][y][z];
                graph->SetPoint(y, zValue, value);
                graph->SetPointError(y, 0.0, error); 
            }
            graph->SetTitle(("R vs z, pt2=" + std::to_string(pt2Value)).c_str());
            graph->GetXaxis()->SetTitle("z");
            graph->GetYaxis()->SetTitle("R");
            graph->SetMarkerStyle(20);
            graph->Draw("AP");
        }

        canvas.SaveAs(pdfFileName.c_str());
    }
}

void Ratio::PlotRatio(const std::string filename) {
    TCanvas Rcanv("Ratio canvas", "Ratio Plots");
    Rcanv.Divide(1, 2);
    Rcanv.cd(1);
    h_nu_z_pt2A->Draw();
    Rcanv.cd(2);
    graph_rat->Draw();
    Rcanv.Print((filename + ".pdf").c_str());

}




    //outputFile = new TFile("ratio_output.root", "RECREATE");
