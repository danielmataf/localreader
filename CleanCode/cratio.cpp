#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <fstream>
#include <math.h>
#include <TPDF.h>
#include "Event.h" 
#include "cratio.h"
#include "CutSet.h"


cratio::cratio(CutSet cutsD, CutSet cutsA): //: cutsD(cutsD), cutsSn(cutsSn) {
    CratioMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    errorCratioMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))) {

    cutd = cutsD;
    cuta = cutsA;
    }



    void cratio::FillHistograms(const Event& event){
        int targetType = event.GetTargetType();
        if (targetType == 0 && cutd.PassCutsElectrons(event)==true) {
            counter_elLD2 ++;
            //set a counter that increases when electroncuts = passed; in order for it to be called when R is  computed in had variables (?) TBD
            hCratio_nuD->Fill(event.Getnu());
            for (const Particle& hadron : event.GetHadrons()) {
                if (cutd.PassCutsHadrons(hadron)==true){
                    double phiD = hadron.Getphih();
                    h_wD_Cratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), cos(phiD));    //3 arguments and the WEIGHT
                    h_wD_sqCratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), cos(phiD)*cos(phiD));    //3 arguments and the WEIGHT (pt2 squared) 4 variance
                    h_D_Cratio3D->Fill(event.Getxb(), event.GetQ2(), hadron.Getz());    //3 arguments only counts not weight (cphi)
                }
            }
        }
        else if (targetType == 1 && cuta.PassCutsElectrons(event)==true) {
            counter_elSn++; //counter for only electrons for z and pt
            //here change the else if to just else in order to have a generic target 
            hCratio_nuA->Fill(event.Getnu());
            for (const Particle& hadron : event.GetHadrons()) {
                if (cuta.PassCutsHadrons(hadron)==true){
                    double phiA = hadron.Getphih();
                    h_wA_Cratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), cos(phiA));    //3 arguments and the WEIGHT
                    h_wA_sqCratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), cos(phiA)*cos(phiA));    //3 arguments and the WEIGHT (pt2 squared)
                    h_A_Cratio3D->Fill(event.Getxb(), event.GetQ2(), hadron.Getz());    //3 arguments only counts not weight

                }
            }
        }
    }


void cratio::calcCratio(){
    int numBinsX = h_wD_Cratio->GetNbinsX();    //same bins 4 Deut and A 
    int numBinsY = h_wD_Cratio->GetNbinsY(); 
    int numBinsZ = h_wD_Cratio->GetNbinsZ(); 
    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {
        for (int Ybin=1; Ybin <= numBinsY; Ybin++) {
            for (int Zbin= 1; Zbin <= numBinsZ; Zbin++){
                double valD = h_wD_Cratio->GetBinContent(Xbin,Ybin,Zbin);
                double valA = h_wA_Cratio->GetBinContent(Xbin,Ybin,Zbin);
                double countD = h_D_Cratio3D->GetBinContent(Xbin,Ybin,Zbin);    //3D histo no counts
                double countA = h_A_Cratio3D->GetBinContent(Xbin,Ybin,Zbin);
                double sqvalD = h_wD_sqCratio->GetBinContent(Xbin,Ybin,Zbin);
                double sqvalA = h_wA_sqCratio->GetBinContent(Xbin,Ybin,Zbin);
                double wavg_CratioD = (countD > 0) ? valD/countD : 0.0;   //weighted average ok 
                double wavg_CratioA = (countA > 0) ? valA/countA : 0.0;
                double Cratio_point = wavg_CratioA / wavg_CratioD;

                double varianceD = (countD > 0) ? sqvalD/countD - wavg_CratioD*wavg_CratioD : 0.0;
                double varianceA = (countA > 0) ? sqvalA/countA - wavg_CratioA*wavg_CratioA : 0.0;
                double Err_valD = (countD > 0) ? sqrt(varianceD/countD) : 0.0;
                double Err_valA = (countA > 0) ? sqrt(varianceA/countA) : 0.0;
                //double Err_Cratio_point = Cratio_point * sqrt(Err_valD*Err_valD + Err_valA*Err_valA);
                double Err_Cratio_point = Cratio_point * sqrt(pow(Err_valA / wavg_CratioA, 2) + pow(Err_valD / wavg_CratioD, 2));

//                double Err_dpt_point = sqrt(Err_valD*Err_valD + Err_valA*Err_valA);
                CratioMatrix[Xbin-1][Ybin-1][Zbin-1] = Cratio_point;
                errorCratioMatrix[Xbin-1][Ybin-1][Zbin-1] = Err_Cratio_point;

            }
        }
    }  



}




void cratio::writeMatrixToFile(const std::string& filename){
    //not needed unless output in file wanted 
    std::ofstream outputFile(filename);
    for (int x = 0; x < Cratiobin; ++x) {
        double xaxisval = h_wD_Cratio->GetXaxis()->GetBinCenter(x + 1);
        outputFile << "xb = " << xaxisval << std::endl;
        for (int y =0; y< Cratiobin; ++y) {
            for (int z = 0; z < Cratiobin; ++z) {
                double value = CratioMatrix[x][y][z];
                double error = errorCratioMatrix[x][y][z];
                outputFile << value << " +- " << error << "\t\t";
            }
            outputFile << std::endl;
        }
        outputFile << std::endl << std::endl;
    }
    outputFile.close();
}
    



void cratio::multiplotCratio(){
    for (int x = 0; x < Cratiobin; ++x){
        double xValue = h_wD_Cratio->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "multiplotCratio_x" + std::to_string(xValue) + ".pdf";
        TCanvas canvasCratio("c", "Multiplot Cratio", 1200, 800);
        canvasCratio.Divide(3, 2); 
        for (int y=0; y < Cratiobin; ++y) {
            double Q2Value = h_wD_Cratio->GetYaxis()->GetBinCenter(y + 1);
            canvasCratio.cd(y + 1);
            TGraphErrors *graphCratio = new TGraphErrors();
            for (int z=0; z < Cratiobin; ++z) {
                double zValue = h_wD_Cratio->GetZaxis()->GetBinCenter(z + 1);
                double value = CratioMatrix[x][y][z];
                double error = errorCratioMatrix[x][y][z];
                graphCratio->SetPoint(y, zValue, value);
                graphCratio->SetPointError(y, 0.0, error); 
                
            }
            graphCratio->SetTitle(("<cos #phi_{h}>_A / <cos #phi_{h}>_D vs z, Q^{2}=" + std::to_string(Q2Value)).c_str());
            graphCratio->GetXaxis()->SetTitle("z");
            graphCratio->GetYaxis()->SetTitle("<cos #phi_{h}>_A / <cos #phi_{h}>_D");
            graphCratio->SetMarkerStyle(20);
            graphCratio->Draw("AP");

        }
        canvasCratio.SaveAs(pdfFileName.c_str());
 
    }

}
