#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <fstream>
#include <TPDF.h>
#include "Event.h" 
#include "Dpt.h"
#include "CutSet.h"


deltaptsq::deltaptsq(CutSet cutsD, CutSet cutsA): //: cutsD(cutsD), cutsSn(cutsSn) {
    DptMatrix(Dptbin, std::vector<std::vector<double>>(Dptbin, std::vector<double>(Dptbin, 0.0))),
    errorDptMatrix(Dptbin, std::vector<std::vector<double>>(Dptbin, std::vector<double>(Dptbin, 0.0))) {

    cutd = cutsD;
    cuta = cutsA;
    }

// ONLY NEED 3 VARIABLES Q,  Z and nu
// Dont use pt nor pt2 as a variable 

void deltaptsq::FillHistograms(const Event& event) {
    int targetType = event.GetTargetType();
    if (targetType == 0 && cutd.PassCutsElectrons(event)==true) {
        counter_elLD2 ++;
        hDpt_nuD->Fill(event.Getnu());
        //still using nu for reference following Mratio 
        //Q can be used ? (needs to be checked)
        for (const Particle& hadron : event.GetHadrons()) {
            if (cutd.PassCutsHadrons(hadron)==true){
                h_D_onlypt2->Fill(hadron.Getpt2());
                h_wD_pt->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), hadron.Getpt2());    //3 arguments and the WEIGHT
                h_wD_sqpt2->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), hadron.Getpt2()*hadron.Getpt2());    //3 arguments and the WEIGHT (pt2 squared)
                h_D_pt3D->Fill(event.Getxb(), event.GetQ2(), hadron.Getz());    //3 arguments only counts not weight
                hDpt_Q_nu_zD->Fill(event.GetQ2(), event.Getnu(), hadron.Getz()); // uselesss ? 

            }
        }
    }
    else if (targetType == 1 && cuta.PassCutsElectrons(event)==true) {
        counter_elSn++; //counter for only electrons for z
        hDpt_nuA->Fill(event.Getnu());
        for (const Particle& hadron : event.GetHadrons()) {
            if (cuta.PassCutsHadrons(hadron)==true){
                h_A_onlypt2->Fill(hadron.Getpt2());
                h_wA_pt->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), hadron.Getpt2());    //3 arguments and the WEIGHT
                h_wA_sqpt2->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), hadron.Getpt2()*hadron.Getpt2());    //3 arguments and the WEIGHT (pt2 squared)
                h_A_pt3D->Fill(event.Getxb(), event.GetQ2(), hadron.Getz());    //3 arguments only counts not weight
                hDpt_Q_nu_zA->Fill(event.GetQ2(), event.Getnu(), hadron.Getz());    //useless I guess 
            }
        }
    }
}

void deltaptsq::calcDpt(){
    //using here X=nud, Y=z, Z=Q2
    int numBinsX = h_wD_pt->GetNbinsX();    //same bins 4 Deut and A 
    int numBinsY = h_wD_pt->GetNbinsY(); 
    int numBinsZ = h_wD_pt->GetNbinsZ(); 
    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {

        for (int Ybin=1; Ybin <= numBinsY; Ybin++) {
            for (int Zbin= 1; Zbin <= numBinsZ; Zbin++){
                // loop with 125 values (still using 5 per axis-binning))
                double valD = h_wD_pt->GetBinContent(Xbin,Ybin,Zbin);       //weighted histo 3D
                double valA = h_wA_pt->GetBinContent(Xbin,Ybin,Zbin);
                double countD = h_D_pt3D->GetBinContent(Xbin,Ybin,Zbin);    //3D histo same values w/o weight
                double countA = h_A_pt3D->GetBinContent(Xbin,Ybin,Zbin);    //same for A
                double sqvalD = h_wD_sqpt2->GetBinContent(Xbin,Ybin,Zbin);  //square-weighted histo 3D
                double sqvalA = h_wA_sqpt2->GetBinContent(Xbin,Ybin,Zbin);  //same for A
                double wavg_ptD = (countD != 0) ? valD/countD : 0.0;
                double wavg_ptA = (countA != 0) ? valA/countA : 0.0;
                double dpt_point = wavg_ptA - wavg_ptD;
                double varianceD = (countD != 0) ? sqvalD/countD - wavg_ptD*wavg_ptD : 0.0;
                double varianceA = (countA != 0) ? sqvalA/countA - wavg_ptA*wavg_ptA : 0.0;
                double Err_valD = (countD != 0) ? sqrt(varianceD/countD) : 0.0;
                double Err_valA = (countA != 0) ? sqrt(varianceA/countA) : 0.0;
                double Err_dpt_point = sqrt(Err_valD*Err_valD + Err_valA*Err_valA);
                DptMatrix[Xbin-1][Ybin-1][Zbin-1] = dpt_point;
                errorDptMatrix[Xbin-1][Ybin-1][Zbin-1] = Err_dpt_point;
            } 
        }
    }

}
        
void deltaptsq::writeMatrixToFile(const std::string& filename){
    std::ofstream outputFile(filename);
    std::cout<< " Dptbin? -> "<<Dptbin<<std::endl;
    //x=xbj; y= Q2; z=z;
    for (int x = 0; x < Dptbin; ++x) {
        double xaxisval = h_wD_pt->GetXaxis()->GetBinCenter(x + 1);
        outputFile << "xb = " << xaxisval << std::endl;
        for (int y = 0; y < Dptbin; ++y) {
            double yaxisval = h_wD_pt->GetYaxis()->GetBinCenter(y + 1);
            for (int z = 0; z < Dptbin; ++z) {
                double zaxisval = h_wD_pt->GetZaxis()->GetBinCenter(z + 1);
                double value = DptMatrix[x][y][z];
                double error = errorDptMatrix[x][y][z];
                outputFile << "("<< xaxisval<<" ,"<<yaxisval <<" ,"<< zaxisval<< " ) = "<< value << " +- " << error << "\t\t";
                //std::cout << value << " +- " << error << "\t\t";
            }
            outputFile << std::endl;
        }
        outputFile << std::endl << std::endl;
    }
    outputFile.close();
}

void deltaptsq::multiplotDpt(){
    //x= xb, y=Q2, z=z

    for (int x = 0; x < Dptbin; ++x) {
        double xValue = h_wD_pt->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "multiplotDpt_x" + std::to_string(xValue) + ".pdf";
        TCanvas canvasDpt("c", "Multiplot Dpt", 1200, 800);
        canvasDpt.Divide(3, 2); 

        for (int y=0; y < Dptbin; ++y) {
            double Q2Value = h_wD_pt->GetYaxis()->GetBinCenter(y + 1);
            TGraphErrors *graphDpt = new TGraphErrors();
            canvasDpt.cd(y + 1);
            for (int z=0; z < Dptbin; ++z) {
                double zValue = h_wD_pt->GetZaxis()->GetBinCenter(z + 1);
                double value = DptMatrix[x][y][z];
                double error = errorDptMatrix[x][y][z];
                //std::cout << "("<< x<<" ,"<<y <<" ,"<< z<< " ) = "<<std::endl;   // value << " +- " << error << "\t\t";

                graphDpt->SetPoint(z, zValue, value);
                graphDpt->SetPointError(z, 0.0, error); 


                //graph_rat->SetPoint(graph_pointnb, value, error);
                //graph_pointnb++;
            }
            graphDpt->SetTitle(("#Delta <p_{t}^{2}> vs z, Q^{2}=" + std::to_string(Q2Value)).c_str());
            graphDpt->GetXaxis()->SetTitle("z");
            graphDpt->GetYaxis()->SetTitle("#Delta <p_{t}^{2}>");
            graphDpt->SetMarkerStyle(20);
            graphDpt->GetYaxis()->SetRangeUser(-1.0, 1.0); // Set Y axis range from -0.1 to 0.1
            graphDpt->Draw("AP");
            TLine *line = new TLine(graphDpt->GetXaxis()->GetXmin(), 0.0, graphDpt->GetXaxis()->GetXmax(), 0.0);
            line->SetLineStyle(2); // Dotted line
            line->Draw("same");



        }
        canvasDpt.SaveAs(pdfFileName.c_str());
        
    }
    
}