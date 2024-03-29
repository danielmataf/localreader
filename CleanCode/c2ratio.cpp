#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TPad.h>
#include <fstream>
#include <math.h>
#include <TPDF.h>
#include "Event.h" 
#include "c2ratio.h"
#include "CutSet.h"
#include <TMultiGraph.h>
#include <TLegend.h>


c2ratio::c2ratio(CutSet cutsD, CutSet cutsA, const std::string& targetName): //: cutsD(cutsD), cutsSn(cutsSn) {
    targetName(targetName),
    C2ratioMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    errorC2ratioMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),

    hC2ratio_Q_nu_zD ( new TH3F(("C2ratio:nu,z,pt2,D"+targetName).c_str(), ("C2ratio:histo nu,z,pt2 for D"+targetName).c_str(),Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ, Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    hC2ratio_Q_nu_zA ( new TH3F(("C2ratio:nu,z,pt2,A"+targetName).c_str(), ("C2ratio:histo nu,z,pt2 for A"+targetName).c_str(),Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ, Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    hC2ratio_nu_z_pt2A_onlye ( new TH3F(("C2ratio:Q2, nu,z,A onlye"+targetName).c_str(), ("Cratio:histo_e Q2,nu,z for A"+targetName).c_str(), Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ,Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ),
    hC2ratio_nu_z_pt2D_onlye ( new TH3F(("C2ratio:Q2, nu,z,D onlye"+targetName).c_str(), ("Cratio:histo_e Q2,nu,z for A"+targetName).c_str(), Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ,Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ),
    hC2ratio_nuA ( new TH1F(("C2ratio:nu_A"+targetName).c_str(), ("C2ratio:nu_A"+targetName).c_str(), Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio)) ,
    hC2ratio_nuD ( new TH1F(("C2ratio:nu_D"+targetName).c_str(), ("C2ratio:nu_D"+targetName).c_str(), Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio)) ,
    h_wD_C2ratio ( new TH3F(("wD_C2ratio"+targetName).c_str(), ("wD_C2ratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    h_wA_C2ratio ( new TH3F(("wA_C2ratio"+targetName).c_str(), ("wA_C2ratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    h_D_C2ratio3D ( new TH3F(("countC2ratio:D"+targetName).c_str(), ("count:wD_C2ratio"+targetName).c_str(), Constants::Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    h_A_C2ratio3D ( new TH3F(("countC2ratio:A"+targetName).c_str(), ("count:wD_C2ratio"+targetName).c_str(), Constants::Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    h_wD_sqC2ratio ( new TH3F(("wD_sqC2ratio"+targetName).c_str(), ("wD_sqC2ratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ),//definition w/ 3 args
    h_wA_sqC2ratio ( new TH3F(("wA_sqC2ratio"+targetName).c_str(), ("wA_sqC2ratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ){


    cutd = cutsD;
    cuta = cutsA;
    }



    void c2ratio::FillHistograms(const Event& event){
        int targetType = event.GetTargetType();
        if (targetType == 0 && cutd.PassCutsElectrons(event)==true) {
            counter_elLD2 ++;
            //set a counter that increases when electroncuts = passed; in order for it to be called when R is  computed in had variables (?) TBD
            hC2ratio_nuD->Fill(event.Getnu());
            for (const Particle& hadron : event.GetHadrons()) {
                if (cutd.PassCutsHadrons(hadron)==true){
                    double phiD = hadron.Getphih();
                    h_wD_C2ratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), cos(2 * phiD));    //3 arguments and the WEIGHT
                    h_wD_sqC2ratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), cos(2 * phiD)*cos(2 * phiD));    //3 arguments and the WEIGHT (pt2 squared) 4 variance
                    h_D_C2ratio3D->Fill(event.Getxb(), event.GetQ2(), hadron.Getz());    //3 arguments only counts not weight (cphi)
                }
            }
        }
        else if (targetType == 1 && cuta.PassCutsElectrons(event)==true) {
            counter_elSn++; //counter for only electrons for z and pt
            //here change the else if to just else in order to have a generic target 
            hC2ratio_nuA->Fill(event.Getnu());
            for (const Particle& hadron : event.GetHadrons()) {
                if (cuta.PassCutsHadrons(hadron)==true){
                    double phiA = hadron.Getphih();
                    h_wA_C2ratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), cos(2 * phiA));    //3 arguments and the WEIGHT
                    h_wA_sqC2ratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), cos(2 * phiA)*cos(2 * phiA));    //3 arguments and the WEIGHT (pt2 squared)
                    h_A_C2ratio3D->Fill(event.Getxb(), event.GetQ2(), hadron.Getz());    //3 arguments only counts not weight

                }
            }
        }
    }


void c2ratio::calcC2ratio(){
    int numBinsX = h_wD_C2ratio->GetNbinsX();    //same bins 4 Deut and A 
    int numBinsY = h_wD_C2ratio->GetNbinsY(); 
    int numBinsZ = h_wD_C2ratio->GetNbinsZ(); 
    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {
        for (int Ybin=1; Ybin <= numBinsY; Ybin++) {
            for (int Zbin= 1; Zbin <= numBinsZ; Zbin++){
                double valD = h_wD_C2ratio->GetBinContent(Xbin,Ybin,Zbin);
                double valA = h_wA_C2ratio->GetBinContent(Xbin,Ybin,Zbin);
                double countD = h_D_C2ratio3D->GetBinContent(Xbin,Ybin,Zbin);    //3D histo no counts
                double countA = h_A_C2ratio3D->GetBinContent(Xbin,Ybin,Zbin);
                double sqvalD = h_wD_sqC2ratio->GetBinContent(Xbin,Ybin,Zbin);
                double sqvalA = h_wA_sqC2ratio->GetBinContent(Xbin,Ybin,Zbin);
                double wavg_C2ratioD = (countD != 0) ? valD/countD : 0.0;   //weighted average ok 
                double wavg_C2ratioA = (countA != 0) ? valA/countA : 0.0;
                double C2ratio_point = (wavg_C2ratioD != 0) ? wavg_C2ratioA / wavg_C2ratioD : 0.0;


                double varianceD = (countD != 0) ? sqvalD/countD - wavg_C2ratioD*wavg_C2ratioD : 0.0;
                double varianceA = (countA != 0) ? sqvalA/countA - wavg_C2ratioA*wavg_C2ratioA : 0.0;
                double Err_valD = (countD != 0) ? sqrt(varianceD/countD) : 0.0;
                double Err_valA = (countA != 0) ? sqrt(varianceA/countA) : 0.0;
                //double Err_Cratio_point = Cratio_point * sqrt(Err_valD*Err_valD + Err_valA*Err_valA);
                //double Err_Cratio_point = Cratio_point * sqrt(pow(Err_valA / wavg_CratioA, 2) + pow(Err_valD / wavg_CratioD, 2));
                double Err_C2ratio_point = (wavg_C2ratioA != 0 && wavg_C2ratioD != 0) ? C2ratio_point * sqrt(pow(Err_valA / wavg_C2ratioA, 2) + pow(Err_valD / wavg_C2ratioD, 2)) : 0.0;

//                double Err_dpt_point = sqrt(Err_valD*Err_valD + Err_valA*Err_valA);
                C2ratioMatrix[Xbin-1][Ybin-1][Zbin-1] = C2ratio_point;
                errorC2ratioMatrix[Xbin-1][Ybin-1][Zbin-1] = Err_C2ratio_point;

            }
        }
    }  



}




void c2ratio::writeMatrixToFile(const std::string& filename){
    //not needed unless output in file wanted 
    std::ofstream outputFile(filename);
    for (int x = 0; x < Cratiobin; ++x) {
        double xaxisval = h_wD_C2ratio->GetXaxis()->GetBinCenter(x + 1);
        outputFile << "xb = " << xaxisval << std::endl;
        for (int y =0; y< Cratiobin; ++y) {
            for (int z = 0; z < Cratiobin; ++z) {
                double value = C2ratioMatrix[x][y][z];
                double error = errorC2ratioMatrix[x][y][z];
                outputFile << value << " +- " << error << "\t\t";
            }
            outputFile << std::endl;
        }
        outputFile << std::endl << std::endl;
    }
    outputFile.close();
}
    



void c2ratio::multiplotC2ratio(){
    for (int x = 0; x < Cratiobin; ++x){
        double xValue = h_wD_C2ratio->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "multiplotC2ratio_x" + std::to_string(xValue) + ".pdf";
        TCanvas canvasC2ratio("c", "Multiplot C2ratio", 1200, 800);
        canvasC2ratio.Divide(3, 2); 
        for (int y=0; y < Cratiobin; ++y) {
            double Q2Value = h_wD_C2ratio->GetYaxis()->GetBinCenter(y + 1);
            canvasC2ratio.cd(y + 1);
            TGraphErrors *graphC2ratio = new TGraphErrors();
            for (int z=0; z < Cratiobin; ++z) {
                double zValue = h_wD_C2ratio->GetZaxis()->GetBinCenter(z + 1);
                double value = C2ratioMatrix[x][y][z];
                double error = errorC2ratioMatrix[x][y][z];
                graphC2ratio->SetPoint(z, zValue, value);
                graphC2ratio->SetPointError(z, 0.0, error); 

                
            }
            graphC2ratio->SetTitle(("<cos 2 #phi_{h}>_A / <cos 2 #phi_{h}>_D vs z, Q^{2}=" + std::to_string(Q2Value)).c_str());
            graphC2ratio->GetXaxis()->SetTitle("z");
            graphC2ratio->GetYaxis()->SetTitle("<cos 2 #phi_{h}>_A / <cos 2 #phi_{h}>_D");
            graphC2ratio->SetMarkerStyle(20);
            graphC2ratio->GetYaxis()->SetRangeUser(-10.0, 10.0); // Set Y axis range from 0.0 to 2.0
            graphC2ratio->Draw("AP");
            TLine *line = new TLine(graphC2ratio->GetXaxis()->GetXmin(), 1.0, graphC2ratio->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line
            line->Draw("same");




        }
        canvasC2ratio.SaveAs(pdfFileName.c_str());
 
    }

}


void c2ratio::multiplotC2ratio( c2ratio& c2ratioOther, c2ratio& c2ratioThird){
    for (int x = 0; x < Cratiobin; ++x){
        double xValue = h_wD_C2ratio->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "tripleTargetC2ratio_x" + std::to_string(xValue) + ".pdf";
        TCanvas canvasC2ratio("c", "Multiplot C2ratio", 1200, 800);
        canvasC2ratio.Divide(3, 2); 
        for (int y=0; y < Cratiobin; ++y) {
            TMultiGraph *mg = new TMultiGraph();

            double Q2Value = h_wD_C2ratio->GetYaxis()->GetBinCenter(y + 1);
            canvasC2ratio.cd(y + 1);
            TGraphErrors *graphC2ratio = new TGraphErrors();
            TGraphErrors *graphC2ratioOther = new TGraphErrors();
            TGraphErrors *graphC2ratioThird = new TGraphErrors();
            for (int z=0; z < Cratiobin; ++z) {
                double zValue = h_wD_C2ratio->GetZaxis()->GetBinCenter(z + 1);
                double value = C2ratioMatrix[x][y][z];
                double error = errorC2ratioMatrix[x][y][z];
                double valueOther = c2ratioOther.C2ratioMatrix[x][y][z];
                double errorOther = c2ratioOther.errorC2ratioMatrix[x][y][z];
                double valueThird = c2ratioThird.C2ratioMatrix[x][y][z];
                double errorThird = c2ratioThird.errorC2ratioMatrix[x][y][z];


                graphC2ratio->SetPoint(z, zValue, value);
                graphC2ratio->SetPointError(z, 0.0, error);
                graphC2ratioOther->SetPoint(z, zValue+0.01, valueOther);
                graphC2ratioOther->SetPointError(z+0.01, 0.0, errorOther); 
                graphC2ratioThird->SetPoint(z, zValue+0.02, valueThird);
                graphC2ratioThird->SetPointError(z+0.02, 0.0, errorThird);

                
            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << Q2Value;
            std::string formattedQ2Value = ss.str();

            graphC2ratio->SetTitle(("<cos #phi_{h}>_A / <cos #phi_{h}>_D vs z, Q^{2}=" + std::to_string(Q2Value)).c_str());
            graphC2ratio->GetXaxis()->SetTitle("z");
            graphC2ratio->GetYaxis()->SetTitle("<cos #phi_{h}>_A / <cos #phi_{h}>_D");
            graphC2ratio->SetMarkerStyle(20);
            graphC2ratio->GetYaxis()->SetRangeUser(-10.0, 10.0); // Set Y axis range from 0.0 to 2.0
            graphC2ratio->Draw("AP");
            graphC2ratioOther->SetMarkerStyle(20);
            graphC2ratioOther->SetMarkerColor(kBlue);
            //graphCratioOther->Draw("P");
            graphC2ratio->SetMarkerColor(kRed);
            graphC2ratioThird->SetMarkerStyle(20);
            graphC2ratioThird->SetMarkerColor(kGreen);
            //graphCratioThird->Draw("P");

            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graphC2ratio, "Sn", "lp");
            legend->AddEntry(graphC2ratioOther, "Cu", "lp");
            legend->AddEntry(graphC2ratioThird, "CxC", "lp");

            TLine *line = new TLine(graphC2ratio->GetXaxis()->GetXmin(), 1.0, graphC2ratio->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line

            mg->Add(graphC2ratio);
            mg->Add(graphC2ratioOther);
            mg->Add(graphC2ratioThird);
            mg->SetTitle(("<cos 2 #phi_{h}>_A / <cos 2 #phi_{h}>_D vs z, Q^{2}=" + formattedQ2Value).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("<cos 2 #phi_{h}>_A / <cos 2 #phi_{h}>_D");
            mg->Draw("APE1");
            legend->Draw("same");
            line->Draw("same");




        }
        canvasC2ratio.SaveAs(pdfFileName.c_str());
 
    }

}