#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLatex.h>  

#include <TPad.h>
#include <fstream>
#include <math.h>
#include <TPDF.h>
#include "Event.h" 
#include "cratio.h"
#include "CutSet.h"
#include <TMultiGraph.h>
#include <TLegend.h>


cratio::cratio(CutSet cutsD, CutSet cutsA, const std::string& targetName): //: cutsD(cutsD), cutsSn(cutsSn) {
    targetName(targetName),
    CratioMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    errorCratioMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),

    hCratio_Q_nu_zD ( new TH3F(("Cratio:nu,z,pt2,D"+targetName).c_str(), ("Cratio:histo nu,z,pt2 for D"+targetName).c_str(),Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ, Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    hCratio_Q_nu_zA ( new TH3F(("Cratio:nu,z,pt2,A"+targetName).c_str(), ("Cratio:histo nu,z,pt2 for A"+targetName).c_str(),Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ, Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    hCratio_nu_z_pt2A_onlye ( new TH3F(("Cratio:Q2, nu,z,A onlye"+targetName).c_str(), ("Cratio:histo_e Q2,nu,z for A"+targetName).c_str(), Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ,Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ),
    hCratio_nu_z_pt2D_onlye ( new TH3F(("Cratio:Q2, nu,z,D onlye"+targetName).c_str(), ("Cratio:histo_e Q2,nu,z for A"+targetName).c_str(), Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ,Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ),
    hCratio_nuA ( new TH1F(("Cratio:nu_A"+targetName).c_str(), ("Cratio:nu_A"+targetName).c_str(), Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio)) ,
    hCratio_nuD ( new TH1F(("Cratio:nu_D"+targetName).c_str(), ("Cratio:nu_D"+targetName).c_str(), Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio)) ,
    h_wD_Cratio ( new TH3F(("wD_Cratio"+targetName).c_str(), ("wD_Cratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    h_wA_Cratio ( new TH3F(("wA_Cratio"+targetName).c_str(), ("wA_Cratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    h_D_Cratio3D ( new TH3F(("countCratio:D"+targetName).c_str(), ("count:wD_Cratio"+targetName).c_str(), Constants::Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    h_A_Cratio3D ( new TH3F(("countCratio:A"+targetName).c_str(), ("count:wD_Cratio"+targetName).c_str(), Constants::Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    h_wD_sqCratio ( new TH3F(("wD_sqCratio"+targetName).c_str(), ("wD_sqCratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ),//definition w/ 3 args
    h_wA_sqCratio ( new TH3F(("wA_sqCratio"+targetName).c_str(), ("wA_sqCratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ){


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
                double wavg_CratioD = (countD != 0) ? valD/countD : 0.0;   //weighted average ok 
                double wavg_CratioA = (countA != 0) ? valA/countA : 0.0;
                double Cratio_point = (wavg_CratioD != 0) ? wavg_CratioA / wavg_CratioD : 0.0;


                double varianceD = (countD != 0) ? sqvalD/countD - wavg_CratioD*wavg_CratioD : 0.0;
                double varianceA = (countA != 0) ? sqvalA/countA - wavg_CratioA*wavg_CratioA : 0.0;
                double Err_valD = (countD != 0) ? sqrt(varianceD/countD) : 0.0;
                double Err_valA = (countA != 0) ? sqrt(varianceA/countA) : 0.0;
                //double Err_Cratio_point = Cratio_point * sqrt(Err_valD*Err_valD + Err_valA*Err_valA);
                //double Err_Cratio_point = Cratio_point * sqrt(pow(Err_valA / wavg_CratioA, 2) + pow(Err_valD / wavg_CratioD, 2));
                double Err_Cratio_point = (wavg_CratioA != 0 && wavg_CratioD != 0) ? Cratio_point * sqrt(pow(Err_valA / wavg_CratioA, 2) + pow(Err_valD / wavg_CratioD, 2)) : 0.0;

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
                graphCratio->SetPoint(z, zValue, value);
                graphCratio->SetPointError(z, 0.0, error); 

                
            }
            graphCratio->SetTitle(("<cos #phi_{h}>_A / <cos #phi_{h}>_D vs z, Q^{2}=" + std::to_string(Q2Value)).c_str());
            graphCratio->GetXaxis()->SetTitle("z");
            graphCratio->GetYaxis()->SetTitle("<cos #phi_{h}>_A / <cos #phi_{h}>_D");
            graphCratio->SetMarkerStyle(20);
            graphCratio->GetYaxis()->SetRangeUser(-10.0, 10.0); // Set Y axis range from 0.0 to 2.0
            graphCratio->Draw("AP");
            TLine *line = new TLine(graphCratio->GetXaxis()->GetXmin(), 1.0, graphCratio->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line
            line->Draw("same");




        }
        canvasCratio.SaveAs(pdfFileName.c_str());
 
    }

}


void cratio::multiplotCratio( cratio& cratioOther, cratio& cratioThird){
    for (int x = 0; x < Cratiobin; ++x){
        double xValue = h_wD_Cratio->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "tripleTargetCratio_x" + std::to_string(xValue) + ".pdf";
        TCanvas canvasCratio("c", "Multiplot Cratio", 1200, 800);
        canvasCratio.Divide(3, 2); 
        for (int y=0; y < Cratiobin; ++y) {
            TMultiGraph *mg = new TMultiGraph();

            double Q2Value = h_wD_Cratio->GetYaxis()->GetBinCenter(y + 1);
            canvasCratio.cd(y + 1);
            TGraphErrors *graphCratio = new TGraphErrors();
            TGraphErrors *graphCratioOther = new TGraphErrors();
            TGraphErrors *graphCratioThird = new TGraphErrors();
            for (int z=0; z < Cratiobin; ++z) {
                double zValue = h_wD_Cratio->GetZaxis()->GetBinCenter(z + 1);
                double value = CratioMatrix[x][y][z];
                double error = errorCratioMatrix[x][y][z];
                double valueOther = cratioOther.CratioMatrix[x][y][z];
                double errorOther = cratioOther.errorCratioMatrix[x][y][z];
                double valueThird = cratioThird.CratioMatrix[x][y][z];
                double errorThird = cratioThird.errorCratioMatrix[x][y][z];


                graphCratio->SetPoint(z, zValue, value);
                graphCratio->SetPointError(z, 0.0, error);
                graphCratioOther->SetPoint(z, zValue+0.01, valueOther);
                graphCratioOther->SetPointError(z+0.01, 0.0, errorOther); 
                graphCratioThird->SetPoint(z, zValue+0.02, valueThird);
                graphCratioThird->SetPointError(z+0.02, 0.0, errorThird);

                
            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << Q2Value;
            std::string formattedQ2Value = ss.str();

            graphCratio->SetTitle(("<cos #phi_{h}>_A / <cos #phi_{h}>_D vs z, Q^{2}=" + std::to_string(Q2Value)).c_str());
            graphCratio->GetXaxis()->SetTitle("z");
            graphCratio->GetYaxis()->SetTitle("<cos #phi_{h}>_A / <cos #phi_{h}>_D");
            graphCratio->SetMarkerStyle(20);
            graphCratio->GetYaxis()->SetRangeUser(-10.0, 10.0); // Set Y axis range from 0.0 to 2.0
            graphCratio->Draw("AP");
            graphCratioOther->SetMarkerStyle(20);
            graphCratioOther->SetMarkerColor(kBlue);
            //graphCratioOther->Draw("P");
            graphCratio->SetMarkerColor(kRed);
            graphCratioThird->SetMarkerStyle(20);
            graphCratioThird->SetMarkerColor(kGreen);
            //graphCratioThird->Draw("P");

            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graphCratio, "Sn", "lp");
            legend->AddEntry(graphCratioOther, "Cu", "lp");
            legend->AddEntry(graphCratioThird, "CxC", "lp");

            TLine *line = new TLine(graphCratio->GetXaxis()->GetXmin(), 1.0, graphCratio->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line

            mg->Add(graphCratio);
            mg->Add(graphCratioOther);
            mg->Add(graphCratioThird);
            mg->SetTitle(("<cos #phi_{h}>_A / <cos #phi_{h}>_D vs z, Q^{2}=" + formattedQ2Value).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("<cos #phi_{h}>_A / <cos #phi_{h}>_D");
            mg->Draw("APE1");
            legend->Draw("same");
            line->Draw("same");
            TLatex* prelimText = new TLatex();
            prelimText->SetTextSize(0.08);  // Larger text size
            prelimText->SetTextAngle(45);
            prelimText->SetTextColorAlpha(kGray + 1, 0.3);  // Gray color with transparency
            prelimText->SetNDC();
            prelimText->SetTextAlign(22);  // Centered alignment
            prelimText->DrawLatex(0.5, 0.5, "preliminary");




        }
        canvasCratio.SaveAs(pdfFileName.c_str());
 
    }

}