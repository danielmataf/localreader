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
#include "sratio.h"
#include "CutSet.h"
#include <TMultiGraph.h>
#include <TLegend.h>


sratio::sratio(CutSet cutsD, CutSet cutsA, const std::string& targetName): //: cutsD(cutsD), cutsSn(cutsSn) {
    targetName(targetName),
    SratioMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    errorSratioMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),

    hSratio_Q_nu_zD ( new TH3F(("Sratio:nu,z,pt2,D"+targetName).c_str(), ("Sratio:histo nu,z,pt2 for D"+targetName).c_str(),Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ, Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    hSratio_Q_nu_zA ( new TH3F(("Sratio:nu,z,pt2,A"+targetName).c_str(), ("Sratio:histo nu,z,pt2 for A"+targetName).c_str(),Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ, Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    hSratio_nu_z_pt2A_onlye ( new TH3F(("Sratio:Q2, nu,z,A onlye"+targetName).c_str(), ("Sratio:histo_e Q2,nu,z for A"+targetName).c_str(), Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ,Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ),
    hSratio_nu_z_pt2D_onlye ( new TH3F(("Sratio:Q2, nu,z,D onlye"+targetName).c_str(), ("Sratio:histo_e Q2,nu,z for A"+targetName).c_str(), Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ,Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ),
    hSratio_nuA ( new TH1F(("Sratio:nu_A"+targetName).c_str(), ("Sratio:nu_A"+targetName).c_str(), Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio)) ,
    hSratio_nuD ( new TH1F(("Sratio:nu_D"+targetName).c_str(), ("Sratio:nu_D"+targetName).c_str(), Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio)) ,
    h_wD_Sratio ( new TH3F(("wD_Sratio"+targetName).c_str(), ("wD_Sratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    h_wA_Sratio ( new TH3F(("wA_Sratio"+targetName).c_str(), ("wA_Sratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    h_D_Sratio3D ( new TH3F(("countSratio:D"+targetName).c_str(), ("count:wD_Sratio"+targetName).c_str(), Constants::Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    h_A_Sratio3D ( new TH3F(("countSratio:A"+targetName).c_str(), ("count:wD_Sratio"+targetName).c_str(), Constants::Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    h_wD_sqSratio ( new TH3F(("wD_sqSratio"+targetName).c_str(), ("wD_sqSratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ),//definition w/ 3 args
    h_wA_sqSratio ( new TH3F(("wA_sqSratio"+targetName).c_str(), ("wA_sqSratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ){


    cutd = cutsD;
    cuta = cutsA;
    }



    void sratio::FillHistograms(const Event& event){
        int targetType = event.GetTargetType();
        if (targetType == 0 && cutd.PassCutsElectrons(event)==true) {
            counter_elLD2 ++;
            //set a counter that increases when electroncuts = passed; in order for it to be called when R is  computed in had variables (?) TBD
            hSratio_nuD->Fill(event.Getnu());
            for (const Particle& hadron : event.GetHadrons()) {
                if (cutd.PassCutsHadrons(hadron)==true){
                    double phiD = hadron.Getphih();
                    h_wD_Sratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), sin(phiD));    //3 arguments and the WEIGHT
                    h_wD_sqSratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), sin(phiD)*sin(phiD));    //3 arguments and the WEIGHT (pt2 squared) 4 variance
                    h_D_Sratio3D->Fill(event.Getxb(), event.GetQ2(), hadron.Getz());    //3 arguments only counts not weight (cphi)
                }
            }
        }
        else if (targetType == 1 && cuta.PassCutsElectrons(event)==true) {
            counter_elSn++; //counter for only electrons for z and pt
            //here change the else if to just else in order to have a generic target 
            hSratio_nuA->Fill(event.Getnu());
            for (const Particle& hadron : event.GetHadrons()) {
                if (cuta.PassCutsHadrons(hadron)==true){
                    double phiA = hadron.Getphih();
                    h_wA_Sratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), sin(phiA));    //3 arguments and the WEIGHT
                    h_wA_sqSratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), sin(phiA)*sin(phiA));    //3 arguments and the WEIGHT (pt2 squared)
                    h_A_Sratio3D->Fill(event.Getxb(), event.GetQ2(), hadron.Getz());    //3 arguments only counts not weight

                }
            }
        }
    }


void sratio::calcSratio(){
    int numBinsX = h_wD_Sratio->GetNbinsX();    //same bins 4 Deut and A 
    int numBinsY = h_wD_Sratio->GetNbinsY(); 
    int numBinsZ = h_wD_Sratio->GetNbinsZ(); 
    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {
        for (int Ybin=1; Ybin <= numBinsY; Ybin++) {
            for (int Zbin= 1; Zbin <= numBinsZ; Zbin++){
                double valD = h_wD_Sratio->GetBinContent(Xbin,Ybin,Zbin);
                double valA = h_wA_Sratio->GetBinContent(Xbin,Ybin,Zbin);
                double countD = h_D_Sratio3D->GetBinContent(Xbin,Ybin,Zbin);    //3D histo no counts
                double countA = h_A_Sratio3D->GetBinContent(Xbin,Ybin,Zbin);
                double sqvalD = h_wD_sqSratio->GetBinContent(Xbin,Ybin,Zbin);
                double sqvalA = h_wA_sqSratio->GetBinContent(Xbin,Ybin,Zbin);
                double wavg_SratioD = (countD != 0) ? valD/countD : 0.0;   //weighted average ok 
                double wavg_SratioA = (countA != 0) ? valA/countA : 0.0;
                double Sratio_point = (wavg_SratioD != 0) ? wavg_SratioA / wavg_SratioD : 0.0;


                double varianceD = (countD != 0) ? sqvalD/countD - wavg_SratioD*wavg_SratioD : 0.0;
                double varianceA = (countA != 0) ? sqvalA/countA - wavg_SratioA*wavg_SratioA : 0.0;
                double Err_valD = (countD != 0) ? sqrt(varianceD/countD) : 0.0;
                double Err_valA = (countA != 0) ? sqrt(varianceA/countA) : 0.0;
                //double Err_Cratio_point = Cratio_point * sqrt(Err_valD*Err_valD + Err_valA*Err_valA);
                //double Err_Cratio_point = Cratio_point * sqrt(pow(Err_valA / wavg_CratioA, 2) + pow(Err_valD / wavg_CratioD, 2));
                double Err_Sratio_point = (wavg_SratioA != 0 && wavg_SratioD != 0) ? Sratio_point * sqrt(pow(Err_valA / wavg_SratioA, 2) + pow(Err_valD / wavg_SratioD, 2)) : 0.0;


//                double Err_dpt_point = sqrt(Err_valD*Err_valD + Err_valA*Err_valA);
                SratioMatrix[Xbin-1][Ybin-1][Zbin-1] = Sratio_point;
                errorSratioMatrix[Xbin-1][Ybin-1][Zbin-1] = Err_Sratio_point;

            }
        }
    }  



}




void sratio::writeMatrixToFile(const std::string& filename){
    //not needed unless output in file wanted 
    std::ofstream outputFile(filename);
    for (int x = 0; x < Cratiobin; ++x) {
        double xaxisval = h_wD_Sratio->GetXaxis()->GetBinCenter(x + 1);
        outputFile << "xb = " << xaxisval << std::endl;
        for (int y =0; y< Cratiobin; ++y) {
            for (int z = 0; z < Cratiobin; ++z) {
                double value = SratioMatrix[x][y][z];
                double error = errorSratioMatrix[x][y][z];
                outputFile << value << " +- " << error << "\t\t";
            }
            outputFile << std::endl;
        }
        outputFile << std::endl << std::endl;
    }
    outputFile.close();
}
    



void sratio::multiplotSratio(){
    for (int x = 0; x < Cratiobin; ++x){
        double xValue = h_wD_Sratio->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "multiplotSratio_x" + std::to_string(xValue) + ".pdf";
        TCanvas canvasSratio("c", "Multiplot Sratio", 1200, 800);
        canvasSratio.Divide(3, 2); 
        for (int y=0; y < Cratiobin; ++y) {
            double Q2Value = h_wD_Sratio->GetYaxis()->GetBinCenter(y + 1);
            canvasSratio.cd(y + 1);
            TGraphErrors *graphSratio = new TGraphErrors();
            for (int z=0; z < Cratiobin; ++z) {
                double zValue = h_wD_Sratio->GetZaxis()->GetBinCenter(z + 1);
                double value =SratioMatrix[x][y][z];
                double error = errorSratioMatrix[x][y][z];
                graphSratio->SetPoint(z, zValue, value);
                graphSratio->SetPointError(z, 0.0, error); 

                
            }
            graphSratio->SetTitle(("<sin #phi_{h}>_A / <sin #phi_{h}>_D vs z, Q^{2}=" + std::to_string(Q2Value)).c_str());
            graphSratio->GetXaxis()->SetTitle("z");
            graphSratio->GetYaxis()->SetTitle("<sin #phi_{h}>_A / <sin #phi_{h}>_D");
            graphSratio->SetMarkerStyle(20);
            graphSratio->GetYaxis()->SetRangeUser(-10.0, 10.0); // Set Y axis range from 0.0 to 2.0
            graphSratio->Draw("AP");
            TLine *line = new TLine(graphSratio->GetXaxis()->GetXmin(), 1.0, graphSratio->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line
            line->Draw("same");




        }
        canvasSratio.SaveAs(pdfFileName.c_str());
 
    }

}


void sratio::multiplotSratio( sratio& SratioOther, sratio& SratioThird){
    for (int x = 0; x < Cratiobin; ++x){
        double xValue = h_wD_Sratio->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "tripleTargetSratio_x" + std::to_string(xValue) + ".pdf";
        TCanvas canvasSratio("c", "Multiplot Sratio", 1200, 800);
        canvasSratio.Divide(3, 2); 
        for (int y=0; y < Cratiobin; ++y) {
            TMultiGraph *mg = new TMultiGraph();

            double Q2Value = h_wD_Sratio->GetYaxis()->GetBinCenter(y + 1);
            canvasSratio.cd(y + 1);
            TGraphErrors *graphSratio = new TGraphErrors();
            TGraphErrors *graphSratioOther = new TGraphErrors();
            TGraphErrors *graphSratioThird = new TGraphErrors();
            for (int z=0; z < Cratiobin; ++z) {
                double zValue = h_wD_Sratio->GetZaxis()->GetBinCenter(z + 1);
                double value = SratioMatrix[x][y][z];
                double error = errorSratioMatrix[x][y][z];
                double valueOther = SratioOther.SratioMatrix[x][y][z];
                double errorOther = SratioOther.errorSratioMatrix[x][y][z];
                double valueThird = SratioThird.SratioMatrix[x][y][z];
                double errorThird = SratioThird.errorSratioMatrix[x][y][z];


                graphSratio->SetPoint(z, zValue, value);
                graphSratio->SetPointError(z, 0.0, error);
                graphSratioOther->SetPoint(z, zValue+0.01, valueOther);
                graphSratioOther->SetPointError(z+0.01, 0.0, errorOther); 
                graphSratioThird->SetPoint(z, zValue+0.02, valueThird);
                graphSratioThird->SetPointError(z+0.02, 0.0, errorThird);

                
            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << Q2Value;
            std::string formattedQ2Value = ss.str();

            graphSratio->SetTitle(("<sin #phi_{h}>_A / <sin #phi_{h}>_D vs z, Q^{2}=" + std::to_string(Q2Value)).c_str());
            graphSratio->GetXaxis()->SetTitle("z");
            graphSratio->GetYaxis()->SetTitle("<sin #phi_{h}>_A / <sin #phi_{h}>_D");
            graphSratio->SetMarkerStyle(20);
            graphSratio->GetYaxis()->SetRangeUser(-10.0, 10.0); // Set Y axis range from 0.0 to 2.0
            graphSratio->Draw("AP");
            graphSratioOther->SetMarkerStyle(20);
            graphSratioOther->SetMarkerColor(kBlue);
            //graphCratioOther->Draw("P");
            graphSratio->SetMarkerColor(kRed);
            graphSratioThird->SetMarkerStyle(20);
            graphSratioThird->SetMarkerColor(kGreen);
            //graphCratioThird->Draw("P");

            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graphSratio, "Sn", "lp");
            legend->AddEntry(graphSratioOther, "Cu", "lp");
            legend->AddEntry(graphSratioThird, "CxC", "lp");

            TLine *line = new TLine(graphSratio->GetXaxis()->GetXmin(), 1.0, graphSratio->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line

            mg->Add(graphSratio);
            mg->Add(graphSratioOther);
            mg->Add(graphSratioThird);
            mg->SetTitle(("<sin #phi_{h}>_A / <sin #phi_{h}>_D vs z, Q^{2}=" + formattedQ2Value).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("<sin #phi_{h}>_A / <sin #phi_{h}>_D");
            mg->Draw("APE1");
            legend->Draw("same");
            line->Draw("same");




        }
        canvasSratio.SaveAs(pdfFileName.c_str());
 
    }

}