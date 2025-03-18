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
    CosPhiA_Matrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    CosPhiD_Matrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    CratioMatrixbis(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    errorCratioMatrixbis(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    errorCosPhiA_Matrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    errorCosPhiD_Matrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),

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

//commenting below while adding intermediary steps in this function, this is a savestate
/*
void cratio::calcCratio(){
    //Using here X=xb, Y=Q2, Z=z
    int numBinsX = h_wD_Cratio->GetNbinsX();    //same bins 4 Deut and A 
    int numBinsY = h_wD_Cratio->GetNbinsY(); 
    int numBinsZ = h_wD_Cratio->GetNbinsZ(); 
    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {
        for (int Ybin=1; Ybin <= numBinsY; Ybin++) {
            //as long as X and Y are not a variable from hadron we can get bin content at the end of the loops 
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
*/



void cratio::calcCratio() {
    int numBinsX = h_wD_Cratio->GetNbinsX();
    int numBinsY = h_wD_Cratio->GetNbinsY();
    int numBinsZ = h_wD_Cratio->GetNbinsZ();

    std::cout << "\n=== Calculating <cos(phih_h)> for A and D ===" << std::endl;
    std::cout << "X_bin\tY_bin\tZ_bin\t<cos(phih_h)>_A\t<cos(phih_h)>_D\tRatio" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {
        for (int Ybin = 1; Ybin <= numBinsY; Ybin++) {
            for (int Zbin = 1; Zbin <= numBinsZ; Zbin++) {
                double valD = h_wD_Cratio->GetBinContent(Xbin, Ybin, Zbin);
                double valA = h_wA_Cratio->GetBinContent(Xbin, Ybin, Zbin);
                double countD = h_D_Cratio3D->GetBinContent(Xbin, Ybin, Zbin);
                double countA = h_A_Cratio3D->GetBinContent(Xbin, Ybin, Zbin);

                //intermediary <cos(phh)> for A and D separately
                double wavg_CosPhiD = (countD != 0) ? valD / countD : 0.0;
                double wavg_CosPhiA = (countA != 0) ? valA / countA : 0.0;
                CosPhiD_Matrix[Xbin - 1][Ybin - 1][Zbin - 1] = wavg_CosPhiD;
                CosPhiA_Matrix[Xbin - 1][Ybin - 1][Zbin - 1] = wavg_CosPhiA;
                //calc normal ratio 
                double Cratio_point = (wavg_CosPhiD != 0) ? wavg_CosPhiA / wavg_CosPhiD : 0.0;
                CratioMatrix[Xbin - 1][Ybin - 1][Zbin - 1] = Cratio_point;
                double varianceD = (countD != 0) ? h_wD_sqCratio->GetBinContent(Xbin, Ybin, Zbin) / countD - wavg_CosPhiD * wavg_CosPhiD : 0.0;
                double varianceA = (countA != 0) ? h_wA_sqCratio->GetBinContent(Xbin, Ybin, Zbin) / countA - wavg_CosPhiA * wavg_CosPhiA : 0.0;
                double Err_valD = (countD != 0) ? sqrt(varianceD / countD) : 0.0;
                double Err_valA = (countA != 0) ? sqrt(varianceA / countA) : 0.0;
                errorCosPhiD_Matrix[Xbin - 1][Ybin - 1][Zbin - 1] = Err_valD;
                errorCosPhiA_Matrix[Xbin - 1][Ybin - 1][Zbin - 1] = Err_valA;


                double Err_Cratio_point = (wavg_CosPhiA != 0 && wavg_CosPhiD != 0) ?
                    Cratio_point * sqrt(pow(Err_valA / wavg_CosPhiA, 2) + pow(Err_valD / wavg_CosPhiD, 2)) : 0.0;

                errorCratioMatrix[Xbin - 1][Ybin - 1][Zbin - 1] = Err_Cratio_point;

                // print in terminal
                std::cout << Xbin << "\t" << Ybin << "\t" << Zbin << "\t"
                          << std::fixed << std::setprecision(4) << wavg_CosPhiA << "\t"
                          << std::fixed << std::setprecision(4) << Err_valA << "\t"
                          << std::fixed << std::setprecision(4) << wavg_CosPhiD << "\t"
                          << std::fixed << std::setprecision(4) << Err_valD << "\t"
                          << std::fixed << std::setprecision(4) << Cratio_point << "\t"
                          << std::fixed << std::setprecision(4) << Err_Cratio_point << std::endl;
            }
        }
    }
    std::cout << "=== Calculation Complete ===\n" << std::endl;
}




    TH3F* cratio::getHwA(){
        return h_wA_Cratio;
        //h_wA2_Cratio
    } 

    TH3F* cratio::getHA3D(){
        return h_A_Cratio3D;
        //h_A2_Cratio3D
    } 

    TH3* cratio::getHwA2() {
        return h_wA_sqCratio;
        //h_wA2_sqCratio
    }



void cratio::calcCratioCarbon( cratio& cratioOther){
    // Using here X=xb, Y=Q2, Z=z
    TH3F* h_wA2_Cratio = cratioOther.getHwA(); 
    TH3F* h_A2_Cratio3D = cratioOther.getHA3D();
    TH3* h_wA2_sqCratio = cratioOther.getHwA2();

    // Check if histograms are properly initialized
    if (!h_wA2_Cratio || !h_A2_Cratio3D || !h_wA2_sqCratio) {
        std::cerr << "ERROR: One or more histograms from cratioOther are not initialized!" << std::endl;
        return;
    }

    int numBinsX = h_wA_Cratio->GetNbinsX();
    int numBinsY = h_wA_Cratio->GetNbinsY();
    int numBinsZ = h_wA_Cratio->GetNbinsZ();

    // Ensure the matrix has correct dimensions
    CratioMatrixbis.resize(numBinsX, std::vector<std::vector<double>>(numBinsY, std::vector<double>(numBinsZ, 0.0)));
    errorCratioMatrixbis.resize(numBinsX, std::vector<std::vector<double>>(numBinsY, std::vector<double>(numBinsZ, 0.0)));

    // Debugging Output
    std::cout << "numBinsX: " << numBinsX << ", numBinsY: " << numBinsY << ", numBinsZ: " << numBinsZ << std::endl;
    std::cout << "Matrix size: " << CratioMatrixbis.size() << " x " << CratioMatrixbis[0].size() << " x " << CratioMatrixbis[0][0].size() << std::endl;

    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {
        for (int Ybin=1; Ybin <= numBinsY; Ybin++) {
            for (int Zbin= 1; Zbin <= numBinsZ; Zbin++){
                // Bounds check to prevent out-of-range access
                if (Xbin-1 >= CratioMatrixbis.size() || 
                    Ybin-1 >= CratioMatrixbis[0].size() || 
                    Zbin-1 >= CratioMatrixbis[0][0].size()) {
                    std::cerr << "ERROR: Out-of-bounds access! (" << Xbin-1 << ", " << Ybin-1 << ", " << Zbin-1 << ")" << std::endl;
                    continue;
                }

                double valC1 = h_wA2_Cratio->GetBinContent(Xbin,Ybin,Zbin);
                double valC2 = h_wA_Cratio->GetBinContent(Xbin,Ybin,Zbin);
                double countC1 = h_A2_Cratio3D->GetBinContent(Xbin,Ybin,Zbin);
                double countC2 = h_A_Cratio3D->GetBinContent(Xbin,Ybin,Zbin);
                double sqvalC1 = h_wA2_sqCratio->GetBinContent(Xbin,Ybin,Zbin);
                double sqvalC2 = h_wA_sqCratio->GetBinContent(Xbin,Ybin,Zbin);

                double wavg_CratioC1 = (countC1 != 0) ? valC1/countC1 : 0.0;
                double wavg_CratioC2 = (countC2 != 0) ? valC2/countC2 : 0.0;
                double Cratio_point = (wavg_CratioC1 != 0) ? wavg_CratioC2 / wavg_CratioC1 : 0.0;

                double varianceC1 = (countC1 != 0) ? sqvalC1/countC1 - wavg_CratioC1*wavg_CratioC1 : 0.0;
                double varianceC2 = (countC2 != 0) ? sqvalC2/countC2 - wavg_CratioC2*wavg_CratioC2 : 0.0;
                double Err_valC1 = (countC1 != 0) ? sqrt(varianceC1/countC1) : 0.0;
                double Err_valC2 = (countC2 != 0) ? sqrt(varianceC2/countC2) : 0.0;
                double Err_Cratio_point = (wavg_CratioC2 != 0 && wavg_CratioC1 != 0) 
                    ? Cratio_point * sqrt(pow(Err_valC2 / wavg_CratioC2, 2) + pow(Err_valC1 / wavg_CratioC1, 2)) 
                    : 0.0;

                CratioMatrixbis[Xbin-1][Ybin-1][Zbin-1] = Cratio_point;
                errorCratioMatrixbis[Xbin-1][Ybin-1][Zbin-1] = Err_Cratio_point;
            }
        }
    }    
}



void cratio::multiplotCratioBis() {
    std::cout << "starting Cratio plots for x bin..." << std::endl;
    for (int x = 0; x < Cratiobin; ++x) {
        double xValue = h_wA_Cratio->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "sameTargetCratio_x" + std::to_string(xValue) + ".pdf";

        TCanvas canvas("c", "Multiplot Cratio (C/C)", 1200, 800);
        canvas.Divide(3, 2);

        for (int y = 0; y < Cratiobin; ++y) {
            double Q2Value = h_wA_Cratio->GetYaxis()->GetBinCenter(y + 1);
            canvas.cd(y + 1);

            TGraphErrors *graph = new TGraphErrors();
            for (int z = 0; z < Cratiobin; ++z) {
                double zValue = h_wA_Cratio->GetZaxis()->GetBinCenter(z + 1);
                double value = CratioMatrixbis[x][y][z];
                double error = errorCratioMatrixbis[x][y][z];

                graph->SetPoint(z, zValue, value);
                graph->SetPointError(z, 0.0, error);
            }

            // Graph appearance and labels
            graph->SetTitle(("<cos #phi_{h}>_C1 / <cos #phi_{h}>_C2 vs z, Q^{2}=" + std::to_string(Q2Value)).c_str());
            graph->GetXaxis()->SetTitle("z");
            graph->GetYaxis()->SetTitle("<cos #phi_{h}>_C1 / <cos #phi_{h}>_C2");
            graph->SetMarkerStyle(20);
            graph->SetMarkerColor(kRed);
            graph->GetYaxis()->SetRangeUser(0.0, 2.0);
            graph->Draw("AP");

            // Reference line at 1.0
            TLine *line = new TLine(graph->GetXaxis()->GetXmin(), 1.0, graph->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2);
            line->Draw("same");

            // "Preliminary" label
            TLatex* prelimText = new TLatex();
            prelimText->SetTextSize(0.08);
            prelimText->SetTextAngle(45);
            prelimText->SetTextColorAlpha(kGray + 1, 0.3);
            prelimText->SetNDC();
            prelimText->SetTextAlign(22);
            prelimText->DrawLatex(0.5, 0.5, "preliminary");
        }

        // Save the figure
        canvas.SaveAs(pdfFileName.c_str());
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

void cratio::multiplotCratio(cratio& cratioOther, cratio& cratioThird) {
    for (int x = 0; x < Cratiobin; ++x) {
        double xValue = h_wD_Cratio->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "tripleTargetCratio_x" + std::to_string(xValue) + ".pdf";

        // Open PDF file
        TCanvas canvasCratio("c", "Multiplot Cratio", 1200, 800);
        canvasCratio.Divide(3, 2);
        canvasCratio.Print((pdfFileName + "[").c_str());

        //cos(phi_h)_A / cos(phi_h)_D
        for (int y = 0; y < Cratiobin; ++y) {
            canvasCratio.cd(y + 1);
            TMultiGraph *mg = new TMultiGraph();

            double Q2Value = h_wD_Cratio->GetYaxis()->GetBinCenter(y + 1);
            TGraphErrors *graphCratio = new TGraphErrors();
            TGraphErrors *graphCratioOther = new TGraphErrors();
            TGraphErrors *graphCratioThird = new TGraphErrors();

            for (int z = 0; z < Cratiobin; ++z) {
                double zValue = h_wD_Cratio->GetZaxis()->GetBinCenter(z + 1);
                double value = CratioMatrix[x][y][z];
                double error = errorCratioMatrix[x][y][z];
                double valueOther = cratioOther.CratioMatrix[x][y][z];
                double errorOther = cratioOther.errorCratioMatrix[x][y][z];
                double valueThird = cratioThird.CratioMatrix[x][y][z];
                double errorThird = cratioThird.errorCratioMatrix[x][y][z];

                graphCratio->SetPoint(z, zValue, value);
                graphCratio->SetPointError(z, 0.0, error);
                graphCratioOther->SetPoint(z, zValue + 0.01, valueOther);
                graphCratioOther->SetPointError(z + 0.01, 0.0, errorOther);
                graphCratioThird->SetPoint(z, zValue + 0.02, valueThird);
                graphCratioThird->SetPointError(z + 0.02, 0.0, errorThird);
            }

            graphCratio->SetMarkerStyle(20);
            graphCratio->SetMarkerColor(kRed);
            graphCratio->SetMarkerSize(1.2);
            graphCratioOther->SetMarkerStyle(20);
            graphCratioOther->SetMarkerColor(kBlue);
            graphCratioOther->SetMarkerSize(1.2);
            graphCratioThird->SetMarkerStyle(20);
            graphCratioThird->SetMarkerColor(kGreen);
            graphCratioThird->SetMarkerSize(1.2);

            mg->Add(graphCratio);
            mg->Add(graphCratioOther);
            mg->Add(graphCratioThird);
            mg->SetTitle(("Ratio <cos #phi_{h}>_A / <cos #phi_{h}>_D vs z, Q^{2}=" + std::to_string(Q2Value)).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("<cos #phi_{h}>_A / <cos #phi_{h}>_D");
            mg->GetYaxis()->SetRangeUser(-10.0, 10.0);
            mg->Draw("APE1");

            TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend->AddEntry(graphCratio, "Sn", "lp");
            legend->AddEntry(graphCratioOther, "Cu", "lp");
            legend->AddEntry(graphCratioThird, "CxC", "lp");
            legend->Draw("same");

            TLine *line = new TLine(graphCratio->GetXaxis()->GetXmin(), 1.0, graphCratio->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2);
            line->Draw("same");
        }
        canvasCratio.Print(pdfFileName.c_str());

        // cos(phi_h)_A
        canvasCratio.Clear();
        canvasCratio.Divide(3, 2);
        for (int y = 0; y < Cratiobin; ++y) {
            canvasCratio.cd(y + 1);
            TMultiGraph *mgA = new TMultiGraph();
            TGraphErrors *graphCosPhiA = new TGraphErrors();
            TGraphErrors *graphCosPhiA_Other = new TGraphErrors();
            TGraphErrors *graphCosPhiA_Third = new TGraphErrors();

            for (int z = 0; z < Cratiobin; ++z) {
                double zValue = h_wD_Cratio->GetZaxis()->GetBinCenter(z + 1);
                double valueA = CosPhiA_Matrix[x][y][z];
                double errorA = errorCosPhiA_Matrix[x][y][z];
                double valueA_Other = cratioOther.CosPhiA_Matrix[x][y][z];
                double errorA_Other = cratioOther.errorCosPhiA_Matrix[x][y][z];
                double valueA_Third = cratioThird.CosPhiA_Matrix[x][y][z];
                double errorA_Third = cratioThird.errorCosPhiA_Matrix[x][y][z];
                graphCosPhiA->SetPoint(z, zValue, valueA);
                graphCosPhiA->SetPointError(z, 0.0, errorA);
                graphCosPhiA_Other->SetPoint(z, zValue + 0.01, valueA_Other);
                //graphCosPhiA_Other->SetPointError(z, 0.0, errorA_Other);
                graphCosPhiA_Third->SetPoint(z, zValue + 0.02, valueA_Third);
                graphCosPhiA_Third->SetPointError(z, 0.0, errorA_Third);
            }

            graphCosPhiA->SetMarkerStyle(20);
            graphCosPhiA->SetMarkerColor(kRed);
            graphCosPhiA->SetMarkerSize(1.2);
            graphCosPhiA_Other->SetMarkerStyle(20);
            graphCosPhiA_Other->SetMarkerColor(kBlue);
            graphCosPhiA_Other->SetMarkerSize(1.2);
            graphCosPhiA_Third->SetMarkerStyle(20);
            graphCosPhiA_Third->SetMarkerColor(kGreen);
            graphCosPhiA_Third->SetMarkerSize(1.2);

            mgA->Add(graphCosPhiA);
            mgA->Add(graphCosPhiA_Other);
            mgA->Add(graphCosPhiA_Third);
            mgA->SetTitle("cos(phi_h)_A");
            mgA->GetYaxis()->SetRangeUser(-1, 1);
            mgA->Draw("APE1");

            mgA->Draw("APE1");
        }
        canvasCratio.Print(pdfFileName.c_str());

        // cos(phi_h)_D 
        canvasCratio.Clear();
        canvasCratio.Divide(3, 2);
        for (int y = 0; y < Cratiobin; ++y) {
            canvasCratio.cd(y + 1);
            TMultiGraph *mgD = new TMultiGraph();
            TGraphErrors *graphCosPhiD = new TGraphErrors();
            TGraphErrors *graphCosPhiD_Other = new TGraphErrors();
            TGraphErrors *graphCosPhiD_Third = new TGraphErrors();

            for (int z = 0; z < Cratiobin; ++z) {
                double zValue = h_wD_Cratio->GetZaxis()->GetBinCenter(z + 1);
                double valueD = CosPhiD_Matrix[x][y][z];
                double valueD_Other = cratioOther.CosPhiD_Matrix[x][y][z];
                double valueD_Third = cratioThird.CosPhiD_Matrix[x][y][z];

                graphCosPhiD->SetPoint(z, zValue, valueD);
                graphCosPhiD_Other->SetPoint(z, zValue + 0.01, valueD_Other);
                graphCosPhiD_Third->SetPoint(z, zValue + 0.02, valueD_Third);
            }

            mgD->Draw("APE1");
        }
        canvasCratio.Print(pdfFileName.c_str());

        canvasCratio.Print((pdfFileName + "]").c_str());
    }
}



void cratio::multiCratsimus( cratio& cratioCu,  cratio& cratioC1 , cratio& cratioC2 ){
    for (int x = 0; x < Cratiobin; ++x){
        double xValue = h_wD_Cratio->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "multiCratsimus_x" + std::to_string(xValue) + ".pdf";
        TCanvas canvasCratio("c", "Multiplot Cratio", 1200, 800);
        canvasCratio.Divide(3, 2); 
        for (int y=0; y < Cratiobin; ++y) {
            TMultiGraph *mg = new TMultiGraph();

            double Q2Value = h_wD_Cratio->GetYaxis()->GetBinCenter(y + 1);
            canvasCratio.cd(y + 1);
            TGraphErrors *graphCratioSn = new TGraphErrors();
            TGraphErrors *graphCratioCu = new TGraphErrors();
            TGraphErrors *graphCratioC1 = new TGraphErrors();
            TGraphErrors *graphCratioC2 = new TGraphErrors();
            for (int z=0; z < Cratiobin; ++z) {
                double zValue = h_wD_Cratio->GetZaxis()->GetBinCenter(z + 1);
                double valueSn = CratioMatrix[x][y][z];
                double errorSn = errorCratioMatrix[x][y][z];
                double valueCu = cratioCu.CratioMatrix[x][y][z];
                double errorCu = cratioCu.errorCratioMatrix[x][y][z];
                double valueC1 = cratioC1.CratioMatrix[x][y][z];
                double errorC1 = cratioC1.errorCratioMatrix[x][y][z];
                double valueC2 = cratioC2.CratioMatrix[x][y][z];
                double errorC2 = cratioC2.errorCratioMatrix[x][y][z];
                graphCratioSn->SetPoint(z, zValue, valueSn);
                graphCratioSn->SetPointError(z, 0.0, errorSn);
                graphCratioCu->SetPoint(z, zValue+0.01, valueCu);
                graphCratioCu->SetPointError(z+0.01, 0.0, errorCu); 
                graphCratioC1->SetPoint(z, zValue+0.02, valueC1);
                graphCratioC1->SetPointError(z+0.02, 0.0, errorC1);
                graphCratioC2->SetPoint(z, zValue+0.03, valueC2);
                graphCratioC2->SetPointError(z+0.03, 0.0, errorC2);
            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << Q2Value;
            std::string formattedQ2Value = ss.str();
            std::string title = "<cos #phi_{h}>_A / <cos #phi_{h}>_D vs z, Q^{2}= " + formattedQ2Value + " GeV^{2}";

            graphCratioSn->SetTitle(title.c_str());
            graphCratioSn->GetXaxis()->SetTitle("z");
            graphCratioSn->GetYaxis()->SetRangeUser(-10.0, 10.0); // Set Y axis range from 0.0 to 2.0
            graphCratioSn->SetMarkerStyle(20);
            graphCratioCu->SetMarkerStyle(20);
            graphCratioC1->SetMarkerStyle(20);
            graphCratioC2->SetMarkerStyle(20);

            graphCratioSn->SetMarkerColor(kGreen);
            graphCratioCu->SetMarkerColor(kRed);
            graphCratioC1->SetMarkerColor(kBlue);
            graphCratioC2->SetMarkerColor(kBlack);

            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graphCratioSn, "Sn sim", "lp");
            legend->AddEntry(graphCratioCu, "Cu sim", "lp");
            legend->AddEntry(graphCratioC1, "C1 sim", "lp");
            legend->AddEntry(graphCratioC2, "C2 sim", "lp");

            TLine *line = new TLine(graphCratioSn->GetXaxis()->GetXmin(), 1.0, graphCratioSn->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line

            mg->Add(graphCratioSn);
            mg->Add(graphCratioCu);
            mg->Add(graphCratioC1);
            mg->Add(graphCratioC2);

            mg->SetTitle(("<cos #phi_{h}>_A / <cos #phi_{h}>_D vs z, Q^{2}=" + formattedQ2Value).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("<cos #phi_{h}>_A / <cos #phi_{h}>_D");
            mg->SetTitle((title).c_str());
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


void cratio::multiCrattrue( cratio& cratioCu,  cratio& cratioC1 , cratio& cratioC2 ){
    for (int x = 0; x < Cratiobin; ++x){
        double xValue = h_wD_Cratio->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "multiCratsimus_x" + std::to_string(xValue) + ".pdf";
        TCanvas canvasCratio("c", "Multiplot Cratio", 1200, 800);
        canvasCratio.Divide(3, 2); 
        for (int y=0; y < Cratiobin; ++y) {
            TMultiGraph *mg = new TMultiGraph();

            double Q2Value = h_wD_Cratio->GetYaxis()->GetBinCenter(y + 1);
            canvasCratio.cd(y + 1);
            TGraphErrors *graphCratioSn = new TGraphErrors();
            TGraphErrors *graphCratioCu = new TGraphErrors();
            TGraphErrors *graphCratioC1 = new TGraphErrors();
            TGraphErrors *graphCratioC2 = new TGraphErrors();
            for (int z=0; z < Cratiobin; ++z) {
                double zValue = h_wD_Cratio->GetZaxis()->GetBinCenter(z + 1);
                double valueSn = CratioMatrix[x][y][z];
                double errorSn = errorCratioMatrix[x][y][z];
                double valueCu = cratioCu.CratioMatrix[x][y][z];
                double errorCu = cratioCu.errorCratioMatrix[x][y][z];
                double valueC1 = cratioC1.CratioMatrix[x][y][z];
                double errorC1 = cratioC1.errorCratioMatrix[x][y][z];
                double valueC2 = cratioC2.CratioMatrix[x][y][z];
                double errorC2 = cratioC2.errorCratioMatrix[x][y][z];
                graphCratioSn->SetPoint(z, zValue, valueSn);
                graphCratioSn->SetPointError(z, 0.0, errorSn);
                graphCratioCu->SetPoint(z, zValue+0.01, valueCu);
                graphCratioCu->SetPointError(z+0.01, 0.0, errorCu); 
                graphCratioC1->SetPoint(z, zValue+0.02, valueC1);
                graphCratioC1->SetPointError(z+0.02, 0.0, errorC1);
                graphCratioC2->SetPoint(z, zValue+0.03, valueC2);
                graphCratioC2->SetPointError(z+0.03, 0.0, errorC2);
            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << Q2Value;
            std::string formattedQ2Value = ss.str();
            std::string title = "<cos #phi_{h}>_A / <cos #phi_{h}>_D vs z, Q^{2}= " + formattedQ2Value + " GeV^{2}";

            graphCratioSn->SetTitle(title.c_str());
            graphCratioSn->GetXaxis()->SetTitle("z");
            graphCratioSn->GetYaxis()->SetRangeUser(-10.0, 10.0); // Set Y axis range from 0.0 to 2.0
            graphCratioSn->SetMarkerStyle(20);
            graphCratioCu->SetMarkerStyle(20);
            graphCratioC1->SetMarkerStyle(20);
            graphCratioC2->SetMarkerStyle(20);

            graphCratioSn->SetMarkerColor(kGreen);
            graphCratioCu->SetMarkerColor(kRed);
            graphCratioC1->SetMarkerColor(kBlue);
            graphCratioC2->SetMarkerColor(kBlack);

            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graphCratioSn, "Sn", "lp");
            legend->AddEntry(graphCratioCu, "Cu", "lp");
            legend->AddEntry(graphCratioC1, "C1", "lp");
            legend->AddEntry(graphCratioC2, "C2", "lp");

            TLine *line = new TLine(graphCratioSn->GetXaxis()->GetXmin(), 1.0, graphCratioSn->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line

            mg->Add(graphCratioSn);
            mg->Add(graphCratioCu);
            mg->Add(graphCratioC1);
            mg->Add(graphCratioC2);

            mg->SetTitle(("<cos #phi_{h}>_A / <cos #phi_{h}>_D vs z, Q^{2}=" + formattedQ2Value).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("<cos #phi_{h}>_A / <cos #phi_{h}>_D");
            mg->SetTitle((title).c_str());
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

void cratio::multiCratall( cratio& cratioCu,  cratio& cratioC1 , cratio& cratioC2 , cratio& cratioSnsim , cratio& cratioCusim,  cratio& cratioC1sim , cratio& cratioC2sim  ){
    for (int x = 0; x < Cratiobin; ++x){
        double xValue = h_wD_Cratio->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "multiCratall_x" + std::to_string(xValue) + ".pdf";
        TCanvas canvasCratio("c", "Multiplotall Cratio", 1200, 800);
        canvasCratio.Divide(3, 2); 
        for (int y=0; y < Cratiobin; ++y) {
            TMultiGraph *mg = new TMultiGraph();

            double Q2Value = h_wD_Cratio->GetYaxis()->GetBinCenter(y + 1);
            canvasCratio.cd(y + 1);
            TGraphErrors *graphCratioSn = new TGraphErrors();
            TGraphErrors *graphCratioCu = new TGraphErrors();
            TGraphErrors *graphCratioC1 = new TGraphErrors();
            TGraphErrors *graphCratioC2 = new TGraphErrors();
            TGraphErrors *graphCratioSnsim = new TGraphErrors();
            TGraphErrors *graphCratioCusim = new TGraphErrors();
            TGraphErrors *graphCratioC1sim = new TGraphErrors();
            TGraphErrors *graphCratioC2sim = new TGraphErrors();
            for (int z=0; z < Cratiobin; ++z) {
                double zValue = h_wD_Cratio->GetZaxis()->GetBinCenter(z + 1);
                double valueSn = CratioMatrix[x][y][z];
                double errorSn = errorCratioMatrix[x][y][z];
                double valueCu = cratioCu.CratioMatrix[x][y][z];
                double errorCu = cratioCu.errorCratioMatrix[x][y][z];
                double valueC1 = cratioC1.CratioMatrix[x][y][z];
                double errorC1 = cratioC1.errorCratioMatrix[x][y][z];
                double valueC2 = cratioC2.CratioMatrix[x][y][z];
                double errorC2 = cratioC2.errorCratioMatrix[x][y][z];
                double valueSnsim = cratioSnsim.CratioMatrix[x][y][z];
                double errorSnsim = cratioSnsim.errorCratioMatrix[x][y][z];
                double valueCusim = cratioCusim.CratioMatrix[x][y][z];
                double errorCusim = cratioCusim.errorCratioMatrix[x][y][z];
                double valueC1sim = cratioC1sim.CratioMatrix[x][y][z];
                double errorC1sim = cratioC1sim.errorCratioMatrix[x][y][z];
                double valueC2sim = cratioC2sim.CratioMatrix[x][y][z];
                double errorC2sim = cratioC2sim.errorCratioMatrix[x][y][z];

                graphCratioSn->SetPoint(z, zValue, valueSn);
                graphCratioSn->SetPointError(z, 0.0, errorSn);
                graphCratioCu->SetPoint(z, zValue+0.01, valueCu);
                graphCratioCu->SetPointError(z+0.01, 0.0, errorCu); 
                graphCratioC1->SetPoint(z, zValue+0.02, valueC1);
                graphCratioC1->SetPointError(z+0.02, 0.0, errorC1);
                graphCratioC2->SetPoint(z, zValue+0.03, valueC2);
                graphCratioC2->SetPointError(z+0.03, 0.0, errorC2);
                graphCratioSnsim->SetPoint(z, zValue+0.005, valueSnsim);
                graphCratioSnsim->SetPointError(z+0.005, 0.0, errorSnsim);
                graphCratioCusim->SetPoint(z, zValue+0.015, valueCusim);
                graphCratioCusim->SetPointError(z+0.015, 0.0, errorCusim);
                graphCratioC1sim->SetPoint(z, zValue+0.025, valueC1sim);
                graphCratioC1sim->SetPointError(z+0.025, 0.0, errorC1sim);
                graphCratioC2sim->SetPoint(z, zValue+0.035, valueC2sim);
                graphCratioC2sim->SetPointError(z+0.035, 0.0, errorC2sim);

            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << Q2Value;
            std::string formattedQ2Value = ss.str();
            std::string title = "<cos #phi_{h}>_A / <cos #phi_{h}>_D vs z, Q^{2}= " + formattedQ2Value + " GeV^{2}";

            graphCratioSn->SetTitle(title.c_str());
            graphCratioSn->GetXaxis()->SetTitle("z");
            graphCratioSn->GetYaxis()->SetRangeUser(-10.0, 10.0); // Set Y axis range from 0.0 to 2.0
            graphCratioSn->SetMarkerStyle(20);
            graphCratioCu->SetMarkerStyle(20);
            graphCratioC1->SetMarkerStyle(20);
            graphCratioC2->SetMarkerStyle(20);
            graphCratioSnsim->SetMarkerStyle(20);
            graphCratioCusim->SetMarkerStyle(20);
            graphCratioC1sim->SetMarkerStyle(20);
            graphCratioC2sim->SetMarkerStyle(20);


            graphCratioSn->SetMarkerColor(kGreen);
            graphCratioCu->SetMarkerColor(kRed);
            graphCratioC1->SetMarkerColor(kBlue);
            graphCratioC2->SetMarkerColor(kBlack);
            graphCratioSnsim->SetMarkerColor(kSpring-7);
            graphCratioCusim->SetMarkerColor(kOrange+8);
            graphCratioC1sim->SetMarkerColor(kCyan+3);
            graphCratioC2sim->SetMarkerColor(kGray);



            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graphCratioSn, "Sn", "lp");
            legend->AddEntry(graphCratioCu, "Cu", "lp");
            legend->AddEntry(graphCratioC1, "C1", "lp");
            legend->AddEntry(graphCratioC2, "C2", "lp");
            legend->AddEntry(graphCratioSnsim, "Sn sim", "lp");
            legend->AddEntry(graphCratioCusim, "Cu sim", "lp");
            legend->AddEntry(graphCratioC1sim, "C1 sim", "lp");
            legend->AddEntry(graphCratioC2sim, "C2 sim", "lp");


            TLine *line = new TLine(graphCratioSn->GetXaxis()->GetXmin(), 1.0, graphCratioSn->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line

            mg->Add(graphCratioSn);
            mg->Add(graphCratioCu);
            mg->Add(graphCratioC1);
            mg->Add(graphCratioC2);
            mg->Add(graphCratioSnsim);
            mg->Add(graphCratioCusim);
            mg->Add(graphCratioC1sim);
            mg->Add(graphCratioC2sim);


            mg->SetTitle(("<cos #phi_{h}>_A / <cos #phi_{h}>_D vs z, Q^{2}=" + formattedQ2Value).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("<cos #phi_{h}>_A / <cos #phi_{h}>_D");
            mg->SetTitle((title).c_str());
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
void cratio::ValidateHistograms() {
    std::cout << "Validating histograms for Cratio calculations..." << std::endl;
    std::vector<std::pair<TH3F*, std::string>> histograms = {
        {h_wD_Cratio, "h_wD_Cratio"},
        {h_wA_Cratio, "h_wA_Cratio"},
        {h_D_Cratio3D, "h_D_Cratio3D"},
        {h_A_Cratio3D, "h_A_Cratio3D"}
    };

    for (const auto& [hist, name] : histograms) {
        if (!hist || hist->GetEntries() == 0) {
            std::cerr << "Warning: Histogram " << name << " is empty or not initialized!" << std::endl;
        }
    }
}
void cratio::LogBinContent() {
    auto logHistogram = [](TH3F* hist, const std::string& name) {
        int numBinsX = hist->GetNbinsX();
        int numBinsY = hist->GetNbinsY();
        int numBinsZ = hist->GetNbinsZ();

        std::cout << "Contents of histogram " << name << ":" << std::endl;

        for (int x = 1; x <= numBinsX; ++x) {
            double xCenter = hist->GetXaxis()->GetBinCenter(x);
            std::cout << "x = " << xCenter << ":" << std::endl;

            for (int y = 1; y <= numBinsY; ++y) {
                double yCenter = hist->GetYaxis()->GetBinCenter(y);
                std::cout << "  y = " << yCenter << ": ";

                for (int z = 1; z <= numBinsZ; ++z) {
                    double content = hist->GetBinContent(x, y, z);
                    std::cout << content << " ";
                }
                std::cout << std::endl;
            }
        }
    };

    logHistogram(h_wD_Cratio, "h_wD_Cratio");
    logHistogram(h_wA_Cratio, "h_wA_Cratio");
}

cratio::~cratio() {
    delete hCratio_Q_nu_zD;
    delete hCratio_Q_nu_zA;
    delete hCratio_nu_z_pt2A_onlye;
    delete hCratio_nu_z_pt2D_onlye;
    delete hCratio_nuA;
    delete hCratio_nuD;
    delete h_wD_Cratio;
    delete h_wA_Cratio;
    delete h_D_Cratio3D;
    delete h_A_Cratio3D;
    delete h_wD_sqCratio;
    delete h_wA_sqCratio;
    delete graph_crat;
}
