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
#include "sratio.h"
#include "CutSet.h"
#include <TMultiGraph.h>
#include <TLegend.h>


sratio::sratio(CutSet cutsD, CutSet cutsA, const std::string& targetName): //: cutsD(cutsD), cutsSn(cutsSn) {
    targetName(targetName),
    SratioMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    errorSratioMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    AsymmMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    errorAsymmMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),

    hSratio_Q_nu_zD ( new TH3F(("Sratio:nu,z,pt2,D"+targetName).c_str(), ("Sratio:histo nu,z,pt2 for D"+targetName).c_str(),Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ, Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    hSratio_Q_nu_zA ( new TH3F(("Sratio:nu,z,pt2,A"+targetName).c_str(), ("Sratio:histo nu,z,pt2 for A"+targetName).c_str(),Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ, Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    hSratio_nu_z_pt2A_onlye ( new TH3F(("Sratio:Q2, nu,z,A onlye"+targetName).c_str(), ("Sratio:histo_e Q2,nu,z for A"+targetName).c_str(), Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ,Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ),
    hSratio_nu_z_pt2D_onlye ( new TH3F(("Sratio:Q2, nu,z,D onlye"+targetName).c_str(), ("Sratio:histo_e Q2,nu,z for A"+targetName).c_str(), Constants::Cratiobin_Q, Constants::RcutminQ, Constants::RcutmaxQ,Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ),
    hSratio_nuA ( new TH1F(("Sratio:nu_A"+targetName).c_str(), ("Sratio:nu_A"+targetName).c_str(), Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio)) ,
    hSratio_nuD ( new TH1F(("Sratio:nu_D"+targetName).c_str(), ("Sratio:nu_D"+targetName).c_str(), Constants::Cratiobin_nu,Constants::numinCratio,numaxCratio)) ,
    hSratio_zA ( new TH1F(("Sratio:z_A"+targetName).c_str(), ("Sratio:z_A"+targetName).c_str(), Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    hSratio_zD ( new TH1F(("Sratio:z_D"+targetName).c_str(), ("Sratio:z_D"+targetName).c_str(), Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    hSratio_zA_w ( new TH1F(("Sratio:z_A_w"+targetName).c_str(), ("Sratio:z_A_w"+targetName).c_str(), Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    hSratio_zD_w ( new TH1F(("Sratio:z_D_w"+targetName).c_str(), ("Sratio:z_D_w"+targetName).c_str(), Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    h_wD_Sratio ( new TH3F(("wD_Sratio"+targetName).c_str(), ("wD_Sratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    h_wA_Sratio ( new TH3F(("wA_Sratio"+targetName).c_str(), ("wA_Sratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    h_D_Sratio3D ( new TH3F(("countSratio:D"+targetName).c_str(), ("count:wD_Sratio"+targetName).c_str(), Constants::Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    h_A_Sratio3D ( new TH3F(("countSratio:A"+targetName).c_str(), ("count:wD_Sratio"+targetName).c_str(), Constants::Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ  )),
    h_wD_sqSratio ( new TH3F(("wD_sqSratio"+targetName).c_str(), ("wD_sqSratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ),//definition w/ 3 args
    h_wA_sqSratio ( new TH3F(("wA_sqSratio"+targetName).c_str(), ("wA_sqSratio"+targetName).c_str(),Constants::Cratiobin_x  , xminCratio, xmaxCratio,Constants::Cratiobin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Cratiobin_z,Constants::RcutminZ, Constants::RcutmaxZ ) ),
    h_phiMonA ( new TH1F(("phiMonA"+targetName).c_str(), ("phiMonA"+targetName).c_str(), 100, 0, 360)),
    h_phiMonD ( new TH1F(("phiMonD"+targetName).c_str(), ("phiMonD"+targetName).c_str(), 100, 0, 360)),

    h_Dminus (new TH3F (( " count3Dm"+targetName).c_str(), ( "count3Dm"+targetName).c_str(),5, Constants::RcutminZ, Constants::RcutmaxZ, 5 , Constants::RcutminPt2, Constants::RcutmaxPt2, 5, 0,360)),
    h_Dplus (new TH3F (( " count3Dp"+targetName).c_str(), ( "count3Dp"+targetName).c_str(), 5, Constants::RcutminZ, Constants::RcutmaxZ, 5 , Constants::RcutminPt2, Constants::RcutmaxPt2, 5, 0,360)),
    h_Aminus (new TH3F (( " count3Am"+targetName).c_str(), ( "count3Am"+targetName).c_str(),5, Constants::RcutminZ, Constants::RcutmaxZ, 5 , Constants::RcutminPt2, Constants::RcutmaxPt2, 5, 0,360)),
    h_Aplus (new TH3F (( " count3Ap"+targetName).c_str(), ( "count3Ap"+targetName).c_str(), 5, Constants::RcutminZ, Constants::RcutmaxZ, 5 , Constants::RcutminPt2, Constants::RcutmaxPt2, 5, 0,360)),

    h_xQA ( new TH2F(("xQA"+targetName).c_str(), ("xQA"+targetName).c_str(), 100, xminCratio, xmaxCratio, 100, Constants::RcutminQ, Constants::RcutmaxQ)),
    h_xQD ( new TH2F(("xQD"+targetName).c_str(), ("xQD"+targetName).c_str(), 100, xminCratio, xmaxCratio, 100, Constants::RcutminQ, Constants::RcutmaxQ)),
    h_xphiA ( new TH2F(("xphiA"+targetName).c_str(), ("xphiA"+targetName).c_str(),  100, 0, 360, 100, xminCratio, xmaxCratio)),
    h_xphiD ( new TH2F(("xphiD"+targetName).c_str(), ("xphiD"+targetName).c_str(),  100, 0, 360, 100, xminCratio, xmaxCratio)),
    h_QphiA ( new TH2F(("QphiA"+targetName).c_str(), ("QphiA"+targetName).c_str(),  100, 0, 360, 100, Constants::RcutminQ, Constants::RcutmaxQ)),
    h_QphiD ( new TH2F(("QphiD"+targetName).c_str(), ("QphiD"+targetName).c_str(),  100, 0, 360, 100, Constants::RcutminQ, Constants::RcutmaxQ)),
    h_zphiA ( new TH2F(("zphiA"+targetName).c_str(), ("zphiA"+targetName).c_str(),  100, 0, 360, 100, Constants::RcutminZ, Constants::RcutmaxZ)),
    h_zphiD ( new TH2F(("zphiD"+targetName).c_str(), ("zphiD"+targetName).c_str(),  100, 0, 360, 100, Constants::RcutminZ, Constants::RcutmaxZ)){



    cutd = cutsD;
    cuta = cutsA;
    }



    void sratio::FillHistograms(const Event& event){
        int targetType = event.GetTargetType();
        int helicity = event.GetHel();

        //retrieve hel to consider "loss" when value is -1. aka gain when value = +1
        //for 'loss' and 'gain' we'll add additionnal weight in count histograms. And adding a helicity factor to the histos already weighted. 
            //std::cout << " S " << std::endl;

        if (targetType == 0 && cutd.PassCutsElectrons(event)==true && cutd.PassCutsDetectors(event)==true) {
            //std::cout << " D " << std::endl;
            counter_elLD2 ++;
            //set a counter that increases when electroncuts = passed; in order for it to be called when R is  computed in had variables (?) TBD
            hSratio_nuD->Fill(event.electron.Getnu(), helicity); //only counts. Weighting with helicity 
            // I think this is useless. We dont care about electron count anymore
            //should be X ????


            for (const Particle& hadron : event.GetHadrons()) {

                if (cutd.PassCutsHadrons(hadron)==true){
                    if (hadron.GetPID() == Constants::PION_PLUS_PID  ){   //adding this condition for pion+ and erasing the condit at evtprocessr

                        double phiD = hadron.Getphih();
                        h_wD_Sratio->Fill(event.electron.Getxb(), event.electron.GetQ2(), hadron.Getz(), sin(phiD)*helicity) ;    //3 arguments and the WEIGHT
                        h_wD_sqSratio->Fill(event.electron.Getxb(), event.electron.GetQ2(), hadron.Getz(), sin(phiD)*sin(phiD)*helicity);    //3 arguments and the WEIGHT (pt2 squared) 4 variance
                        h_D_Sratio3D->Fill(event.electron.Getxb(), event.electron.GetQ2(), hadron.Getz());    //3 arguments only counts not weight (cphi)
                        //std::cout << " D -> x value" <<event.electron.Getxb()<< std::endl;
                        //std::cout << " D -> Q2 value" << event.electron.GetQ2()<< std::endl;
                        //std::cout << " D -> z value" << hadron.Getz() << std::endl;
                        h_xphiD->Fill(phiD, event.electron.Getxb()  );
                        h_QphiD->Fill(phiD, event.electron.GetQ2()  );
                        h_zphiD->Fill(phiD, hadron.Getz() );        
                        h_xQD->Fill(event.electron.Getxb(), event.electron.GetQ2()); 
                        h_phiMonD->Fill(phiD);
                        hSratio_zD->Fill(hadron.Getz());    //*helicity
                        hSratio_zD_w->Fill(hadron.Getz(), sin(phiD)*helicity);
                        if (helicity == -1){
                            h_Dminus->Fill(hadron.Getz(), hadron.Getpt2(), phiD);

                        }
                        else if (helicity == 1){

                            h_Dplus->Fill(hadron.Getz(), hadron.Getpt2(), phiD);
                        }       
                    }
                }
            }
        }
        else if (targetType == 1 && cuta.PassCutsElectrons(event)==true) {
            counter_elSn++; //counter for only electrons for z and pt
            //here change the else if to just else in order to have a generic target 
            hSratio_nuA->Fill(event.Getnu());
            for (const Particle& hadron : event.GetHadrons()) {
                if (cuta.PassCutsHadrons(hadron)==true){
                    if (hadron.GetPID() == Constants::PION_PLUS_PID  ){   //adding this condition for pion+ and erasing the condit at evtprocessr

                        double phiA = hadron.Getphih();
                        h_wA_Sratio->Fill(event.electron.Getxb(), event.electron.GetQ2(), hadron.Getz(), sin(phiA)*helicity);    //3 arguments and the WEIGHT
                        h_wA_sqSratio->Fill(event.electron.Getxb(), event.electron.GetQ2(), hadron.Getz(), sin(phiA)*sin(phiA)*helicity);    //3 arguments and the WEIGHT (pt2 squared)
                        h_A_Sratio3D->Fill(event.electron.Getxb(), event.electron.GetQ2(), hadron.Getz());    //3 arguments only counts not weight

                        h_xphiA->Fill(phiA, event.Getxb()  );
                        h_QphiA->Fill(phiA, event.GetQ2()  );
                        h_zphiA->Fill(phiA, hadron.Getz()  );
                        h_xQA->Fill(event.Getxb(), event.GetQ2());
                        h_phiMonA->Fill(phiA);
                        hSratio_zA->Fill(hadron.Getz() );   //*helicity
                        hSratio_zA_w->Fill(hadron.Getz(), sin(phiA)*helicity);
                        if (helicity == -1){
                            h_Aminus->Fill(hadron.Getz(), hadron.Getpt2(), phiA);
                        }
                        else if (helicity == 1){
                            h_Aplus->Fill(hadron.Getz(), hadron.Getpt2(), phiA);
                        }       

                    }
                }
            }
        }
    }

void sratio::calcSratio1D(){
    int numBins1D = hSratio_zD->GetNbinsX();    //supp: same bins 4 Deut and A
    for (int bin1D=1 ; bin1D <= numBins1D; bin1D++){
        double valD = hSratio_zD_w->GetBinContent(bin1D);
        double valA = hSratio_zA_w->GetBinContent(bin1D);
        double countD = hSratio_zA->GetBinContent(bin1D);    //3D histo no counts
        double countA = hSratio_zA->GetBinContent(bin1D);
        double wavg_SratioD = (countD != 0) ? valD/countD : 0.002;   //weighted average ok 
        double wavg_SratioA = (countA != 0) ? valA/countA : 0.001;
        double Sratio_point = (wavg_SratioD != 0) ? wavg_SratioA / wavg_SratioD : 0.003;
        //std::cout << "wavA = valA/countA = "<< wavg_SratioA << " = " << valA << " / "<< countA << std::endl; 
        //std::cout << "wavA = valD/countD = "<< wavg_SratioD << " = " << valD << " / "<< countD << std::endl; 
        //std::cout << "Sratio_point = wavA/wavD =  " << Sratio_point << " = "<<  wavg_SratioA << " /  " <<wavg_SratioD  <<std::endl;
        
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
                double wavg_SratioD = (countD != 0) ? valD/countD : 0.01;   //weighted average ok 
                double wavg_SratioA = (countA != 0) ? valA/countA : 0.001;
                double Sratio_point = (wavg_SratioD != 0) ? wavg_SratioA / wavg_SratioD : 0.01;
                //double Sratio_point =  wavg_SratioD ;
                //std::cout << "Sratio_point = " << Sratio_point << std::endl;

                double varianceD = (countD != 0) ? sqvalD/countD - wavg_SratioD*wavg_SratioD : 0.0;
                double varianceA = (countA != 0) ? sqvalA/countA - wavg_SratioA*wavg_SratioA : 0.0;
                double Err_valD =  (countD != 0) ? sqrt(varianceD/countD) : 0.0;
                double Err_valA =  (countA != 0) ? sqrt(varianceA/countA) : 0.0;
                //double Err_Cratio_point = Cratio_point * sqrt(Err_valD*Err_valD + Err_valA*Err_valA);
                //double Err_Cratio_point = Cratio_point * sqrt(pow(Err_valA / wavg_CratioA, 2) + pow(Err_valD / wavg_CratioD, 2));
                double Err_Sratio_point = (wavg_SratioA != 0 && wavg_SratioD != 0) ? Sratio_point * sqrt(pow(Err_valA / wavg_SratioA, 2) + pow(Err_valD / wavg_SratioD, 2)) : 0.5;


                double Err_dpt_point = sqrt(Err_valD*Err_valD + Err_valA*Err_valA);
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
            graphSratio->GetYaxis()->SetRangeUser(-5.0, 5.0); // Set Y axis range from 0.0 to 2.0
            graphSratio->Draw("AP");
            TLine *line = new TLine(graphSratio->GetXaxis()->GetXmin(), 1.0, graphSratio->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line
            line->Draw("same");




        }
        canvasSratio.SaveAs(pdfFileName.c_str());
 
    }

}

void sratio::calcAsymmetries(){
    //TCanvas canvasAsy("c", "asy Sratiomon", 1200, 800);
    //canvasAsy.Divide(1, 2);
    //TH3 *h_num = (TH3*) h_Dplus->Clone("h_num");
    //TH3 *h_den = (TH3*) h_Dplus->Clone("h_den");
    //TH3 *h_A;
    //h_num->Add(h_Dminus, -1);
    //h_den ->Add(h_Dminus);
    //h_A = (TH3*) h_num->Clone("h_as"); 
    //h_A->Divide(h_den);
//
    //delete h_num;
    //delete h_den;
    //canvasAsy.cd(1);
    //h_A->SetTitle("AsyPhiMon");
    //h_A->SetLineColor(kRed);
    //h_A->SetAxisRange(-0.6,0.6,"Y");
    //h_A->Draw();
    //canvasAsy.Print("asyphi.pdf)");
    
    int numBinsX = h_Aplus->GetNbinsX();    //same bins 4 Dplus and Dminus 
    int numBinsY = h_Aplus->GetNbinsY(); 
    int numBinsZ = h_Aplus->GetNbinsZ(); 


    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {
        //double xValue = h_Dplus->GetXaxis()->GetBinCenter(Xbin+1);
        //std::string pdfFileName = "calcPlotAsymphih_z" + std::to_string(xValue) + ".pdf";
        ////this is one pdf for each x value ( which is z )
        //TCanvas canvasphi("c", "Multiplot phih", 1200, 800);
        ////one canvas per pdf, which be divided in 6 (we ll use 5 pads)
        //canvasphi.Divide(3, 2);
        for (int Ybin = 1; Ybin <= numBinsY; Ybin++) {
            ////both X,Y are binned the same for h_Dplus and h_Dminus
            //double YValue = h_Dplus->GetYaxis()->GetBinCenter(Ybin+1);
            //canvasphi.cd(Ybin+ 1);
            //TGraphErrors* graphphih = new TGraphErrors();
            ////create one graph per pad, so per pt2 value. 
            ////following loop allows to have graphs as a fct of Phih
            for (int Zbin = 1; Zbin <= numBinsZ; Zbin++) {
                //same here, Z is binned the same for h_Dplus and h_Dminus
                double ZValue = h_Aplus->GetZaxis()->GetBinCenter(Zbin );   //useless here
                double valueplus = h_Aplus->GetBinContent(Xbin, Ybin, Zbin);
                double valueminus = h_Aminus->GetBinContent(Xbin, Ybin, Zbin);
                double errorplus = (valueplus > 0) ? sqrt(valueplus) : 0.0;
                double errorminus = (valueminus > 0) ? sqrt(valueminus) : 0.0; 
                double AsymPoint =  (valueplus - valueminus) / (valueplus + valueminus) ;   //add ward here ?
                double Err_AsymPoint = (2 * sqrt((valueminus * errorplus) * (valueminus * errorplus) + (valueplus * errorminus) * (valueplus * errorminus))) /((valueplus + valueminus) * (valueplus + valueminus));

                //std::cout << "Asymmetry = " << value << " +- " << error << std::endl;
                //graphphih->SetPoint(Zbin, ZValue, AsymPoint);
                //std::cout << "phih = " << ZValue << std::endl;
                //graphphih->SetPointError(z, 0.0, Err_AsymPoint);
                AsymmMatrix[Xbin-1][Ybin-1][Zbin-1] = AsymPoint;
                errorAsymmMatrix[Xbin-1][Ybin-1][Zbin-1] = Err_AsymPoint;

     
            }
            //graphphih->SetTitle(("Asymmetry (#phi_{h}>_D) , pt2=" + std::to_string(YValue)).c_str());
            //graphphih->GetXaxis()->SetTitle("#phi_{h}");
            //graphphih->GetYaxis()->SetTitle("A(#phi_{h})");
            //graphphih->SetMarkerStyle(20);
            //graphphih->GetYaxis()->SetRangeUser(-1.5, 1.5);            graphphih->Draw("AP");
            //TLine* line = new TLine(graphphih->GetXaxis()->GetXmin(), 0.0, graphphih->GetXaxis()->GetXmax(), 0.0);
            //line->SetLineStyle(2);
            //line->Draw("same");
        }
        //canvasphi.SaveAs(pdfFileName.c_str());
    }
    
}
//if ((valueplus + valueminus) != 0) {
//        double AsymPoint = (valueplus - valueminus) / (valueplus + valueminus);
//        double Err_AsymPoint = (2 * sqrt((valueminus * errorplus) * (valueminus * errorplus) + (valueplus * errorminus) * (valueplus * errorminus))) /
//                               ((valueplus + valueminus) * (valueplus + valueminus));
//
//        // Set the point and error in the TGraphErrors
//        graphphih->SetPoint(Zbin, ZValue, AsymPoint);
//        graphphih->SetPointError(Zbin, 0.0, Err_AsymPoint);
//    } else {
//        // Handle division by zero if the denominator is zero
//        graphphih->SetPoint(Zbin, ZValue, 0);
//        graphphih->SetPointError(Zbin, 0.0, 0.0);
//    }

void sratio::writeAsymmToFile(const std::string& filename){
    //not needed unless output in file wanted for verification  
    std::ofstream outputFile(filename);
    for (int x = 0; x < Cratiobin; ++x) { //cratiobin should be 5 
        double xaxisval = h_Dplus->GetXaxis()->GetBinCenter(x + 1);
        outputFile << "z = " << xaxisval << std::endl;
        for (int y =0; y< Cratiobin; ++y) {
            for (int z = 0; z < Cratiobin; ++z) {
                double value = AsymmMatrix[x][y][z];
                double error = errorAsymmMatrix[x][y][z];
                outputFile << value << " +- " << error << "\t\t";
            }
            outputFile << std::endl;
        }
        outputFile << std::endl << std::endl;
    }
    outputFile.close();
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


void sratio::multiplotSratio( sratio& SratioCu , sratio& SratioC1 , sratio& SratioC2){
    for (int x = 0; x < Constants::Cratiobin; ++x){
        double xValue = h_wD_Sratio->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "multiSrat_x" + std::to_string(xValue) + ".pdf";
        TCanvas canvasSratio("c", "Multiplot Sratio", 1200, 800);
        canvasSratio.Divide(3, 2); 
        for (int y=0; y < Constants::Cratiobin; ++y) {
            TMultiGraph *mg = new TMultiGraph();

            double Q2Value = h_wD_Sratio->GetYaxis()->GetBinCenter(y + 1);
            canvasSratio.cd(y + 1);
            TGraphErrors *graphSratioSn = new TGraphErrors();
            TGraphErrors *graphSratioCu = new TGraphErrors();
            TGraphErrors *graphSratioC1 = new TGraphErrors();
            TGraphErrors *graphSratioC2 = new TGraphErrors();
            for (int z=0; z < Constants::Cratiobin; ++z) {
                double zValue = h_wD_Sratio->GetZaxis()->GetBinCenter(z + 1);
                double valueSn = SratioMatrix[x][y][z];
                double errorSn = errorSratioMatrix[x][y][z];
                double valueCu = SratioCu.SratioMatrix[x][y][z];
                double errorCu = SratioCu.errorSratioMatrix[x][y][z];
                double valueC1 = SratioC1.SratioMatrix[x][y][z];
                double errorC1 = SratioC1.errorSratioMatrix[x][y][z];
                double valueC2 = SratioC2.SratioMatrix[x][y][z];
                double errorC2 = SratioC2.errorSratioMatrix[x][y][z];
                graphSratioSn->SetPoint(z, zValue, valueSn);
                graphSratioSn->SetPointError(z, 0.0, errorSn);
                graphSratioCu->SetPoint(z, zValue+0.01, valueCu);
                graphSratioCu->SetPointError(z+0.01, 0.0, errorCu); 
                graphSratioC1->SetPoint(z, zValue+0.02, valueC1);
                graphSratioC1->SetPointError(z+0.02, 0.0, errorC1);
                graphSratioC2->SetPoint(z, zValue+0.03, valueC2);
                graphSratioC2->SetPointError(z+0.03, 0.0, errorC2);

                
            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << Q2Value;
            std::string formattedQ2Value = ss.str();
            std::string title = "<sin #phi_{h}>_A / <sin #phi_{h}>_D vs z, Q^{2}=" + formattedQ2Value + "GeV^{2}"; 

            graphSratioSn->SetTitle(title.c_str());
            graphSratioSn->GetXaxis()->SetTitle("z");
            graphSratioSn->GetYaxis()->SetTitle("<sin #phi_{h}>_A / <sin #phi_{h}>_D");
            graphSratioSn->GetYaxis()->SetRangeUser(-10.0, 10.0); // Set Y axis range from 0.0 to 2.0

            graphSratioSn->SetMarkerStyle(20);
            graphSratioCu->SetMarkerStyle(20);
            graphSratioC1->SetMarkerStyle(20);
            graphSratioC2->SetMarkerStyle(20);
            //graphCratioOther->Draw("P");
            
            graphSratioSn->SetMarkerColor(kGreen);
            graphSratioCu->SetMarkerColor(kRed);
            graphSratioC1->SetMarkerColor(kBlue);
            graphSratioC2->SetMarkerColor(kBlack);
            //graphCratioThird->Draw("P");

            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graphSratioSn, "Sn", "lp");
            legend->AddEntry(graphSratioCu, "Cu", "lp");
            legend->AddEntry(graphSratioC1, "C2", "lp");
            legend->AddEntry(graphSratioC2, "C1", "lp");

            TLine *line = new TLine(graphSratioSn->GetXaxis()->GetXmin(), 1.0, graphSratioSn->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line

            mg->Add(graphSratioSn);
            mg->Add(graphSratioCu);
            mg->Add(graphSratioC1);
            mg->Add(graphSratioC2);

            mg->SetTitle(("<sin #phi_{h}>_A / <sin #phi_{h}>_D vs z, Q^{2}=" + formattedQ2Value).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("<sin #phi_{h}>_A / <sin #phi_{h}>_D");
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
        canvasSratio.SaveAs(pdfFileName.c_str());
 
    }

}



void sratio::multiSratsimus(sratio& SratioCu ,sratio& SratioC1,sratio& SratioC2){
        for (int x = 0; x < Constants::Cratiobin; ++x){
        double xValue = h_wD_Sratio->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "MultiSratiosim_x" + std::to_string(xValue) + ".pdf";
        TCanvas canvasSratiosim("c", "Multiplot Sratio", 1200, 800);
        canvasSratiosim.Divide(3, 2); 
        for (int y=0; y < Constants::Cratiobin; ++y) {
            TMultiGraph *mg = new TMultiGraph();

            double Q2Value = h_wD_Sratio->GetYaxis()->GetBinCenter(y + 1);
            canvasSratiosim.cd(y + 1);
            TGraphErrors *graphSratioSn = new TGraphErrors();
            TGraphErrors *graphSratioCu = new TGraphErrors();
            TGraphErrors *graphSratioC1 = new TGraphErrors();
            TGraphErrors *graphSratioC2 = new TGraphErrors();
            for (int z=0; z < Constants::Cratiobin; ++z) {
                double zValue = h_wD_Sratio->GetZaxis()->GetBinCenter(z + 1);
                double valueSn = SratioMatrix[x][y][z];
                double errorSn = errorSratioMatrix[x][y][z];
                double valueCu = SratioCu.SratioMatrix[x][y][z];
                double errorCu = SratioCu.errorSratioMatrix[x][y][z];
                double valueC1 = SratioC1.SratioMatrix[x][y][z];
                double errorC1 = SratioC1.errorSratioMatrix[x][y][z];
                double valueC2 = SratioC2.SratioMatrix[x][y][z];
                double errorC2 = SratioC2.errorSratioMatrix[x][y][z];
                graphSratioSn->SetPoint(z, zValue, valueSn);
                graphSratioSn->SetPointError(z, 0.0, errorSn);
                graphSratioCu->SetPoint(z, zValue+0.01, valueCu);
                graphSratioCu->SetPointError(z+0.01, 0.0, errorCu); 
                graphSratioC1->SetPoint(z, zValue+0.02, valueC1);
                graphSratioC1->SetPointError(z+0.02, 0.0, errorC1);
                graphSratioC2->SetPoint(z, zValue+0.03, valueC2);
                graphSratioC2->SetPointError(z+0.03, 0.0, errorC2);

                
            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << Q2Value;
            std::string formattedQ2Value = ss.str();
            std::string title = "<sin #phi_{h}>_A / <sin #phi_{h}>_D vs z, Q^{2}=" + formattedQ2Value + "GeV^{2}"; 

            graphSratioSn->SetTitle(title.c_str());
            graphSratioSn->GetXaxis()->SetTitle("z");
            graphSratioSn->GetYaxis()->SetTitle("<sin #phi_{h}>_A / <sin #phi_{h}>_D");
            graphSratioSn->GetYaxis()->SetRangeUser(-10.0, 10.0); // Set Y axis range from 0.0 to 2.0

            graphSratioSn->SetMarkerStyle(20);
            graphSratioCu->SetMarkerStyle(20);
            graphSratioC1->SetMarkerStyle(20);
            graphSratioC2->SetMarkerStyle(20);
            //graphCratioOther->Draw("P");
            
            graphSratioSn->SetMarkerColor(kGreen);
            graphSratioCu->SetMarkerColor(kRed);
            graphSratioC1->SetMarkerColor(kBlue);
            graphSratioC2->SetMarkerColor(kBlack);
            //graphCratioThird->Draw("P");

            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graphSratioSn, "Sn sim", "lp");
            legend->AddEntry(graphSratioCu, "Cu sim", "lp");
            legend->AddEntry(graphSratioC1, "C1 sim", "lp");
            legend->AddEntry(graphSratioC2, "C2 sim", "lp");

            TLine *line = new TLine(graphSratioSn->GetXaxis()->GetXmin(), 1.0, graphSratioSn->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line

            mg->Add(graphSratioSn);
            mg->Add(graphSratioCu);
            mg->Add(graphSratioC1);
            mg->Add(graphSratioC2);

            mg->SetTitle(("<sin #phi_{h}>_A / <sin #phi_{h}>_D vs z, Q^{2}=" + formattedQ2Value).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("<sin #phi_{h}>_A / <sin #phi_{h}>_D");
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
        canvasSratiosim.SaveAs(pdfFileName.c_str());
 
    }

}




void sratio::DrawMonSinrat(const std::string& outputname){
    TCanvas canvasSmon("c", "Mon Sratio", 1200, 800);
    canvasSmon.Divide(2, 2); 
    canvasSmon.cd(1);
    h_xphiA->Draw("colz");
    canvasSmon.cd(2);
    h_QphiA->Draw("colz");
    canvasSmon.cd(3);
    h_zphiA->Draw("colz");
    canvasSmon.cd(4);
    //h_xQA->Draw("colz");
    h_phiMonA->Draw();

    canvasSmon.Print((outputname + ".pdf(").c_str());
    TCanvas canvasSmon2("c2", "Mon Sratio2", 1200, 800);
    canvasSmon2.Divide(2, 2); 
    canvasSmon2.cd(1);
    h_xphiD->Draw("colz");
    canvasSmon2.cd(2);
    h_QphiD->Draw("colz");
    canvasSmon2.cd(3);
    h_zphiD->Draw("colz");
    canvasSmon2.cd(4);
    //h_xQD->Draw("colz");
    h_phiMonD->Draw();
    canvasSmon2.Print((outputname + ".pdf)").c_str());


}

//void multiSratRGD(sratio&,sratio&,sratio&); //for 4 targets in RGD data (real)
//void multiSratall(sratio&,sratio&,sratio&,sratio&,sratio&,sratio&,sratio&); //for 4 targets in sim and RGD (separation of C1 & C2)
//void multiSratall2(sratio&,sratio&,sratio&,sratio&,sratio&); // three ragets in sim and RGD ( no separation of CxC)
void sratio::ValidateHistograms() {
    std::cout << "Validating histograms for Sratio calculations..." << std::endl;

    // List of 1D histograms to validate
    std::vector<std::pair<TH1F*, std::string>> hist1D = {
        {hSratio_nuD, "hSratio_nuD"},
        {hSratio_nuA, "hSratio_nuA"},
        {hSratio_zD, "hSratio_zD"},
        {hSratio_zA, "hSratio_zA"}
    };

    // List of 3D histograms to validate
    std::vector<std::pair<TH3F*, std::string>> hist3D = {
        {h_wD_Sratio, "h_wD_Sratio"},
        {h_wA_Sratio, "h_wA_Sratio"},
        {h_D_Sratio3D, "h_D_Sratio3D"},
        {h_A_Sratio3D, "h_A_Sratio3D"}
    };

    // Validate and track 1D histograms
    for (const auto& [hist, name] : hist1D) {
        if (hist) {
            std::cout << "Histogram " << name << " has " << hist->GetEntries() << " entries." << std::endl;
            if (hist->GetEntries() == 0) {
                std::cerr << "Warning: Histogram " << name << " is empty!" << std::endl;
            }
        } else {
            std::cerr << "Error: Histogram " << name << " is not initialized!" << std::endl;
        }
    }

    // Validate and track 3D histograms
    for (const auto& [hist, name] : hist3D) {
        if (hist) {
            std::cout << "Histogram " << name << " has " << hist->GetEntries() << " entries." << std::endl;
            if (hist->GetEntries() == 0) {
                std::cerr << "Warning: Histogram " << name << " is empty!" << std::endl;
            }
        } else {
            std::cerr << "Error: Histogram " << name << " is not initialized!" << std::endl;
        }
    }

    // Check consistency of binning between corresponding histograms
    if (h_wD_Sratio && h_wA_Sratio) {
        if (h_wD_Sratio->GetNbinsX() != h_wA_Sratio->GetNbinsX() ||
            h_wD_Sratio->GetNbinsY() != h_wA_Sratio->GetNbinsY() ||
            h_wD_Sratio->GetNbinsZ() != h_wA_Sratio->GetNbinsZ()) {
            std::cerr << "Error: Binning mismatch between h_wD_Sratio and h_wA_Sratio!" << std::endl;
        }
    } else {
        std::cerr << "Error: h_wD_Sratio or h_wA_Sratio is not initialized!" << std::endl;
    }

    // Add further specific validations if needed
    std::cout << "Histogram validation completed." << std::endl;
}

void sratio::LogBinContent() {
    std::cout << "Logging bin content for 3D histograms in Sratio..." << std::endl;

    // Function to log a 3D histogram in a table-like format
    auto logHistogram = [](TH3F* hist, const std::string& name) {
        if (!hist) {
            std::cerr << "Error: Histogram " << name << " is not initialized!" << std::endl;
            return;
        }

        int numBinsX = hist->GetNbinsX();
        int numBinsY = hist->GetNbinsY();
        int numBinsZ = hist->GetNbinsZ();

        std::cout << "BIN CONTENT FOR " << name << " (3D):" << std::endl;

        for (int x = 1; x <= numBinsX; ++x) { // Loop over 'xb' bins
            double xValue = hist->GetXaxis()->GetBinCenter(x);
            std::cout << "xb = " << xValue << " (Layer " << x << "):" << std::endl;

            // Print column headers (z values)
            std::cout << std::setw(10) << "pt2 \\ z";
            for (int y = 1; y <= numBinsY; ++y) {
                double zValue = hist->GetYaxis()->GetBinCenter(y);
                std::cout << std::setw(10) << zValue;
            }
            std::cout << std::endl;

            // Print matrix content for fixed 'xb' layer
            for (int z = 1; z <= numBinsZ; ++z) {
                double pt2Value = hist->GetZaxis()->GetBinCenter(z);
                std::cout << std::setw(10) << pt2Value; // Row header (pt2)

                for (int y = 1; y <= numBinsY; ++y) {
                    double content = hist->GetBinContent(x, y, z);
                    std::cout << std::setw(10) << content; // Bin content
                }
                std::cout << std::endl;
            }
            std::cout << std::endl; // Separate layers for clarity
        }
    };

    // List of 3D histograms to log
    std::vector<std::pair<TH3F*, std::string>> histograms = {
        {h_wD_Sratio, "h_wD_Sratio"},
        {h_wA_Sratio, "h_wA_Sratio"},
        {h_D_Sratio3D, "h_D_Sratio3D"},
        {h_A_Sratio3D, "h_A_Sratio3D"}
    };

    // Log the contents of each histogram
    for (const auto& [hist, name] : histograms) {
        logHistogram(hist, name);
    }

    std::cout << "Finished logging bin contents for 3D histograms in Sratio." << std::endl;
}
