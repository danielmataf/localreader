#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <fstream>
#include <TPDF.h>
#include "Event.h" 
#include "Dpt.h"
#include "CutSet.h"
#include <iomanip> // including <iomanip> for formatting digits for pt2 precision
#include <TMultiGraph.h>
#include <TLegend.h>


deltaptsq::deltaptsq(CutSet cutsD, CutSet cutsA, const std::string& targetName) //: cutsD(cutsD), cutsSn(cutsSn) {
    :  cut2(cutsA),cut1(cutsD),targetName(targetName),
    DptMatrix(Dptbin, std::vector<std::vector<double>>(Dptbin, std::vector<double>(Dptbin, 0.0))),
    errorDptMatrix(Dptbin,  std::vector<std::vector<double>>(Dptbin, std::vector<double>(Dptbin, 0.0))), 
    DptMatrixbis(Dptbin, std::vector<std::vector<double>>(Dptbin, std::vector<double>(Dptbin, 0.0))),
    errorDptMatrixbis(Dptbin,  std::vector<std::vector<double>>(Dptbin, std::vector<double>(Dptbin, 0.0))), 
    hDpt_Q_nu_zD(new TH3F (("Dpt:nu,z,pt2,D_"+targetName).c_str(),("h_Dpt:nu,z,pt2,D_"+targetName).c_str(), Constants::Dptbin_Q, Constants::RcutminQ, Constants::RcutmaxQ, Constants::Dptbin_nu,Constants::numinDpt,Constants::numaxDpt, Constants::Dptbin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    hDpt_Q_nu_zA(new TH3F (("Dpt:nu,z,pt2,A_"+targetName).c_str(),("h_Dpt:nu,z,pt2,A_"+targetName).c_str(), Constants::Dptbin_Q, Constants::RcutminQ, Constants::RcutmaxQ, Constants::Dptbin_nu,Constants::numinDpt,Constants::numaxDpt, Constants::Dptbin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
 
    hDpt_nu_z_pt2A_onlye(new TH3F (("Dpt:Q2, nu,z,A onlye_"+targetName).c_str(),("h_Dpt:Q2, nu,z,A onlye_"+targetName).c_str(), Constants::Dptbin_Q, Constants::RcutminQ, Constants::RcutmaxQ, Constants::Dptbin_nu,numinDpt,numaxDpt,Constants::Dptbin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    hDpt_nu_z_pt2D_onlye(new TH3F (("Dpt:Q2, nu,z,D onlye_"+targetName).c_str(),("h_Dpt:Q2, nu,z,D onlye_"+targetName).c_str(), Constants::Dptbin_Q, Constants::RcutminQ, Constants::RcutmaxQ, Constants::Dptbin_nu,numinDpt,numaxDpt,Constants::Dptbin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    hDpt_nuA(new TH1F (("Dpt:nu_A_"+targetName).c_str(),("h_Dpt:nu_A_"+targetName).c_str(), Constants::Dptbin_nu,Constants::numinDpt,Constants::numaxDpt)),
    hDpt_nuD(new TH1F (("Dpt:nu_D_"+targetName).c_str(),("h_Dpt:nu_D_"+targetName).c_str(), Constants::Dptbin_nu,Constants::numinDpt,Constants::numaxDpt)),

    h_wD_pt(new TH3F (("wD_pt_"+targetName).c_str(), ("h_wD_pt_"+targetName).c_str(), Constants::Dptbin_x  , Constants::xminDpt, Constants::xmaxDpt,Constants::Dptbin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Dptbin_z,Constants::RcutminZ, Constants::RcutmaxZ   )),
    h_wA_pt(new TH3F (("wA_pt_"+targetName).c_str(), ("h_wA_pt_"+targetName).c_str(), Constants::Dptbin_x  , Constants::xminDpt, Constants::xmaxDpt,Constants::Dptbin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Dptbin_z,Constants::RcutminZ, Constants::RcutmaxZ   )),

    h_D_pt3D(new TH3F (("count:D_"+targetName).c_str(), ("h_count:D_"+targetName).c_str(),Constants::Dptbin_x  , Constants::xminDpt, Constants::xmaxDpt,Constants::Dptbin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Dptbin_z,Constants::RcutminZ, Constants::RcutmaxZ    )),
    h_A_pt3D(new TH3F (("count:A_"+targetName).c_str(), ("h_count:A_"+targetName).c_str(),Constants::Dptbin_x  , Constants::xminDpt, Constants::xmaxDpt,Constants::Dptbin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Dptbin_z,Constants::RcutminZ, Constants::RcutmaxZ    )),
    



    h_D_onlypt2(new TH1F (("count:D_pt2_"+targetName).c_str(), ("h_count:D_pt2_"+targetName).c_str(),5  , Constants::RcutminPt2, Constants::RcutmaxPt2    )),
    h_A_onlypt2(new TH1F (("count:A_pt2_"+targetName).c_str(), ("h_count:A_pt2_"+targetName).c_str(),5  , Constants::RcutminPt2, Constants::RcutmaxPt2    )),
    
    h_wD_sqpt2(new TH3F (("wD_sqpt2_"+targetName).c_str(), ("h_wD_sqpt2_"+targetName).c_str(),Constants::Dptbin_x,Constants::xminDpt, Constants::xmaxDpt,Constants::Dptbin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Dptbin_z,Constants::RcutminZ, Constants::RcutmaxZ)),
    h_wA_sqpt2(new TH3F (("wA_sqpt2_"+targetName).c_str(), ("h_wA_sqpt2_"+targetName).c_str(),Constants::Dptbin_x,Constants::xminDpt, Constants::xmaxDpt,Constants::Dptbin_Q,Constants::RcutminQ,Constants::RcutmaxQ,Constants::Dptbin_z,Constants::RcutminZ, Constants::RcutmaxZ)),

    h_z_A_had(new TH1F (("zMonDpt_A_"+targetName).c_str(), ("h_zMonDpt_A_"+targetName).c_str() , 50, 0 , 1)),
    h_z_D_had(new TH1F (("zMonDpt_D_"+targetName).c_str(), ("h_zMonDpt_D_"+targetName).c_str() , 50, 0 , 1)),
    h_pt2_A_had(new TH1F (("ptMonDpt_A_"+targetName).c_str(), ("h_ptMonDpt_A_"+targetName).c_str() , 50, 0 , 3)),
    h_pt2_D_had(new TH1F (("ptMonDpt_D_"+targetName).c_str(), ("h_ptMonDpt_D_"+targetName).c_str() , 50, 0 , 3)) {

    
    //cutd(cutsD),
    //cuta(cutsA);
    }

// ONLY NEED 3 VARIABLES Q,  Z and nu
// Dont use pt nor pt2 as a variable 

void deltaptsq::FillOnlyptandz(const Event& event) {
    int targetType = event.GetTargetType();
    if (targetType == 1 && cut2.PassCutsElectrons(event)==true) {
        for (const Particle& hadron : event.GetHadrons()) {
            if (cut2.PassCutsHadrons(hadron)==true){
                //if (hadron.Getz()>=0.3 && hadron.Getz()<=0.7) {
                h_z_D_had->Fill(hadron.Getz());
                h_pt2_D_had->Fill(hadron.Getpt2());
                //} 
            }
        }
    }
    else if (targetType == 0 && cut1.PassCutsElectrons(event)==true) {
        for (const Particle& hadron : event.GetHadrons()) {
            if (cut1.PassCutsHadrons(hadron)==true){
                //if (hadron.Getz()>=0.3 && hadron.Getz()<=0.7) {
                h_z_A_had->Fill(hadron.Getz());
                h_pt2_A_had->Fill(hadron.Getpt2());
                //}
            }
        }
    }

}

void deltaptsq::DrawOnlyptandz(const std::string& filename) {
    TCanvas canvas("c", "DptMon", 1200, 800);
    canvas.Divide(2, 2);
    canvas.cd(1);
    h_z_A_had->Draw();

    canvas.cd(2);
    h_z_D_had->Draw();
    canvas.cd(3);
    h_pt2_A_had->Draw();
    canvas.cd(4);
    h_pt2_D_had->Draw();
    canvas.SaveAs(filename.c_str());
}




void deltaptsq::FillHistograms(const Event& event) {
    int targetType = event.GetTargetType();
    if (targetType == 0 && cut1.PassCutsElectrons(event)==true) {
        counter_elLD2 ++;
        hDpt_nuD->Fill(event.Getnu());
        //still using nu for reference following Mratio 
        //Q can be used ? (needs to be checked)
        for (const Particle& hadron : event.GetHadrons()) {

            if (cut1.PassCutsHadrons(hadron)==true){
                if (hadron.Getz()>=0.3 && hadron.Getz()<=0.7) {
                //std::cout<<" zvalue LD2 "<<hadron.Getz()<<std::endl;
                h_D_onlypt2->Fill(hadron.Getpt2());
                h_wD_pt->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), hadron.Getpt2());    //3 arguments and the WEIGHT
                h_wD_sqpt2->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), hadron.Getpt2()*hadron.Getpt2());    //3 arguments and the WEIGHT (pt2 squared)
                h_D_pt3D->Fill(event.Getxb(), event.GetQ2(), hadron.Getz());    //3 arguments only counts not weight
                hDpt_Q_nu_zD->Fill(event.GetQ2(), event.Getnu(), hadron.Getz()); // uselesss ? 
                //h_z_A_had->Fill(hadron.Getz());
                //h_pt2_A_had->Fill(hadron.Getpt2());
                }
            }
        }
    }
    else if (targetType == 1 && cut2.PassCutsElectrons(event)==true) {
        counter_elSn++; //counter for only electrons for z
        hDpt_nuA->Fill(event.Getnu());
        for (const Particle& hadron : event.GetHadrons()) {
            if (cut2.PassCutsHadrons(hadron)==true){
                if (hadron.Getz()>=0.3 && hadron.Getz()<=0.7) {
                //std::cout<<" zvalue Sn "<<hadron.Getz()<<std::endl;
                h_A_onlypt2->Fill(hadron.Getpt2());
                h_wA_pt->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), hadron.Getpt2());    //3 arguments and the WEIGHT
                h_wA_sqpt2->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), hadron.Getpt2()*hadron.Getpt2());    //3 arguments and the WEIGHT (pt2 squared)
                h_A_pt3D->Fill(event.Getxb(), event.GetQ2(), hadron.Getz());    //3 arguments only counts not weight
                hDpt_Q_nu_zA->Fill(event.GetQ2(), event.Getnu(), hadron.Getz());    //useless I guess
                //h_z_D_had->Fill(hadron.Getz());
                //h_pt2_D_had->Fill(hadron.Getpt2()); 
                }
            }
        }
    }
}

void deltaptsq::calcDpt(){
    //using here X=nu, Y=z, Z=Q2
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
                        std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << Q2Value;
            std::string formattedQ2Value = ss.str();
            std::string title = "#Delta <p_{t}^{2}> vs z, Q^{2}=" + formattedQ2Value + " GeV^{2}";


            graphDpt->SetTitle((title).c_str());
            graphDpt->GetXaxis()->SetTitle("z");
            graphDpt->GetYaxis()->SetTitle("#Delta <p_{t}^{2}> (GeV^{2})");
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
TH3F* deltaptsq::getHwA() {
    return h_wA_pt;
}
TH3F* deltaptsq::getHwA_sqpt2() {
    return h_wA_sqpt2;
}
TH3F* deltaptsq::getHA3D() {
    return h_A_pt3D;
}

void deltaptsq::calcDptCarbon( deltaptsq& deltaptsqOther){
    int graph_pointnb=0;
    //x= xb, y=Q2, z=z
    TH3F* h_wA2_pt    = deltaptsqOther.getHwA();
    TH3F* h_A2_pt3D   = deltaptsqOther.getHA3D();
    TH3F* h_wA2_sqpt2 = deltaptsqOther.getHwA_sqpt2();
    int counter_3D = 0;
    int numBinsX = h_wA_pt->GetNbinsX();    //same bins 4 Deut and A 
    int numBinsY = h_wA_pt->GetNbinsY(); 
    int numBinsZ = h_wA_pt->GetNbinsZ(); 
    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {  
        for (int Ybin = 1; Ybin <= numBinsY; Ybin++ ){
            for (int Zbin =1; Zbin<=numBinsZ; Zbin++  ){
                double valC2 = h_wA2_pt->GetBinContent(Xbin,Ybin,Zbin);
                double valC1 = h_wA_pt->GetBinContent(Xbin,Ybin,Zbin);
                double countC1 = h_A_pt3D->GetBinContent(Xbin,Ybin,Zbin);
                double countC2 = h_A2_pt3D->GetBinContent(Xbin,Ybin,Zbin);
                double sqvalC1 = h_wA_sqpt2->GetBinContent(Xbin,Ybin,Zbin);
                double sqvalC2 = h_wA2_sqpt2->GetBinContent(Xbin,Ybin,Zbin);
                double wavg_ptC1 = (countC1 != 0) ? valC1/countC1 : 0.0;
                double wavg_ptC2 = (countC2 != 0) ? valC2/countC2 : 0.0;
                double dpt_point = wavg_ptC2 - wavg_ptC1;
                double varianceC1 = (countC1 != 0) ? sqvalC1/countC1 - wavg_ptC1*wavg_ptC1 : 0.0;
                double varianceC2 = (countC2 != 0) ? sqvalC2/countC2 - wavg_ptC2*wavg_ptC2 : 0.0;
                double Err_valC1 = (countC1 != 0) ? sqrt(varianceC1/countC1) : 0.0;
                double Err_valC2 = (countC2 != 0) ? sqrt(varianceC2/countC2) : 0.0;
                double Err_dpt_point = sqrt(Err_valC1*Err_valC1 + Err_valC2*Err_valC2);
                DptMatrixbis[Xbin-1][Ybin-1][Zbin-1] = dpt_point;
                errorDptMatrixbis[Xbin-1][Ybin-1][Zbin-1] = Err_dpt_point;
            }
        }
    }
}

void deltaptsq::multiplotDpt(deltaptsq& dptother, deltaptsq& dptthird){
    for (int x = 0; x < Dptbin; ++x) {
        double xValue = h_wD_pt->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "tripleTargetDpt_x" + std::to_string(xValue) + ".pdf";
        TCanvas canvasDpt("c", "Multiplot Dpt", 1200, 800);
        canvasDpt.Divide(3, 2); 

        for (int y=0; y < Dptbin; ++y) {
            TMultiGraph *mg = new TMultiGraph();

            double Q2Value = h_wD_pt->GetYaxis()->GetBinCenter(y + 1);
            TGraphErrors *graphDpt = new TGraphErrors();
            canvasDpt.cd(y + 1);
            TGraphErrors *graph = new TGraphErrors();
            TGraphErrors *graphOther = new TGraphErrors();
            TGraphErrors *graphThird = new TGraphErrors();
            for (int z=0; z < Dptbin; ++z) {
                double zValue = h_wD_pt->GetZaxis()->GetBinCenter(z + 1);
                double value = DptMatrix[x][y][z];
                double error = errorDptMatrix[x][y][z];
                double valueOther = dptother.getDptMatrix()[x][y][z];
                double errorOther = dptother.getErrorMatrix()[x][y][z];
                double valueThird = dptthird.getDptMatrix()[x][y][z];
                double errorThird = dptthird.getErrorMatrix()[x][y][z];
                //std::cout << "("<< x<<" ,"<<y <<" ,"<< z<< " ) = "<<std::endl;   // value << " +- " << error << "\t\t";

                graphDpt->SetPoint(z, zValue, value);
                graphDpt->SetPointError(z, 0.0, error); 
                graphOther->SetPoint(z, zValue+0.01, valueOther);
                graphOther->SetPointError(z+0.01, 0.0, errorOther);
                graphThird->SetPoint(z, zValue+0.02, valueThird);
                graphThird->SetPointError(z+0.02, 0.0, errorThird);


                //graph_rat->SetPoint(graph_pointnb, value, error);
                //graph_pointnb++;
            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << Q2Value;
            std::string formattedQ2Value = ss.str();
            std::string title = "#Delta <p_{t}^{2}> vs z, Q^{2}=" + formattedQ2Value + " GeV^{2}";


            graphDpt->SetTitle((title).c_str());
            graphDpt->GetXaxis()->SetTitle("z");
            graphDpt->GetYaxis()->SetTitle("#Delta <p_{t}^{2}>");
            graphDpt->SetMarkerStyle(20);
            graphDpt->GetYaxis()->SetRangeUser(-1.0, 1.0); // Set Y axis range from -0.1 to 0.1
            graphDpt->Draw("AP");
            graphDpt->SetMarkerStyle(20);
            graphOther->SetMarkerStyle(20);
            graphThird->SetMarkerStyle(20);
            graph->GetYaxis()->SetRangeUser(0.0, 2.0); // Set Y axis range from 0.0 to 2.0
            //graph->Draw("AP");
            graphOther->SetMarkerColor(kBlue);

            graphDpt->SetMarkerColor(kRed);
            graphThird->SetMarkerColor(kGreen);
            //graphOther->Draw("P");
            
            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graphDpt, "Sn", "lp");
            legend->AddEntry(graphOther, "Cu", "lp");
            legend->AddEntry(graphThird, "C", "lp");

            TLine *line = new TLine(graph->GetXaxis()->GetXmin(), 0.0, graphDpt->GetXaxis()->GetXmax(), 0.0);
            line->SetLineStyle(2); // Dotted line

            mg->Add(graphDpt);
            mg->Add(graphOther);
            mg->Add(graphThird);
            mg->SetTitle((title).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("#Delta <p_{t}^{2}>");

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
        canvasDpt.SaveAs(pdfFileName.c_str());
        
    }
}

void deltaptsq::multiplotDptbis(){
   for (int x = 0; x< Dptbin; ++x){
       double xbValue = h_wA_pt->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "sameTargetR_x" + std::to_string(xbValue) + ".pdf";
        TCanvas canvasDpt("c", "Multiplot Dpt", 1200, 800);
        canvasDpt.Divide(3, 2);
        for (int y=0; y < Dptbin; ++y){
            double Q2value = h_wA_pt->GetYaxis()->GetBinCenter(y + 1);
            canvasDpt.cd(y + 1);
            TGraphErrors *graphDpt = new TGraphErrors();
            for (int z=0; z < Dptbin; ++z) {
                double zValue = h_wA_pt->GetZaxis()->GetBinCenter(z + 1);
                double value = DptMatrixbis[x][y][z];
                double error = errorDptMatrixbis[x][y][z];
                graphDpt->SetPoint(z, zValue, value);
                graphDpt->SetPointError(z, 0.0, error);
            }
            graphDpt->SetTitle(("#Delta <p_{t}^{2}> vs z, Q^{2}=" + std::to_string(Q2value)).c_str());
            graphDpt->GetXaxis()->SetTitle("z");
            graphDpt->GetYaxis()->SetTitle("#Delta <p_{t}^{2}>");
            graphDpt->SetMarkerStyle(20);
            graphDpt->GetYaxis()->SetRangeUser(-1.0, 1.0); // Set Y axis range from -0.1 to 0.1
            graphDpt->Draw("AP");

            TLine *line = new TLine(graphDpt->GetXaxis()->GetXmin(), 0.0, graphDpt->GetXaxis()->GetXmax(), 0.0);
            line->SetLineStyle(2); // Dotted line
            line->Draw("same");

                        TLatex* prelimText = new TLatex();
            prelimText->SetTextSize(0.08);  // Larger text size
            prelimText->SetTextAngle(45);
            prelimText->SetTextColorAlpha(kGray + 1, 0.3);  // Gray color with transparency
            prelimText->SetNDC();
            prelimText->SetTextAlign(22);  // Centered alignment
            prelimText->DrawLatex(0.5, 0.5, "preliminary");

        }
        canvasDpt.SaveAs(pdfFileName.c_str()); 
   }
}  