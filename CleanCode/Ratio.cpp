#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TPad.h>
#include <fstream>
#include <TPDF.h>
#include "Event.h" 
#include "Ratio.h"
#include "CutSet.h"
#include <iomanip> // including <iomanip> for formatting digits for pt2 precision

Ratio::Ratio(CutSet cutsD, CutSet cutsA,const std::string& targetName): //: cutsD(cutsD), cutsSn(cutsSn) {
    targetName(targetName),
    ratMatrix(Rbin, std::vector<std::vector<double>>(Rbin, std::vector<double>(Rbin, 0.0))),
    errorMatrix(Rbin, std::vector<std::vector<double>>(Rbin, std::vector<double>(Rbin, 0.0))),
    ratMatrixbis(Rbin, std::vector<std::vector<double>>(Rbin, std::vector<double>(Rbin, 0.0))),
    errorMatrixbis(Rbin, std::vector<std::vector<double>>(Rbin, std::vector<double>(Rbin, 0.0))),

    //histos after passcuthadrons 
    //h_nu_z_pt2D(new TH3F("nu,z,pt2,D", "histo nu,z,pt2 for D", Constants::Rbin_nu , Constants::Rcutminnu , Constants::Rcutmaxnu , Constants::Rbin_z ,Constants::RcutminZ, Constants::RcutmaxZ, Constants::Rbin_pt2 , Constants::RcutminPt2, Constants::RcutmaxPt2  )),
    //h_nu_z_pt2A(new TH3F("nu,z,pt2,A", "histo nu,z,pt2 for A", Constants::Rbin_nu , Constants::Rcutminnu , Constants::Rcutmaxnu , Constants::Rbin_z ,Constants::RcutminZ, Constants::RcutmaxZ, Constants::Rbin_pt2 , Constants::RcutminPt2, Constants::RcutmaxPt2  )),
    h_nu_z_pt2D(new TH3F(("nu,z,pt2_D_"+targetName).c_str(), ("histo nu,z,pt2 for D"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu , Constants::Rcutmaxnu , Constants::Rbin_z ,Constants::RcutminZ, Constants::RcutmaxZ, Constants::Rbin_pt2 , Constants::RcutminPt2, Constants::RcutmaxPt2  )),
    h_nu_z_pt2A(new TH3F(("nu,z,pt2_A_"+targetName).c_str(), ("histo nu,z,pt2 for A"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu , Constants::Rcutmaxnu , Constants::Rbin_z ,Constants::RcutminZ, Constants::RcutmaxZ, Constants::Rbin_pt2 , Constants::RcutminPt2, Constants::RcutmaxPt2  )),
    //histo after passcutelectrons only e for nu unse only
    h_nu_z_pt2A_onlye(new TH3F(("nu,z,pt2,A onlye_"+targetName).c_str(), ("histo_e nu,z,pt2 for A"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu ,Constants::Rcutmaxnu ,Rbin_z,Constants::RcutminZ, Constants::RcutmaxZ,Rbin_pt2, Constants::RcutminPt2, Constants::RcutmaxPt2  )),
    h_nu_z_pt2D_onlye(new TH3F(("nu,z,pt2,D onlye_"+targetName).c_str(), ("histo_e nu,z,pt2 for A"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu ,Constants::Rcutmaxnu ,Rbin_z,Constants::RcutminZ, Constants::RcutmaxZ,Rbin_pt2, Constants::RcutminPt2, Constants::RcutmaxPt2  )),
    h_nuA(new TH1F(("nu_A"+targetName).c_str(), ("nu_A"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu , Constants::Rcutmaxnu)),
    h_nuD(new TH1F(("nu_D"+targetName).c_str(), ("nu_D"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu , Constants::Rcutmaxnu)),
    h_nu_A_had(new TH1F(("nu_A_had"+targetName).c_str(), ("nu_A_had"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu , Constants::Rcutmaxnu)),
    h_z_A_had(new TH1F(("z_A_had"+targetName).c_str(), ("z_A_had"+targetName).c_str(), Constants::Rbin_z , Constants::RcutminZ , Constants::RcutmaxZ)),
    h_pt2_A_had(new TH1F(("pt2_A_had"+targetName).c_str(), ("pt2_A_had"+targetName).c_str(), Constants::Rbin_pt2 , Constants::RcutminPt2 , Constants::RcutmaxPt2)),
    h_nuC1(new TH1F(("nu_C1"+targetName).c_str(), ("nu_C1"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu , Constants::Rcutmaxnu)), 
    h_nuC2(new TH1F(("nu_C2"+targetName).c_str(), ("nu_C2"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu , Constants::Rcutmaxnu))

    
     {
//Ratio::Ratio(CutSet a)
  
    //cut1 = a;
    cutd = cutsD;
    cuta = cutsA;
    
    //class takes a set of cuts anyway 
    //in order to determine Ratio from the histos with these cuts 
}

//~ add deletes 


void Ratio::FillHistograms(const Event& event) {
        
    int targetType = event.GetTargetType();
        
    //add passcut el onlyu once ten cut on target. inverse to =false
    //then add  ALu o cut needed + pass cut hadrons 


                                                //using a flag for targets 
    //forget condition on taarget, that should directly be a cut !!!! TBD 
    if (targetType == 0 && cutd.PassCutsElectrons(event)==true && cutd.PassCutsDetectors(event)==true) {
        counter_elLD2 ++;
        //set a counter that increases when electroncuts = passed; in order for it to be called when R is  computed in had variables (?) TBD
        h_nuD->Fill(event.Getnu());
        //std::cout << "nu = " << event.Getnu() << std::endl;
        for (const Particle& hadron : event.GetHadrons()) {
            if (cutd.PassCutsHadrons(hadron)==true){
                ////not using the if (==false) return statement 

                h_nu_z_pt2D->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
            //std::cout << "nimporte quoi" << event.Getnu()<< ";" << hadron.Getz()<<","<< hadron.Getpt2() <<std::endl;

            }
        }
    }
    else if (targetType == 1 && cuta.PassCutsElectrons(event)==true && cuta.PassCutsDetectors(event)==true) {
        counter_elSn++; //counter for only electrons for z and pt
        //here change the else if to just else in order to have a generic target 
        h_nuA->Fill(event.Getnu()); //can be plotted just like this 
        //if (targetName == "C1"){
        //    h_nuC1->Fill(event.Getnu());
        //}
        //else if (targetName == "C2"){
        //    h_nuC2->Fill(event.Getnu());
        //}
        for (const Particle& hadron : event.GetHadrons()) {


            if (cuta.PassCutsHadrons(hadron)==true){
                h_nu_z_pt2A->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
                h_nu_A_had->Fill(event.Getnu() );   //these 3 histos were added to monitor CxC self Ratio
                h_z_A_had->Fill(hadron.Getz()); //these 3 histos were added to monitor CxC self Ratio
                h_pt2_A_had->Fill(hadron.Getpt2()); //these 3 histos were added to monitor CxC self Ratio
                //if (targetName == "C1" ) {
                //    h_nu_z_pt2C1->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
                //}
                //else if (targetName == "C2") {
                //    h_nu_z_pt2C2->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
                //}
            //std::cout << "nimporte quoi N" << event.Getnu()<< ";" << hadron.Getz()<<","<< hadron.Getpt2() <<std::endl;

            }
        }
        //Add Here Cu && change cut for Sn here. 1st cut is passelectrons 
    }

    //Add CxC

}

void Ratio::DrawHistos(Ratio& ratioOther ){
    //this is only to monitor the nu histograms for self ratio 
    TCanvas *chR = new TCanvas("c", "c");
    chR->Divide(2, 2);
    TH1F* h_nuA2 = ratioOther.getHNuA();
    TH3F* h_nu_z_pt2A2 = ratioOther.getHNuzptA();
    
    chR->cd(1);
    h_nuA->Draw();
    h_nuA->SetTitle("nu_A(C1)");
    chR->cd(2);
    h_nuA2->Draw();
    h_nuA2->SetTitle("nu_A(C2)");
    chR->cd(3);
    h_nu_z_pt2A->Draw();
    h_nu_z_pt2A->SetTitle("nu_z_pt2_A(C1)");
    chR->cd(4);
    h_nu_z_pt2A2->Draw();
    h_nu_z_pt2A2->SetTitle("nu_z_pt2_A(C2)");


    chR->SaveAs("nu_histos.pdf");
    delete chR;
    //TCanvas *c = new TCanvas("c", "c", 800, 600);
    //h_nu_z_pt2A->Draw();
    //h_nu_z_pt2D->Draw("same");
    //c->SaveAs("nu_z_pt2_histos.pdf");
    //delete c;

}

void Ratio::DrawSelfHistos(Ratio& ratioOther ){
    //this is to draw the histograms used for self ratio
    TCanvas *cSelf = new TCanvas("mon4self", "mon4self");
    cSelf->Divide(2,2);
    cSelf->cd(1);
    h_nuA->Draw();
    h_nuA->SetTitle("nu_C");
    cSelf->cd(2);
    h_nu_A_had->Draw();
    h_nu_A_had->SetTitle("nu_C_had");
    cSelf->cd(3);
    h_z_A_had->Draw();
    h_z_A_had->SetTitle("z_C_had");
    cSelf->cd(4);
    h_pt2_A_had->Draw();
    h_pt2_A_had->SetTitle("pt2_C_had");
    cSelf->SaveAs("self_histos.pdf");


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
        //std :: cout << "val_nuelD = " << val_nuelD << std::endl;  
        double val_nuelA = h_nuA->GetBinContent(Xbin); 
          
        for (int Ybin = 1; Ybin <= numBinsY; Ybin++ ){
            for (int Zbin = 1; Zbin <= numBinsZ; Zbin++ ){

                // loop with 125 values (5 per axis)
                double valD = h_nu_z_pt2D->GetBinContent(Xbin,Ybin,Zbin);   //seems to properly recover value 
                double valA = h_nu_z_pt2A->GetBinContent(Xbin,Ybin,Zbin);
                double interm1nu = (val_nuelA > 0) ? valA / val_nuelA : 0.0;
                double interm2nu = (val_nuelD > 0) ? valD / val_nuelD : 0.0;
                double ratvalue  = (interm2nu > 0) ?interm1nu/interm2nu : 0.0;
                //std::cout << " valeurs = " << val_nuelA << ";  " << val_nuelD << ";  " << interm2nu << " ; "<< interm1nu<< std::endl;
                
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

TH1F* Ratio::getHNuA() {
    return h_nuA;
}
TH3F* Ratio::getHNuzptA() {
    return h_nu_z_pt2A;
}

void Ratio::calcRcarbon( Ratio& ratioOther){
    int graph_pointnb = 0;
    //TH1D* h_nuA2 = dynamic_cast<TH1D*>(ratioOther.h_nuA);
    //TH1D* h_nu_z_pt2A2 = dynamic_cast<TH1D*>(ratioOther.h_nu_z_pt2A);
    TH1F* h_nuA2 = ratioOther.getHNuA();
    TH3F* h_nu_z_pt2A2 = ratioOther.getHNuzptA();
    //IN 3d HISTO x=nu; y= z; z=pt2;    
    int counter_3D = 0;
    int numBinsX = h_nu_z_pt2A->GetNbinsX();    //same bins 4 Deut and A 
    int numBinsY = h_nu_z_pt2A->GetNbinsY(); 
    int numBinsZ = h_nu_z_pt2A->GetNbinsZ(); 
    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {  
        double val_nuelC2 = h_nuA2->GetBinContent(Xbin); 
        double val_nuelC1 = h_nuA->GetBinContent(Xbin); 
        for (int Ybin = 1; Ybin <= numBinsY; Ybin++ ){
            for (int Zbin = 1; Zbin <= numBinsZ; Zbin++ ){
                // loop with 125 values (5 per axis)
                double valC2 = h_nu_z_pt2A2->GetBinContent(Xbin,Ybin,Zbin);   //seems to properly recover value 
                double valC1 = h_nu_z_pt2A->GetBinContent(Xbin,Ybin,Zbin);
                double interm1nu = (val_nuelC1 > 0) ? valC1 / val_nuelC1 : 0.0;
                double interm2nu = (val_nuelC2 > 0) ? valC2 / val_nuelC2 : 0.0;
                double ratvalue  = (interm2nu > 0) ?interm1nu/interm2nu : 0.0;
                double raterr = ratvalue * sqrt(1/valC1 + 1/valC2 + 1/val_nuelC1 + 1/val_nuelC2);
                ratMatrixbis[Xbin - 1][Ybin - 1][Zbin - 1] = ratvalue;
                errorMatrixbis[Xbin - 1][Ybin - 1][Zbin - 1] = raterr;
                //std::cout << "ratvalue = " << ratvalue << std::endl;
                //std::cout << "ratvalue = " << interm1nu <<" / " <<interm2nu <<" = " << ratvalue << "+-"<< raterr<< std::endl;


               //!!!!
                //This doesnt seem to be working as planned. Try a function that plot all the 4 histos used here to see if we are properly recovering them from the other class
                //tested on self deuterium seems to work ok. but Zeros at low pt2 values.
                //suggesting to plot monitoring (1D) histos for nu, z, pt2 for each target (detailed)

            }
        }
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
            graph->SetMarkerStyle(20);
            graph->GetYaxis()->SetRangeUser(0.0, 2.0); // Set Y axis range from 0.0 to 2.0
            graph->Draw("AP");
            
            TLine *line = new TLine(graph->GetXaxis()->GetXmin(), 1.0, graph->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line
            line->Draw("same");
        }

        canvas.SaveAs(pdfFileName.c_str());
    }
}

void Ratio::multiplotR( Ratio& ratioOther, Ratio& ratiothird){
    for (int x = 0; x < Rbin; ++x) {

        double nuValue = h_nu_z_pt2D->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "tripleTargetR_nu" + std::to_string(nuValue) + ".pdf";
        TCanvas canvas("c", "Multiplot double R", 1200, 800);
        canvas.Divide(3, 2); 
        for (int z = 0; z < Rbin; ++z) {
            TMultiGraph *mg = new TMultiGraph();

            double pt2Value = h_nu_z_pt2D->GetZaxis()->GetBinCenter(z + 1);
            canvas.cd(z + 1);
            TGraphErrors *graph = new TGraphErrors();
            TGraphErrors *graphOther = new TGraphErrors();
            TGraphErrors *graphThird = new TGraphErrors();
            
            for (int y = 0; y < Rbin; ++y) {
                double zValue = h_nu_z_pt2D->GetYaxis()->GetBinCenter(y + 1);
                double value = ratMatrix[x][y][z];
                double error = errorMatrix[x][y][z];
                double valueOther = ratioOther.getRatMatrix()[x][y][z];
                double errorOther = ratioOther.getErrorMatrix()[x][y][z];
                double valueThird = ratiothird.getRatMatrix()[x][y][z];
                double errorThird = ratiothird.getErrorMatrix()[x][y][z];

                //std::cout << "Sn= " << value << "; Cu= " << valueOther << std::endl;
                graph->SetPoint(y, zValue, value);
                graph->SetPointError(y, 0.0, error); 
                graphOther->SetPoint(y, zValue+0.01, valueOther);
                graphOther->SetPointError(y+0.01, 0.0, errorOther);
                graphThird->SetPoint(y, zValue+0.02, valueThird);
                graphThird->SetPointError(y+0.02, 0.0, errorThird);
                

            }

            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << pt2Value;
            std::string formattedPt2Value = ss.str();
            //std::cout << "formattedPt2Value = " << formattedPt2Value << std::endl;

            graph->SetTitle(("R vs z, pt2=" + std::to_string(pt2Value)).c_str());
            graph->GetXaxis()->SetTitle("z");
            graph->GetYaxis()->SetTitle("R");
            graph->SetMarkerStyle(20);
            graphOther->SetMarkerStyle(20);
            graphThird->SetMarkerStyle(20);
            graph->GetYaxis()->SetRangeUser(0.0, 2.0); // Set Y axis range from 0.0 to 2.0
            //graph->Draw("AP");
            graphOther->SetMarkerColor(kBlue);

            graph->SetMarkerColor(kRed);
            graphThird->SetMarkerColor(kGreen);
            //graphOther->Draw("P");
            
            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graph, "Sn", "lp");
            legend->AddEntry(graphOther, "Cu", "lp");
            legend->AddEntry(graphThird, "CxC", "lp");

            TLine *line = new TLine(graph->GetXaxis()->GetXmin(), 1.0, graph->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line

            mg->Add(graph);
            mg->Add(graphOther);
            mg->Add(graphThird);
            mg->SetTitle(("R vs z, pt2=" + formattedPt2Value).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("R");

            mg->Draw("APE1");
            legend->Draw("same");
            line->Draw("same");
        
        }

        canvas.SaveAs(pdfFileName.c_str());
    }
}


void Ratio::multiplotR( Ratio& ratioOther, Ratio& ratiothird, Ratio& ratiosimone,   Ratio& ratiosimtwo){
    for (int x = 0; x < Rbin; ++x) {

        double nuValue = h_nu_z_pt2D->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "triplesimTargetR_nu" + std::to_string(nuValue) + ".pdf";
        TCanvas canvas("c", "Multiplot sim double R", 1200, 800);
        canvas.Divide(3, 2); 
        for (int z = 0; z < Rbin; ++z) {
            TMultiGraph *mg = new TMultiGraph();

            double pt2Value = h_nu_z_pt2D->GetZaxis()->GetBinCenter(z + 1);
            canvas.cd(z + 1);
            TGraphErrors *graph = new TGraphErrors();
            TGraphErrors *graphOther = new TGraphErrors();
            TGraphErrors *graphThird = new TGraphErrors();
            TGraphErrors *graphSimone = new TGraphErrors();
            TGraphErrors *graphSimtwo = new TGraphErrors();
            
            for (int y = 0; y < Rbin; ++y) {
                double zValue = h_nu_z_pt2D->GetYaxis()->GetBinCenter(y + 1);
                double value = ratMatrix[x][y][z];
                double error = errorMatrix[x][y][z];
                double valueOther = ratioOther.getRatMatrix()[x][y][z];
                double errorOther = ratioOther.getErrorMatrix()[x][y][z];
                double valueThird = ratiothird.getRatMatrix()[x][y][z];
                double errorThird = ratiothird.getErrorMatrix()[x][y][z];
                double valueSimone = ratiosimone.getRatMatrix()[x][y][z];
                double errorSimone = ratiosimone.getErrorMatrix()[x][y][z];
                double valueSimtwo = ratiosimtwo.getRatMatrix()[x][y][z];
                double errorSimtwo = ratiosimtwo.getErrorMatrix()[x][y][z];

                //std::cout << "Sn= " << value << "; Cu= " << valueOther << std::endl;
                graph->SetPoint(y, zValue, value);
                graph->SetPointError(y, 0.0, error); 
                graphOther->SetPoint(y, zValue+0.01, valueOther);
                graphOther->SetPointError(y+0.01, 0.0, errorOther);
                graphThird->SetPoint(y, zValue+0.02, valueThird);
                graphThird->SetPointError(y+0.02, 0.0, errorThird);
                graphSimone->SetPoint(y, zValue+0.005, valueSimone); //+0.005 to avoid overlap and putting it next its correcponding REC data 
                graphSimone->SetPointError(y+0.005, 0.0, errorSimone);
                graphSimtwo->SetPoint(y, zValue+0.015, valueSimtwo);
                graphSimtwo->SetPointError(y+0.015, 0.0, errorSimtwo);
                

            }

            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << pt2Value;
            std::string formattedPt2Value = ss.str();
            //std::cout << "formattedPt2Value = " << formattedPt2Value << std::endl;

            graph->SetTitle(("R vs z, pt2=" + std::to_string(pt2Value)).c_str());
            graph->GetXaxis()->SetTitle("z");
            graph->GetYaxis()->SetTitle("R");
            graph->SetMarkerStyle(20);
            graphOther->SetMarkerStyle(20);
            graphThird->SetMarkerStyle(20);
            graph->GetYaxis()->SetRangeUser(0.0, 2.0); // Set Y axis range from 0.0 to 2.0
            //graph->Draw("AP");
            graphOther->SetMarkerColor(kBlue);

            graph->SetMarkerColor(kRed);
            graphThird->SetMarkerColor(kGreen);
            //graphOther->Draw("P");
            graphSimone->SetMarkerColor(kRed);
            graphSimtwo->SetMarkerColor(kBlue+1);

            
            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graph, "Sn", "lp");
            legend->AddEntry(graphOther, "Cu", "lp");
            legend->AddEntry(graphThird, "CxC", "lp");
            legend->AddEntry(graphSimone, "Sn (sim)", "lp");
            legend->AddEntry(graphSimtwo, "Cu (sim)", "lp");

            TLine *line = new TLine(graph->GetXaxis()->GetXmin(), 1.0, graph->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line

            mg->Add(graph);
            mg->Add(graphOther);
            mg->Add(graphThird);
            mg->Add(graphSimone);
            mg->Add(graphSimtwo);
            mg->SetTitle(("R vs z, pt2=" + formattedPt2Value).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("R");

            mg->Draw("APE1");
            legend->Draw("same");
            line->Draw("same");
        
        }

        canvas.SaveAs(pdfFileName.c_str());
    }
}

void Ratio::multiplotRbis() {
    for (int x = 0; x < Rbin; ++x) {

        double nuValue = h_nu_z_pt2D->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "sameTargetR_nu" + std::to_string(nuValue) + ".pdf";
        TCanvas canvas("c", "Multiplot R", 1200, 800);
        canvas.Divide(3, 2); 
        for (int z = 0; z < Rbin; ++z) {
            double pt2Value = h_nu_z_pt2D->GetZaxis()->GetBinCenter(z + 1);
            canvas.cd(z + 1);
            TGraphErrors *graph = new TGraphErrors();
            for (int y = 0; y < Rbin; ++y) {
                double zValue = h_nu_z_pt2D->GetYaxis()->GetBinCenter(y + 1);
                double value = ratMatrixbis[x][y][z];
                double error = errorMatrixbis[x][y][z];
                graph->SetPoint(y, zValue, value);
                graph->SetPointError(y, 0.0, error); 
            }
            graph->SetTitle(("R vs z, pt2=" + std::to_string(pt2Value)).c_str());
            graph->GetXaxis()->SetTitle("z");
            graph->GetYaxis()->SetTitle("R");
            graph->SetMarkerStyle(20);
            graph->SetMarkerStyle(20);
            graph->GetYaxis()->SetRangeUser(0.0, 2.0); // Set Y axis range from 0.0 to 2.0
            graph->Draw("AP");
            
            TLine *line = new TLine(graph->GetXaxis()->GetXmin(), 1.0, graph->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line
            line->Draw("same");
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
