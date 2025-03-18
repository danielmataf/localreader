#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TLatex.h>  
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
#include "THnSparse.h"

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
    //Debugging z histos
    h_z_A(new TH1F(("z_A_hadbis"+targetName).c_str(), ("z_A_hadbis"+targetName).c_str(), Constants::Rbin_z , Constants::RcutminZ , Constants::RcutmaxZ)),
    h_z_D(new TH1F(("z_D_hadbis"+targetName).c_str(), ("z_D_hadbis"+targetName).c_str(), Constants::Rbin_z , Constants::RcutminZ , Constants::RcutmaxZ)),
    //histo after passcutelectrons only e for nu unse only
    h_nu_z_pt2A_onlye(new TH3F(("nu,z,pt2,A onlye_"+targetName).c_str(), ("histo_e nu,z,pt2 for A"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu ,Constants::Rcutmaxnu ,Rbin_z,Constants::RcutminZ, Constants::RcutmaxZ,Rbin_pt2, Constants::RcutminPt2, Constants::RcutmaxPt2  )),
    h_nu_z_pt2D_onlye(new TH3F(("nu,z,pt2,D onlye_"+targetName).c_str(), ("histo_e nu,z,pt2 for A"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu ,Constants::Rcutmaxnu ,Rbin_z,Constants::RcutminZ, Constants::RcutmaxZ,Rbin_pt2, Constants::RcutminPt2, Constants::RcutmaxPt2  )),
    h_nuA(new TH1F(("nu_A"+targetName).c_str(), ("nu_A"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu , Constants::Rcutmaxnu)),
    h_nuD(new TH1F(("nu_D"+targetName).c_str(), ("nu_D"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu , Constants::Rcutmaxnu)),
    h_nu_A_had(new TH1F(("nu_A_had"+targetName).c_str(), ("nu_A_had"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu , Constants::Rcutmaxnu)),
    h_z_A_had(new TH1F(("z_A_had"+targetName).c_str(), ("z_A_had"+targetName).c_str(), Constants::Rbin_z , Constants::RcutminZ , Constants::RcutmaxZ)),
    h_nu_D_had(new TH1F(("nu_D_had"+targetName).c_str(), ("nu_D_had"+targetName).c_str(), Constants::Rbin_nu , Constants::Rcutminnu , Constants::Rcutmaxnu)),
    h_pt2_D_had(new TH1F(("pt2_D_had"+targetName).c_str(), ("pt2_D_had"+targetName).c_str(), Constants::Rbin_pt2 , Constants::RcutminPt2 , Constants::RcutmaxPt2)),
    h_pt2_A_had(new TH1F(("pt2_A_had"+targetName).c_str(), ("pt2_A_had"+targetName).c_str(), Constants::Rbin_pt2 , Constants::RcutminPt2 , Constants::RcutmaxPt2)),
    //double binEdges[7];  // 6 bins means 7 edges
    
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
//~ add deletes complete
Ratio::~Ratio() {
    delete h_nu_z_pt2D;
    delete h_nu_z_pt2A;
    delete h_nu_z_pt2D_onlye;
    delete h_nu_z_pt2A_onlye;
    delete h_nuA;
    delete h_nuD;
    delete h_nu_A_had;
    delete h_z_A_had;
    delete h_pt2_A_had;
    delete graph_rat;
}


void Ratio::FillHistograms(const Event& event) {
        
    int targetType = event.GetTargetType();
        
    //add passcut el onlyu once ten cut on target. inverse to =false
    //then add  ALu o cut needed + pass cut hadrons 

                                                //using a flag for targets 
    //forget condition on taarget, that should directly be a cut !!!! TBD 
    if (targetType == 0 && cutd.PassCutsElectrons(event)==true && cutd.PassCutsDetectors(event)==true ) { 
        counter_elLD2 ++;
        //set a counter that increases when electroncuts = passed; in order for it to be called when R is  computed in had variables (?) TBD
        h_nuD->Fill(event.Getnu());
        //h_3D_D_e->Fill(event.GetQ2(), event.Getxb(), event.Getnu());   //filling a 3D histo for 5D w/ only ele vars
        //std::cout << "nu = " << event.Getnu() << std::endl; /bump
        for (const Particle& hadron : event.GetHadrons()) {
            if (cutd.PassCutsHadrons(hadron)==true){
                //if (hadron.GetPID() == Constants::PION_PLUS_PID  ){
                ////not using the if (==false) return statement 
                h_z_D->Fill(hadron.Getz()); //debugging Ybins jan25
                h_nu_z_pt2D->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
                h_nu_D_had->Fill(event.Getnu() );   //these  histos  added to to track binning switch
                h_pt2_D_had->Fill(hadron.Getpt2()); //these  histos  added to to track binning switch
                //h_5D_D_had->Fill(event.GetQ2(), event.Getxb(), event.Getnu(), hadron.Getz(), hadron.Getpt2());    //filling a 5D histo for 5D calc inside hadron loop

            //std::cout << "nimporte quoi" << event.Getnu()<< ";" << hadron.Getz()<<","<< hadron.Getpt2() <<std::endl;
                //}
            }
        }
    }
    else if (targetType == 1 && cuta.PassCutsElectrons(event)==true && cuta.PassCutsDetectors(event)==true) {
        counter_elSn++; //counter for only electrons for z and pt
        //here change the else if to just else in order to have a generic target 
        h_nuA->Fill(event.Getnu()); //can be plotted just like this 
        //h_3D_A_e->Fill(event.GetQ2(), event.Getxb(), event.Getnu());   //filling a 3D histo for 5D w/ only ele vars
        //if (targetName == "C1"){
        //    h_nuC1->Fill(event.Getnu());
        //}
        //else if (targetName == "C2"){
        //    h_nuC2->Fill(event.Getnu());
        //}
        for (const Particle& hadron : event.GetHadrons()) {


            if (cuta.PassCutsHadrons(hadron)==true){
                //if (hadron.GetPID() == Constants::PION_PLUS_PID  ){
                h_nu_z_pt2A->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
                h_nu_A_had->Fill(event.Getnu() );   //these 3 histos were added to monitor CxC self Ratio
                h_z_A_had->Fill(hadron.Getz()); //these 3 histos were added to monitor CxC self Ratio
                h_pt2_A_had->Fill(hadron.Getpt2()); //these 3 histos were added to monitor CxC self Ratio
                h_z_A->Fill(hadron.Getz()); //debugging Ybins jan25
                //h_5D_A_had->Fill(event.GetQ2(), event.Getxb(), event.Getnu(), hadron.Getz(), hadron.Getpt2());    //filling a 5D histo for 5D calc inside hadron loop
                //if (targetName == "C1" ) {
                //    h_nu_z_pt2C1->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
                //}
                //else if (targetName == "C2") {
                //    h_nu_z_pt2C2->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
                //}
            //std::cout << "nimporte quoi N" << event.Getnu()<< ";" << hadron.Getz()<<","<< hadron.Getpt2() <<std::endl;
                //}
            }
        }
        //Add Here Cu && change cut for Sn here. 1st cut is passelectrons 
    }
    //add below target type 2 in order to have CxC with vertex separation then including
    //make a new fct in cutset for pass cut electron w/o vertex, then add a fct vertex cut taht get vzmin,max for both C1 et C2 
    //else if (targetType == 2 && cuta.PassCutsElectrons(event)==true && cuta.PassCutsDetectors(event)==true) {
    //    h_nuA->Fill(event.Getnu()); //can be plotted just like this 
    //    for (const Particle& hadron : event.GetHadrons()) {
    //    }
    //}
    //Add CxC

}

void Ratio::DrawHistos(Ratio& ratioOther ){
    //this is only to monitor the nu histograms for self ratio 
    TCanvas *chR = new TCanvas("c", "c");
    chR->Divide(3, 3);
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

    chR->cd(5);
    h_z_D->Draw();
    h_z_D->SetTitle("check z bins D");
    chR->cd(6);
    h_z_A->Draw();
    h_z_A->SetTitle("check z bins D");

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
    cSelf->Divide(3,3);
    cSelf->cd(1);
    h_nuA->SetTitle("nu_C");
    //h_nuD->SetLineColor(kRed);
    //h_nuD->Draw();
    h_nuA->Draw();

    cSelf->cd(2);
    h_nu_A_had->Draw();
    h_nu_A_had->SetTitle("nu_C_had");

    cSelf->cd(3);
    h_z_A_had->Draw();
    h_z_A_had->SetTitle("z_C_had");
    
    cSelf->cd(4);
    h_pt2_A_had->Draw();
    h_pt2_A_had->SetTitle("pt2_C_had");

    cSelf->cd(5);
    h_nuD->SetTitle("nu_D");
    //h_nuD->SetLineColor(kRed);
    //h_nuD->Draw();
    h_nuD->Draw();

    cSelf->cd(6);
    h_nu_D_had->Draw();
    h_nu_D_had->SetTitle("nu_D_had");

    cSelf->cd(7);
    h_z_D->Draw();
    h_z_D->SetTitle("z_D_had");
    
    cSelf->cd(8);
    h_pt2_D_had->Draw();
    h_pt2_D_had->SetTitle("pt2_D_had");


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
    std::cout << "numBinsX = " << numBinsX << std::endl;
    std::cout << "numBinsY = " << numBinsY << std::endl;
    std::cout << "numBinsZ = " << numBinsZ << std::endl;
    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {  
        double val_nuelD = h_nuD->GetBinContent(Xbin); 
        //std :: cout << "val_nuelD = " << val_nuelD << std::endl;  
        double val_nuelA = h_nuA->GetBinContent(Xbin); 
          
        for (int Ybin = 1; Ybin <= numBinsY; Ybin++ ){
            for (int Zbin = 1; Zbin <= numBinsZ; Zbin++ ){
                counter_3D ++;
                std ::cout << "counter_3D = " << counter_3D << std::endl;

                // loop with 125 values (5 per axis)
                // if we use with bins = 6 then a total of 216 for counter_3D
                double valD = h_nu_z_pt2D->GetBinContent(Xbin,Ybin,Zbin);   //seems to properly recover value 
                double valA = h_nu_z_pt2A->GetBinContent(Xbin,Ybin,Zbin);
                double interm1nu = (val_nuelA > 0) ? valA / val_nuelA : 0.0;
                double interm2nu = (val_nuelD > 0) ? valD / val_nuelD : 0.0;
                double ratvalue  = (interm2nu > 0) ?interm1nu/interm2nu : 0.0;
                //std::cout << " valeurs = " << val_nuelA << ";  " << val_nuelD << ";  " << interm2nu << " ; "<< interm1nu<< std::endl;
                //std::cout << " ratvalue = " << ratvalue << "= interm1nu/interm2nu =  " << interm1nu << " /   " << interm2nu << std::endl;
                //std::cout << " interm1nu = " << interm1nu << "= valA / val_nuelA =  " << valA  << " /   " << val_nuelA << std::endl;
                //std::cout << " interm2nu = " << interm2nu << "= valD / val_nuelD =  " << valD  << " /   " << val_nuelD << std::endl;
                //std::cout << " valA = " << std::endl;
                //std::cout << "          " << std::endl;
                //
                double raterr = ratvalue * sqrt(1/valA + 1/valD + 1/val_nuelA + 1/val_nuelD);
                //std::cout << "ratvalue = " << ratvalue << "+-"<< raterr<< std::endl;
                ratMatrix[Xbin - 1][Ybin - 1][Zbin - 1] = ratvalue;
                //std::cout << "Xbin = "<< Xbin << " Ybin = " << Ybin << " Zbin = " << Zbin << std::endl;

                errorMatrix[Xbin - 1][Ybin - 1][Zbin - 1] = raterr;
                //std::cout << "raterr = " << raterr << std::endl;
            }
        }

    //    Rq_err= interm3 * sqrt(1/nu_A_had + 1/nu_D_had + 1/nu_A_had + 1/nu_D_had);
    //    R_v->SetPoint(bin-1, x_axis, interm3 );
	//	//	
    }   
    //std::cout << "total counter_3D = " << counter_3D << std::endl;
}

/*
void calcBinEdges(double minVal, double maxVal, int nBins, double binEdges[]) {
    //using cutlow and cut high of whatever variable then the binnbr and the edges desired (nbins+1 ? ) 2
    double binWidth = (maxVal - minVal) / nBins;  //bin width
    for (int i = 0; i <= nBins; i++) {
        binEdges[i] = minVal + i * binWidth;  //min/max
    }
}

void Ratio::calcRin5D() {
    //loop on Q2, xB, and nu  for 
    for (int xBin = 1; xBin <= Constants::Rbin_nu; xBin++) {
        for (int q2Bin = 1; q2Bin <= Constants::Rbin_nu; q2Bin++) {
            for (int nuBin = 1; nuBin <= Constants::Rbin_nu; nuBin++) {
                //get electron counts from TH3F histograms
                double NA_e = h_3D_A_e->GetBinContent(xBin, q2Bin, nuBin);
                double ND_e = h_3D_D_e->GetBinContent(xBin, q2Bin, nuBin);

                //ward
                if (NA_e == 0 || ND_e == 0) continue;

                //loop 
                for (int zBin = 1; zBin <= Constants::Rbin_nu; zBin++) {
                    for (int pt2Bin = 1; pt2Bin <= Constants::Rbin_nu; pt2Bin++) {
                        double NA_had = h_5D_A_had->GetBinContent(xBin, q2Bin, nuBin, zBin, pt2Bin);
                        double ND_had = h_5D_D_had->GetBinContent(xBin, q2Bin, nuBin, zBin, pt2Bin);
                        //ward
                        if (ND_had == 0) {
                            ratMatrix5D[xBin - 1][q2Bin - 1][nuBin - 1][zBin - 1][pt2Bin - 1] = 0.0;
                            errorMatrix5D[xBin - 1][q2Bin - 1][nuBin - 1][zBin - 1][pt2Bin - 1] = 0.0;
                            continue;
                        }

                        double ratio = (NA_had / NA_e) / (ND_had / ND_e);
                        ratMatrix5D[xBin - 1][q2Bin - 1][nuBin - 1][zBin - 1][pt2Bin - 1] = ratio;

                        double error = ratio * sqrt(
                            (1.0 / NA_had) + (1.0 / ND_had) + (1.0 / NA_e) + (1.0 / ND_e)
                        );
                        errorMatrix5D[xBin - 1][q2Bin - 1][nuBin - 1][zBin - 1][pt2Bin - 1] = error;
                    }
                }
            }
        }
    }
}
*/


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


void Ratio::calcRwCC( Ratio& ratioOther){
    int graph_pointnb = 0;
    //TH1D* h_nuA2 = dynamic_cast<TH1D*>(ratioOther.h_nuA);
    //TH1D* h_nu_z_pt2A2 = dynamic_cast<TH1D*>(ratioOther.h_nu_z_pt2A);
    TH1F* h_nuA2 = ratioOther.getHNuA();
    TH3F* h_nu_z_pt2A2 = ratioOther.getHNuzptA();
    //h_nuA.Add(ratioOther.getHNuA);
    //h_nu_z_pt2A.Add(ratioOther.getHNuzptA);
    //h_nu_z_pt2A.Add(ratioOther.get)   //need to figure this one out 
    //h_nu_z_pt2A.Add(ratioOther.get)   //need to figure this one out 
    //h_nu_z_pt2A.Add(ratioOther.get)   //need to figure this one out 
    //IN 3d HISTO x=nu; y= z; z=pt2;    
    int counter_3D = 0;
    int numBinsX = h_nu_z_pt2A->GetNbinsX();    //same bins 4 Deut and A 
    int numBinsY = h_nu_z_pt2A->GetNbinsY(); 
    int numBinsZ = h_nu_z_pt2A->GetNbinsZ(); 
//    for (int Xbin = 1; Xbin <= numBinsX; Xbin++) {  
//        double val_nuelC2 = h_nuA2->GetBinContent(Xbin); 
//        double val_nuelC1 = h_nuA->GetBinContent(Xbin); 
//        for (int Ybin = 1; Ybin <= numBinsY; Ybin++ ){
//            for (int Zbin = 1; Zbin <= numBinsZ; Zbin++ ){
//                // loop with 125 values (5 per axis)
//                double valC2 = h_nu_z_pt2A2->GetBinContent(Xbin,Ybin,Zbin);   //seems to properly recover value 
//                double valC1 = h_nu_z_pt2A->GetBinContent(Xbin,Ybin,Zbin);
//                double interm1nu = (val_nuelC1 > 0) ? valC1 / val_nuelC1 : 0.0;
//                double interm2nu = (val_nuelC2 > 0) ? valC2 / val_nuelC2 : 0.0;
//                double ratvalue  = (interm2nu > 0) ?interm1nu/interm2nu : 0.0;
//                double raterr = ratvalue * sqrt(1/valC1 + 1/valC2 + 1/val_nuelC1 + 1/val_nuelC2);
//                ratMatrixbis[Xbin - 1][Ybin - 1][Zbin - 1] = ratvalue;
//                errorMatrixbis[Xbin - 1][Ybin - 1][Zbin - 1] = raterr;
//
//            }
//        }
//    }   
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
    //this function pltos only one target.
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
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << pt2Value;
            std::string formattedPt2Value = ss.str();
            //std::cout << "formattedPt2Value = " << formattedPt2Value << std::endl;
            std::string title = "R vs z, p_{t}^{2}=" + formattedPt2Value + " GeV^{2}";

            graph->SetTitle((title).c_str() );
            graph->GetXaxis()->SetTitle("z");
            graph->GetYaxis()->SetTitle("R");
            graph->SetMarkerStyle(20);
            graph->GetYaxis()->SetRangeUser(0, 0.5);
            graph->Draw("AP");

            

            TLine *line = new TLine(graph->GetXaxis()->GetXmin(), 1.0, graph->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2);
            line->Draw("same");

            delete line;  // Clean up
        }

        canvas.SaveAs(pdfFileName.c_str());
    }
}
void Ratio::multiplotR(Ratio& ratioSecond) {
    // THIS FUNCTION COMPARES 2 TARGETS, USUALLY C AND SN
    for (int x = 0; x < Rbin; ++x) {
        double nuValue = h_nu_z_pt2D->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "doubleTargetR_nu" + std::to_string(nuValue) + ".pdf";
        TCanvas canvas("c", "Multiplot R", 1200, 800);
        canvas.Divide(3, 2);

        for (int z = 0; z < Rbin; ++z) {
            TMultiGraph *mg = new TMultiGraph();
            double pt2Value = h_nu_z_pt2D->GetZaxis()->GetBinCenter(z + 1);
            canvas.cd(z + 1);
            
            TGraphErrors *graphC = new TGraphErrors(); // Carbon (C)
            TGraphErrors *graphSn = new TGraphErrors(); // Tin (Sn)
            
            for (int y = 0; y < Rbin; ++y) {
                double zValue = h_nu_z_pt2D->GetYaxis()->GetBinCenter(y + 1);
                double valueC = ratMatrix[x][y][z];
                double errorC = errorMatrix[x][y][z];
                double valueSn = ratioSecond.getRatMatrix()[x][y][z];
                double errorSn = ratioSecond.getErrorMatrix()[x][y][z];
                
                graphC->SetPoint(y, zValue, valueC);
                graphC->SetPointError(y, 0.0, errorC);
                graphSn->SetPoint(y, zValue + 0.01, valueSn);
                graphSn->SetPointError(y + 0.01, 0.0, errorSn);
            }
            
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << pt2Value;
            std::string formattedPt2Value = ss.str();
            std::string title = "R vs z, p_{t}^{2}=" + formattedPt2Value + " GeV^{2}";
            
            graphC->SetTitle(title.c_str());
            graphC->GetXaxis()->SetTitle("z");
            graphC->GetYaxis()->SetTitle("R");
            graphC->SetMarkerStyle(20);
            graphC->SetMarkerColor(kBlack);
            
            graphSn->SetMarkerStyle(20);
            graphSn->SetMarkerColor(kOrange);
            
            // Legend: Smaller size and bottom-left position
            TLegend *legend = new TLegend(0.15, 0.15, 0.35, 0.30); // Bottom-left corner
            legend->SetTextSize(0.03);
            legend->SetBorderSize(0);  // No border
            legend->SetFillStyle(0);   // Transparent background
            legend->AddEntry(graphC, "C", "lp");
            legend->AddEntry(graphSn, "Sn", "lp");

            TLine *line = new TLine(graphC->GetXaxis()->GetXmin(), 1.0, graphC->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line

            mg->Add(graphC);
            mg->Add(graphSn);
            mg->SetTitle(title.c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("R");
            mg->GetYaxis()->SetRangeUser(0.0, 1.5);
            
            mg->Draw("APE1");
            legend->Draw("same");
            line->Draw("same");

            TLatex* prelimText = new TLatex();
            prelimText->SetTextSize(0.08);
            prelimText->SetTextAngle(45);
            prelimText->SetTextColorAlpha(kGray + 1, 0.3);
            prelimText->SetNDC();
            prelimText->SetTextAlign(22);
            prelimText->DrawLatex(0.5, 0.5, "very preliminary");
        }
        
        canvas.SaveAs(pdfFileName.c_str());
    }
}
void Ratio::multiplotR(Ratio& ratioOther, Ratio& ratiothird) {
    //THIS FUNCTION COMPARES 3 TARGETS, USUALLY C, CU, AND TIN
    for (int x = 0; x < Rbin; ++x) {
        double nuValue = h_nu_z_pt2D->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "tripleTargetR_nu" + std::to_string(nuValue) + ".pdf";
        TCanvas canvas("c", "Multiplot double R", 1200, 800);
        canvas.Divide(3, 2);

        for (int z = 0; z < Rbin; ++z) {
            double pt2Value = h_nu_z_pt2D->GetZaxis()->GetBinCenter(z + 1);
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

                graph->SetPoint(y, zValue, value);
                graph->SetPointError(y, 0.0, error);
                graphOther->SetPoint(y, zValue + 0.01, valueOther);
                graphOther->SetPointError(y, 0.0, errorOther);
                graphThird->SetPoint(y, zValue + 0.02, valueThird);
                graphThird->SetPointError(y, 0.0, errorThird);
            }

            graph->GetXaxis()->SetTitle("z");
            graph->GetYaxis()->SetTitle("R");
            graph->SetMarkerStyle(20);
            graphOther->SetMarkerStyle(20);
            graphThird->SetMarkerStyle(20);
            graph->GetYaxis()->SetRangeUser(0.0, 2.0);

            graph->SetMarkerColor(kOrange);
            graphOther->SetMarkerColor(kGreen);
            graphThird->SetMarkerColor(kBlack);

            // Create individual canvas for PNG output
            TCanvas pngCanvas("pngCanvas", "", 800, 600);
            pngCanvas.cd();

            TMultiGraph *mg = new TMultiGraph();
            mg->Add(graph);
            mg->Add(graphOther);
            mg->Add(graphThird);
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("R");
            mg->GetXaxis()->SetTitleSize(0.0535);  // Adjusts size of "z"
            mg->GetYaxis()->SetTitleSize(0.0535);  // Adjusts size of "R"
            mg->GetYaxis()->SetRangeUser(0.0, 1.5);

            mg->Draw("APE1");

            // Adjusted legend position and size
            TLegend *legend = new TLegend(0.15, 0.15, 0.35, 0.30);
            legend->SetTextSize(0.04);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->SetTextFont(42);
            legend->AddEntry(graph, "Sn", "lp");
            legend->AddEntry(graphOther, "Cu", "lp");
            legend->AddEntry(graphThird, "C", "lp");
            legend->Draw("same");

            // Nu and Pt² values placed to the **right of the legend**
            TLatex text;
            text.SetTextSize(0.045);
            text.DrawLatexNDC(0.40, 0.22, Form("%.2f < #nu < %.2f", nuValue - 0.8, nuValue + 0.8));
text.DrawLatexNDC(0.40, 0.17, Form("%.2f < p_{t}^{2} < %.2f", pt2Value - 0.2, pt2Value + 0.2));


            TLine *line = new TLine(graph->GetXaxis()->GetXmin(), 1.0, graph->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2);
            line->Draw("same");

            // Save PNG output for each individual plot
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << pt2Value;
            std::string formattedPt2Value = ss.str();
            std::string pngFileName = "tripleTargetR_nu" + std::to_string(nuValue) + "_pt2_" + formattedPt2Value + ".png";
            pngCanvas.SaveAs(pngFileName.c_str());

            // Now draw the same plot inside the PDF layout
            canvas.cd(z + 1);
            mg->Draw("APE1");
            legend->Draw("same");
            text.DrawLatexNDC(0.40, 0.22, Form("#nu = %.2f", nuValue));
            text.DrawLatexNDC(0.40, 0.17, Form("p_{t}^{2} = %.2f GeV^{2}", pt2Value));
            line->Draw("same");
        }

        // Save the multi-plot PDF with 6 graphs
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
            legend->AddEntry(graphThird, "C", "lp");
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
            mg->GetYaxis()->SetRangeUser(-0.5, 1.5); // 

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

        canvas.SaveAs(pdfFileName.c_str());
    }
}

void Ratio::multiplotRbis() {       //this is for two carbons
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
            TLatex* prelimText = new TLatex();
            prelimText->SetTextSize(0.08);  // Larger text size
            prelimText->SetTextAngle(45);
            prelimText->SetTextColorAlpha(kGray + 1, 0.3);  // Gray color with transparency
            prelimText->SetNDC();
            prelimText->SetTextAlign(22);  // Centered alignment
            prelimText->DrawLatex(0.5, 0.5, "preliminary");
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




void Ratio::multiRsimus(  Ratio& ratioCusim, Ratio& ratioC1sim,   Ratio& ratioC2sim){
    //fct for simus only arguments only three ratios Cu, C1, C2. Third is included in the fct-class
    for (int x = 0; x < Rbin; ++x) {
//
        double nuValue = h_nu_z_pt2D->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "AllRsimu_nu" + std::to_string(nuValue) + ".pdf";
        TCanvas canvas("c", "Multiplot sim and four R", 1200, 800);
        canvas.Divide(3, 2); 
        for (int z = 0; z < Rbin; ++z) {
            TMultiGraph *mg = new TMultiGraph();
//
            double pt2Value = h_nu_z_pt2D->GetZaxis()->GetBinCenter(z + 1);
            canvas.cd(z + 1);
            TGraphErrors *graphSn = new TGraphErrors();
            TGraphErrors *graphCu = new TGraphErrors();
            TGraphErrors *graphC1 = new TGraphErrors();
            TGraphErrors *graphC2 = new TGraphErrors();
//            
            for (int y = 0; y < Rbin; ++y) {
                double zValue = h_nu_z_pt2D->GetYaxis()->GetBinCenter(y + 1);
                double valueSn = ratMatrix[x][y][z];
                double errorSn = errorMatrix[x][y][z];
                double valueCu = ratioCusim.getRatMatrix()[x][y][z];
                double errorCu = ratioCusim.getErrorMatrix()[x][y][z];
                double valueC1 = ratioC1sim.getRatMatrix()[x][y][z];
                double errorC1 = ratioC1sim.getErrorMatrix()[x][y][z];
                double valueC2 = ratioC2sim.getRatMatrix()[x][y][z];
                double errorC2 = ratioC2sim.getErrorMatrix()[x][y][z];
                //double valueSimtwo = ratiosimtwo.getRatMatrix()[x][y][z];
                //double errorSimtwo = ratiosimtwo.getErrorMatrix()[x][y][z];
//
                //std::cout << "Sn= " << value << "; Cu= " << valueOther << std::endl;
                graphSn->SetPoint(y, zValue,   valueSn);
                graphSn->SetPointError(y, 0.0, errorSn); 
                graphCu->SetPoint(y, zValue+0.01,   valueCu);
                graphCu->SetPointError(y+0.01, 0.0, errorCu);
                graphC1->SetPoint(y, zValue+0.02,   valueC1);
                graphC1->SetPointError(y+0.02, 0.0, errorC1);
                graphC2->SetPoint(y, zValue+0.03,   valueC2); //+0.005 to avoid overlap and putting it next its correcponding REC data 
                graphC2->SetPointError(y+0.03, 0.0, errorC2);
                //graphSimtwo->SetPoint(y, zValue+0.015, valueSimtwo);
                //graphSimtwo->SetPointError(y+0.015, 0.0, errorSimtwo);
//                
//
            }
//
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << pt2Value;
            std::string formattedPt2Value = ss.str();
            //std::cout << "formattedPt2Value = " << formattedPt2Value << std::endl;
//
            graphSn->SetTitle(("R vs z, pt2=" + std::to_string(pt2Value)).c_str());
            graphSn->GetXaxis()->SetTitle("z");
            graphSn->GetYaxis()->SetTitle("R");
            graphSn->SetMarkerStyle(20);
            graphCu->SetMarkerStyle(20);
            graphC1->SetMarkerStyle(20);
            graphC2->SetMarkerStyle(20);
            graphSn->GetYaxis()->SetRangeUser(0.0, 2.0); // Set Y axis range from 0.0 to 2.0
            //graph->Draw("AP");
            graphCu->SetMarkerColor(kRed);
//
            graphSn->SetMarkerColor(kGreen);
            graphC1->SetMarkerColor(kBlue+1);
            graphC2->SetMarkerColor(kBlack);
            //graphOther->Draw("P");
//
//            
            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graphSn, "Sn sim", "lp");
            legend->AddEntry(graphCu, "Cu sim", "lp");
            legend->AddEntry(graphC1, "C1 sim", "lp");
            legend->AddEntry(graphC2, "C2 sim", "lp");
//
            TLine *line = new TLine(graphSn->GetXaxis()->GetXmin(), 1.0, graphSn->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line
//
            mg->Add(graphSn);
            mg->Add(graphCu);
            mg->Add(graphC1);
            mg->Add(graphC2);
            mg->SetTitle(("R vs z, pt^{2}=" + formattedPt2Value).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("R");
//
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
        canvas.SaveAs(pdfFileName.c_str());
    }
}

void Ratio::multiRtrue(  Ratio& ratioCu, Ratio& ratioC1,   Ratio& ratioC2){
    //fct only for true. arguments only three ratios Cu, C1, C2. Third is included in the fct-class
    for (int x = 0; x < Rbin; ++x) {
        double nuValue = h_nu_z_pt2D->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "AllRtrue_nu" + std::to_string(nuValue) + ".pdf";
        TCanvas canvas("c", "Multiplot true and four R", 1200, 800);
        canvas.Divide(3, 2); 
        for (int z = 0; z < Rbin; ++z) {
            TMultiGraph *mg = new TMultiGraph();
            double pt2Value = h_nu_z_pt2D->GetZaxis()->GetBinCenter(z + 1);
            canvas.cd(z + 1);
            TGraphErrors *graphSn = new TGraphErrors();
            TGraphErrors *graphCu = new TGraphErrors();
            TGraphErrors *graphC1 = new TGraphErrors();
            TGraphErrors *graphC2 = new TGraphErrors();
            for (int y = 0; y < Rbin; ++y) {
                double zValue = h_nu_z_pt2D->GetYaxis()->GetBinCenter(y + 1);
                double valueSn = ratMatrix[x][y][z];
                double errorSn = errorMatrix[x][y][z];
                double valueCu = ratioCu.getRatMatrix()[x][y][z];
                double errorCu = ratioCu.getErrorMatrix()[x][y][z];
                double valueC1 = ratioC1.getRatMatrix()[x][y][z];
                double errorC1 = ratioC1.getErrorMatrix()[x][y][z];
                double valueC2 = ratioC2.getRatMatrix()[x][y][z];
                double errorC2 = ratioC2.getErrorMatrix()[x][y][z];
                graphSn->SetPoint(y, zValue,   valueSn);
                graphSn->SetPointError(y, 0.0, errorSn); 
                graphCu->SetPoint(y, zValue+0.01,   valueCu);
                graphCu->SetPointError(y+0.01, 0.0, errorCu);
                graphC1->SetPoint(y, zValue+0.02,   valueC1);
                graphC1->SetPointError(y+0.02, 0.0, errorC1);
                graphC2->SetPoint(y, zValue+0.03,   valueC2); //+0.005 to avoid overlap and putting it next its correcponding REC data 
                graphC2->SetPointError(y+0.03, 0.0, errorC2);
            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << pt2Value;
            std::string formattedPt2Value = ss.str();
            graphSn->SetTitle(("R vs z, pt^{2}=" + std::to_string(pt2Value)).c_str());
            graphSn->GetXaxis()->SetTitle("z");
            graphSn->GetYaxis()->SetTitle("R");
            graphSn->SetMarkerStyle(20);
            graphCu->SetMarkerStyle(20);
            graphC1->SetMarkerStyle(20);
            graphC2->SetMarkerStyle(20);
            graphSn->GetYaxis()->SetRangeUser(0.0, 2.0); // Set Y axis range from 0.0 to 2.0
            graphCu->SetMarkerColor(kRed);
            graphSn->SetMarkerColor(kGreen);
            graphC1->SetMarkerColor(kBlue+1);
            graphC2->SetMarkerColor(kBlack);
            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graphSn, "Sn", "lp");
            legend->AddEntry(graphCu, "Cu", "lp");
            legend->AddEntry(graphC1, "C1", "lp");
            legend->AddEntry(graphC2, "C2", "lp");
            TLine *line = new TLine(graphSn->GetXaxis()->GetXmin(), 1.0, graphSn->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line
            mg->Add(graphSn);
            mg->Add(graphCu);
            mg->Add(graphC1);
            mg->Add(graphC2);
            mg->SetTitle(("R vs z, pt2=" + formattedPt2Value).c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("R");
            mg->Draw("APE1");
            legend->Draw("same");
            line->Draw("same");
            TLatex* prelimText = new TLatex();
            prelimText->SetTextSize(0.08);  // larger text size
            prelimText->SetTextAngle(45);
            prelimText->SetTextColorAlpha(kGray + 1, 0.3);  // ray color with transparency
            prelimText->SetNDC();
            prelimText->SetTextAlign(22);  // centered alignment
            prelimText->DrawLatex(0.5, 0.5, "preliminary");
        }
        canvas.SaveAs(pdfFileName.c_str());
    }
}



void Ratio::multiRall(  Ratio& ratioCu, Ratio& ratioC1,   Ratio& ratioC2, Ratio& ratioSnsim , Ratio& ratioCusim, Ratio& ratioC1sim, Ratio& ratioC2sim){
    //fct only for true. arguments only three ratios Cu, C1, C2. Third is included in the fct-class
    for (int x = 0; x < Rbin; ++x) {
        double nuValue = h_nu_z_pt2D->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "AllRtrue_nu" + std::to_string(nuValue) + ".pdf";
        TCanvas canvas("c", "Multiplot true and four R", 1200, 800);
        canvas.Divide(3, 2); 
        for (int z = 0; z < Rbin; ++z) {
            TMultiGraph *mg = new TMultiGraph();
            double pt2Value = h_nu_z_pt2D->GetZaxis()->GetBinCenter(z + 1);
            canvas.cd(z + 1);
            TGraphErrors *graphSn = new TGraphErrors();
            TGraphErrors *graphCu = new TGraphErrors();
            TGraphErrors *graphC1 = new TGraphErrors();
            TGraphErrors *graphC2 = new TGraphErrors();
            TGraphErrors *graphSnsim = new TGraphErrors();
            TGraphErrors *graphCusim = new TGraphErrors();
            TGraphErrors *graphC1sim = new TGraphErrors();
            TGraphErrors *graphC2sim = new TGraphErrors();
            for (int y = 0; y < Rbin; ++y) {
                double zValue = h_nu_z_pt2D->GetYaxis()->GetBinCenter(y + 1);
                double valueSn = ratMatrix[x][y][z];
                double errorSn = errorMatrix[x][y][z];
                double valueCu = ratioCu.getRatMatrix()[x][y][z];
                double errorCu = ratioCu.getErrorMatrix()[x][y][z];
                double valueC1 = ratioC1.getRatMatrix()[x][y][z];
                double errorC1 = ratioC1.getErrorMatrix()[x][y][z];
                double valueC2 = ratioC2.getRatMatrix()[x][y][z];
                double errorC2 = ratioC2.getErrorMatrix()[x][y][z];
                double valueSnsim = ratioSnsim.getRatMatrix()[x][y][z];
                double errorSnsim = ratioSnsim.getErrorMatrix()[x][y][z];
                double valueCusim = ratioCusim.getRatMatrix()[x][y][z];
                double errorCusim = ratioCusim.getErrorMatrix()[x][y][z];
                double valueC1sim = ratioC1sim.getRatMatrix()[x][y][z];
                double errorC1sim = ratioC1sim.getErrorMatrix()[x][y][z];
                double valueC2sim = ratioC2sim.getRatMatrix()[x][y][z];
                double errorC2sim = ratioC2sim.getErrorMatrix()[x][y][z];

                graphSn->SetPoint(y, zValue,   valueSn);
                graphSn->SetPointError(y, 0.0, errorSn); 
                graphCu->SetPoint(y, zValue+0.01,   valueCu);
                graphCu->SetPointError(y+0.01, 0.0, errorCu);
                graphC1->SetPoint(y, zValue+0.02,   valueC1);
                graphC1->SetPointError(y+0.02, 0.0, errorC1);
                graphC2->SetPoint(y, zValue+0.03,   valueC2); //+0.005 to avoid overlap and putting it next its correcponding REC data 
                graphC2->SetPointError(y+0.03, 0.0, errorC2);
                graphSnsim->SetPoint(y, zValue+0.005,   valueSnsim);
                graphSnsim->SetPointError(y+0.005, 0.0, errorSnsim); 
                graphCusim->SetPoint(y, zValue+0.015,   valueCusim);
                graphCusim->SetPointError(y+0.015, 0.0, errorCusim);
                graphC1sim->SetPoint(y, zValue+0.025,   valueC1sim);
                graphC1sim->SetPointError(y+0.025, 0.0, errorC1sim);
                graphC2sim->SetPoint(y, zValue+0.035,   valueC2sim); //+0.005 to avoid overlap and putting it next its correcponding REC data 
                graphC2sim->SetPointError(y+0.035, 0.0, errorC2sim);

            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << pt2Value;
            std::string formattedPt2Value = ss.str();
            std::string title = "R vs z, p_{t}^{2}=" + formattedPt2Value + " GeV^{2}";          /// how to display proper title!!!! here///
            graphSn->SetTitle(title.c_str());
            graphSn->GetXaxis()->SetTitle("z");
            graphSn->GetYaxis()->SetTitle("R");
            graphSn->SetMarkerStyle(20);
            graphSnsim->SetMarkerStyle(20);
            graphCu->SetMarkerStyle(20);
            graphCusim->SetMarkerStyle(20);
            graphC1->SetMarkerStyle(20);
            graphC1sim->SetMarkerStyle(20);
            graphC2->SetMarkerStyle(20);
            graphC2sim->SetMarkerStyle(20);
            graphSn->GetYaxis()->SetRangeUser(0.0, 2.0); // Set Y axis range from 0.0 to 2.0
            graphCu->SetMarkerColor(kRed);
            graphCusim->SetMarkerColor(kOrange+8);            
            graphSn->SetMarkerColor(kGreen);
            graphSnsim->SetMarkerColor(kSpring-7);            
            graphC1->SetMarkerColor(kBlue+1);
            graphC1sim->SetMarkerColor(kCyan-3);            
            graphC2->SetMarkerColor(kBlack);
            graphC2sim->SetMarkerColor(kGray);            
            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graphSn, "Sn", "lp");
            legend->AddEntry(graphCu, "Cu", "lp");
            legend->AddEntry(graphC1, "C1", "lp");
            legend->AddEntry(graphC2, "C2", "lp");
            legend->AddEntry(graphSnsim , "Sn(sim)", "lp");
            legend->AddEntry(graphCusim , "Cu(sim)", "lp");
            legend->AddEntry(graphC1sim , "C1(sim)", "lp");
            legend->AddEntry(graphC2sim , "C2(sim)", "lp");
            
            TLine *line = new TLine(graphSn->GetXaxis()->GetXmin(), 1.0, graphSn->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line
            mg->Add(graphSn);
            mg->Add(graphCu);
            mg->Add(graphC1);
            mg->Add(graphC2);
            mg->Add(graphSnsim);
            mg->Add(graphCusim);
            mg->Add(graphC1sim);
            mg->Add(graphC2sim);

            mg->SetTitle(title.c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("R");
            mg->Draw("APE1");
            legend->Draw("same");
            line->Draw("same");
            TLatex* prelimText = new TLatex();
            prelimText->SetTextSize(0.08);  // larger text size
            prelimText->SetTextAngle(45);
            prelimText->SetTextColorAlpha(kGray + 1, 0.3);  // ray color with transparency
            prelimText->SetNDC();
            prelimText->SetTextAlign(22);  // centered alignment
            prelimText->DrawLatex(0.5, 0.5, "preliminary");
        }
        canvas.SaveAs(pdfFileName.c_str());
    }
}

void Ratio::multiRall2(  Ratio& ratioCu,    Ratio& ratioC, Ratio& ratioSnsim , Ratio& ratioCusim,  Ratio& ratioCsim){
    //fct only for true. arguments only three ratios Cu, CxC. with sims 
    for (int x = 0; x < Rbin; ++x) {
        double nuValue = h_nu_z_pt2D->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "AllRtrueCC_nu" + std::to_string(nuValue) + ".pdf";
        TCanvas canvas("c", "Multiplot trueCC and four R", 1200, 800);
        canvas.Divide(3, 2); 
        for (int z = 0; z < Rbin; ++z) {
            TMultiGraph *mg = new TMultiGraph();
            double pt2Value = h_nu_z_pt2D->GetZaxis()->GetBinCenter(z + 1);
            canvas.cd(z + 1);
            TGraphErrors *graphSn = new TGraphErrors();
            TGraphErrors *graphCu = new TGraphErrors();
            TGraphErrors *graphC = new TGraphErrors();
            TGraphErrors *graphSnsim = new TGraphErrors();
            TGraphErrors *graphCusim = new TGraphErrors();
            TGraphErrors *graphCsim = new TGraphErrors();
            for (int y = 0; y < Rbin; ++y) {
                double zValue = h_nu_z_pt2D->GetYaxis()->GetBinCenter(y + 1);
                double valueSn = ratMatrix[x][y][z];
                double errorSn = errorMatrix[x][y][z];
                double valueCu = ratioCu.getRatMatrix()[x][y][z];
                double errorCu = ratioCu.getErrorMatrix()[x][y][z];
                double valueC = ratioC.getRatMatrix()[x][y][z];
                double errorC = ratioC.getErrorMatrix()[x][y][z];
                double valueSnsim = ratioSnsim.getRatMatrix()[x][y][z];
                double errorSnsim = ratioSnsim.getErrorMatrix()[x][y][z];
                double valueCusim = ratioCusim.getRatMatrix()[x][y][z];
                double errorCusim = ratioCusim.getErrorMatrix()[x][y][z];
                double valueCsim = ratioCsim.getRatMatrix()[x][y][z];
                double errorCsim = ratioCsim.getErrorMatrix()[x][y][z];

                graphSn->SetPoint(y, zValue,   valueSn);
                graphSn->SetPointError(y, 0.0, errorSn); 
                graphCu->SetPoint(y, zValue+0.01,   valueCu);
                graphCu->SetPointError(y+0.01, 0.0, errorCu);
                graphC->SetPoint(y, zValue+0.03,   valueC); //+0.005 to avoid overlap and putting it next its correcponding REC data 
                graphC->SetPointError(y+0.03, 0.0, errorC);
                graphSnsim->SetPoint(y, zValue+0.005,   valueSnsim);
                graphSnsim->SetPointError(y+0.005, 0.0, errorSnsim); 
                graphCusim->SetPoint(y, zValue+0.015,   valueCusim);
                graphCusim->SetPointError(y+0.015, 0.0, errorCusim);
                graphCsim->SetPoint(y, zValue+0.035,   valueCsim); //+0.005 to avoid overlap and putting it next its correcponding REC data 
                graphCsim->SetPointError(y+0.035, 0.0, errorCsim);

            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << pt2Value;
            std::string formattedPt2Value = ss.str();
            std::string title = "R vs z, p_{t}^{2}=" + formattedPt2Value + " GeV^{2}";          /// how to display proper title!!!! here///
            graphSn->SetTitle(title.c_str());
            graphSn->GetXaxis()->SetTitle("z");
            graphSn->GetYaxis()->SetTitle("R");
            graphSn->SetMarkerStyle(20);
            graphSnsim->SetMarkerStyle(20);
            graphCu->SetMarkerStyle(20);
            graphCusim->SetMarkerStyle(20);
            graphC->SetMarkerStyle(20);
            graphCsim->SetMarkerStyle(20);
            graphSn->GetYaxis()->SetRangeUser(0.0, 2.0); // Set Y axis range from 0.0 to 2.0
            graphCu->SetMarkerColor(kRed);
            graphCusim->SetMarkerColor(kOrange+8);            
            graphSn->SetMarkerColor(kGreen);
            graphSnsim->SetMarkerColor(kSpring-7);            
            graphC->SetMarkerColor(kBlack);
            graphCsim->SetMarkerColor(kGray);            
            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graphSn, "Sn", "lp");
            legend->AddEntry(graphCu, "Cu", "lp");
            legend->AddEntry(graphC, "CxC", "lp");
            legend->AddEntry(graphSnsim , "Sn (sim)", "lp");
            legend->AddEntry(graphCusim , "Cu (sim)", "lp");
            legend->AddEntry(graphCsim , "CxC(sim)", "lp");
            
            TLine *line = new TLine(graphSn->GetXaxis()->GetXmin(), 1.0, graphSn->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // Dotted line
            mg->Add(graphSn);
            mg->Add(graphCu);
            mg->Add(graphC);
            mg->Add(graphSnsim);
            mg->Add(graphCusim);
            mg->Add(graphCsim);

            mg->SetTitle(title.c_str());
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("R");
            mg->Draw("APE1");
            legend->Draw("same");
            line->Draw("same");
            TLatex* prelimText = new TLatex();
            prelimText->SetTextSize(0.08);  // larger text size
            prelimText->SetTextAngle(45);
            prelimText->SetTextColorAlpha(kGray + 1, 0.3);  // ray color with transparency
            prelimText->SetNDC();
            prelimText->SetTextAlign(22);  // centered alignment
            prelimText->DrawLatex(0.5, 0.5, "preliminary");
        }
        canvas.SaveAs(pdfFileName.c_str());
    }
}

void Ratio::Rtargetsimcomp( Ratio& ratiosim){
    for (int x = 0; x < Rbin; ++x) {

        double nuValue = h_nu_z_pt2D->GetXaxis()->GetBinCenter(x + 1);
        std::string pdfFileName = "singletargetcomp_" + targetName + "_nu" + std::to_string(nuValue) + ".pdf";
        TCanvas canvas("c", "Multiplot double R", 1200, 800);
        canvas.Divide(3, 2); 
        for (int z = 0; z < Rbin; ++z) {
            TMultiGraph *mg = new TMultiGraph();

            double pt2Value = h_nu_z_pt2D->GetZaxis()->GetBinCenter(z + 1);
            canvas.cd(z + 1);
            TGraphErrors *graph = new TGraphErrors();
            TGraphErrors *graphsim = new TGraphErrors();
            
            for (int y = 0; y < Rbin; ++y) {
                double zValue = h_nu_z_pt2D->GetYaxis()->GetBinCenter(y + 1);
                double value = ratMatrix[x][y][z];
                double error = errorMatrix[x][y][z];
                double valuesim = ratiosim.getRatMatrix()[x][y][z];
                double errorsim = ratiosim.getErrorMatrix()[x][y][z];

                graph->SetPoint(y, zValue, value);
                graph->SetPointError(y, 0.0, error); 
                graphsim->SetPoint(y, zValue+0.01, valuesim);
                graphsim->SetPointError(y+0.01, 0.0, errorsim);
                

            }

            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << pt2Value;
            std::string formattedPt2Value = ss.str();
            //std::cout << "formattedPt2Value = " << formattedPt2Value << std::endl;
            std::string title = "R vs z for " + targetName + ", p_{t}^{2}=" + formattedPt2Value + " GeV^{2}";
            //std::string title = "R vs z, p_{t}^{2}=" + formattedPt2Value + " GeV^{2}";          /// how to display proper title!!!! here///

            graph->SetTitle(title.c_str());
            graph->GetXaxis()->SetTitle("z");
            graph->GetYaxis()->SetTitle("R");
            graph->SetMarkerStyle(20);
            graphsim->SetMarkerStyle(20);
            graph->GetYaxis()->SetRangeUser(0.0, 2.0); // this seems useless
            //graph->Draw("AP");
            graphsim->SetMarkerColor(kRed);

            graph->SetMarkerColor(kBlue);
            

            TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(graph, "data", "lp");
            legend->AddEntry(graphsim, "simulated", "lp");

            TLine *line = new TLine(graph->GetXaxis()->GetXmin()+0.02, 1.0, graph->GetXaxis()->GetXmax(), 1.0);
            line->SetLineStyle(2); // dotted line
            line->SetLineColor(kGray + 1); //grayer for aethetic 
            mg->Add(graph);
            mg->Add(graphsim);
            mg->SetTitle((title).c_str() );
            mg->GetXaxis()->SetTitle("z");
            mg->GetYaxis()->SetTitle("R");
            mg->GetYaxis()->SetRangeUser(-0.5, 1.5); // 

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

        canvas.SaveAs(pdfFileName.c_str());
    }
}


void Ratio::ValidateHistograms() {
    // Validate histograms for both deuterium and target nucleus
    std::cout << "RATIO VALIDATION " <<std::endl;
    std::cout << "Histogram h_nuD has " << h_nuD->GetEntries() << " entries." << std::endl;
    std::cout << "Histogram h_nuA has " << h_nuA->GetEntries() << " entries." << std::endl;
    std::cout << "Histogram h_nu_z_pt2D has " << h_nu_z_pt2D->GetEntries() << " entries." << std::endl;
    std::cout << "Histogram h_nu_z_pt2A has " << h_nu_z_pt2A->GetEntries() << " entries." << std::endl;

    if (h_nuD->GetEntries() == 0 || h_nuA->GetEntries() == 0) {
        std::cerr << "Warning: Histograms h_nuD or h_nuA are empty!" << std::endl;
    }
    if (h_nu_z_pt2D->GetEntries() == 0) {
        std::cerr << "Warning: Histogram h_nu_z_pt2D is empty!" << std::endl;
    }
    if (h_nu_z_pt2A->GetEntries() == 0) {
        std::cerr << "Warning: Histogram h_nu_z_pt2A is empty!" << std::endl;
    }
}

//void Ratio::LogBinContent() {
//    // Log bin contents for both h_nu_z_pt2D and h_nu_z_pt2A
//    std::cout << "Logging bin content for h_nu_z_pt2D:" << std::endl;
//    for (int x = 1; x <= h_nu_z_pt2D->GetNbinsX(); ++x) {
//        for (int y = 1; y <= h_nu_z_pt2D->GetNbinsY(); ++y) {
//            for (int z = 1; z <= h_nu_z_pt2D->GetNbinsZ(); ++z) {
//                double content = h_nu_z_pt2D->GetBinContent(x, y, z);
//                if (content > 0) {
//                    std::cout << "Bin (D) (" << x << ", " << y << ", " << z << ") = " << content << std::endl;
//                }
//            }
//        }
//    }
//
//    std::cout << "Logging bin content for h_nu_z_pt2A:" << std::endl;
//    for (int x = 1; x <= h_nu_z_pt2A->GetNbinsX(); ++x) {
//        for (int y = 1; y <= h_nu_z_pt2A->GetNbinsY(); ++y) {
//            for (int z = 1; z <= h_nu_z_pt2A->GetNbinsZ(); ++z) {
//                double content = h_nu_z_pt2A->GetBinContent(x, y, z);
//                if (content > 0) {
//                    std::cout << "Bin (A) (" << x << ", " << y << ", " << z << ") = " << content << std::endl;
//                }
//            }
//        }
//    }
//}

void Ratio::LogBinContent() {
    // generalized histo selectioner
    auto logHistogram = [](TH3F* hist, const std::string& name) {
        //careful w/ auto  ???
        int numBinsX = hist->GetNbinsX();
        int numBinsY = hist->GetNbinsY();
        int numBinsZ = hist->GetNbinsZ();

        std::cout << "RATIO BIN CONTENT FOR " << name << " (3D):" << std::endl;

        for (int x = 1; x <= numBinsX; ++x) { // Loop over 'nu' bins
            double nuValue = hist->GetXaxis()->GetBinCenter(x);
            std::cout << "nu = " << nuValue << " (Layer " << x << "):" << std::endl;

            // Print column headers (z values)
            std::cout << std::setw(10) << "pt2 \\ z";
            for (int y = 1; y <= numBinsY; ++y) {
                double zValue = hist->GetYaxis()->GetBinCenter(y);
                std::cout << std::setw(10) << zValue;
            }
            std::cout << std::endl;

            // Print matrix content for fixed 'nu' layer
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

    // Log both histograms
    logHistogram(h_nu_z_pt2D, "h_nu_z_pt2D");
    logHistogram(h_nu_z_pt2A, "h_nu_z_pt2A");
}
