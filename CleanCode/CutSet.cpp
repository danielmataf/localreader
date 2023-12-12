#include "CutSet.h"
#include "Event.h"
#include <TFile.h>
#include <TCanvas.h>
#include <TPDF.h>
#include <TLine.h>
#include "constants.h"






CutSet::CutSet() {
    // initializing cut values here with exagerated values 
    cutQMin = 0.0;
    cutQMax = 10.0;
    // (Q2>1.5)
    cutYMin = 0.0;
    cutYMax = 10.0;
    // (y<0.85 && y>0.25)
    cutVMin = 0.0;
    cutVMax = 10.0;
    cutWMin = 0.0;
    cutWMax = 10.0;
    //  (W>2)       //careful with squares 
    cutZMin = 0.0;
    cutZMax = 10.0;
    // (z>0.3 && z<0.7)
    cutVzMin = -30.0;
    cutVzMax = 30.0 ;
    cutPt2Min= 0.0;
    cutPt2Max=3.0 ;

    hadPID = Constants::PION_PLUS_PID;


    // target Sn (-3.5, -1.5)
    VzminSn= -3.5;       
    VzmaxSn= -1.5;     
    // target LD2 (-7.5,-2.5)
    VzminLD2= -7.5;    
    VzmaxLD2= -2.5;         
    // target Cu (-8.5, -6.5)   +/-1 position
    VzminCu= -8.5;     
    VzmaxCu= -6.5;     
    //target CxC (?)
    VzminCxC= -10.0;    
    VzmaxCxC= 0.0;    

    //OTHER CUTS 
    // (theta_el*180/PI>6)	// /!\ implementing CUT on theta coordinate for electrons!!!!
    // (mmass>1.2)		    //To remove the proton.
    // (v_scapip->E()<3 && v_scapip->E()>1)
}
CutSet::~CutSet() {
}




//Add cut vertex TBD 
//Add cut pour differentes set up de cible (C vs Sn et/ou Cu)
//Simulation must generate events from a target in the correct position 




void CutSet::SetCutQ(double minQ, double maxQ) {
    //this function can be called in main.cpp but can be used in an eventual eventprocessor.cpp 
    //same for other SetCutVariable function
    cutQMin = minQ;
    cutQMax = maxQ;
}

void CutSet::SetCutY(double minY, double maxY) {
    cutYMin = minY;
    cutYMax = maxY;
}

void CutSet::SetCutV(double minV, double maxV) {
    cutVMin = minV;
    cutVMax = maxV;
}

void CutSet::SetCutW(double minW, double maxW) {
    cutWMin = minW;
    cutWMax = maxW;
}

void CutSet::SetCutZ(double minZ, double maxZ) {
    cutZMin = minZ;
    cutZMax = maxZ;
}
void CutSet::SetCutPt2(double minPt2, double maxPt2){
    cutPt2Min = minPt2;
    cutPt2Max = maxPt2;


}


void CutSet::SetCutVx(double vxmin, double vxmax){
    cutVxMin = vxmin;
    cutVxMax = vxmax;
    
}
void CutSet::SetCutVy(double vymin, double vymax){
    cutVyMin = vymin;
    cutVyMax = vymax;
    
}

void CutSet::SetCutVz(double vzmin, double vzmax){
    cutVzMin = vzmin;
    cutVzMax = vzmax;
    
}


//


void CutSet::SetCutHadPID( int hadronPID){
    hadPID = hadronPID;
}

bool CutSet::PassCutSnTarget(const Event& event){
    double Vz = event.GetVz();
    if (Vz >= VzminSn && Vz <= VzmaxSn ){
        return true;
    }
    return false;
}

bool CutSet::PassCutLD2Target(const Event& event){
    double Vz = event.GetVz();
    if (Vz >= VzminLD2 && Vz <= VzmaxLD2 ){
        return true;
    }
    return false;
}

bool CutSet::PassCutCuTarget(const Event& event){
    double Vz = event.GetVz();
    if (Vz >= VzminCu && Vz <= VzmaxCu ){    
        return true;
    }
    return false;
}


bool CutSet::PassCutsElectrons(const Event& event)  {
    // recover kinematic variables from the event
    double Q2 = event.GetQ2();
    double y  = event.Gety();
    double v  = event.Getnu();
    double w  = event.GetW2();
    double Vz = event.GetVz();
    double Vx = event.GetVx();
    double Vy = event.GetVy();

    //for (const Particle& hadron : event.GetHadrons()) {
    //    double z = hadron.Getz();
    //    std::cout << " z value ok ? = "<< z << std::endl;
    //    
    //} 
    //double z  = event.GetZ();             //hadron variable TBD!!
    // Check if the event's kinematic variables pass the cuts
    //std::cout<< Q2<<" .........Q2......"<<std::endl;
    //std::cout<< y<<" .........y......"<<std::endl;
    //std::cout<< v<<" .........v......"<<std::endl;
    //std::cout<< w<<" .........W......"<<std::endl;
    // h_Vzpre->Fill(Vz);
    if (Vz >= cutVzMin && Vz <= cutVzMax ){
        //h_Vxpre->Fill(Vx);
        //h_Vypre->Fill(Vy);
        //h_Vzpos->Fill(Vz);
        //h_Q2pre->Fill(Q2);
        if (Q2 >= cutQMin && Q2 <= cutQMax ){
            //h_Q2pos->Fill(Q2); 
            //h_xbpre->Fill(event.Getxb());
            //h_Qvx->Fill(event.Getxb(), Q2);
            //h_ypre->Fill(y);
            if (y >= cutYMin && y <= cutYMax){
                //h_nupre->Fill(v);
                if (v >= cutVMin && v <= cutVMax ){
                    //h_W2pre->Fill(w);
                    if (w >= cutWMin && w <= cutWMax ){
                        //h_Q4R->Fill(Q2);
                        //h_nu4R->Fill(v);

                        return true;
                    }
                }
            }
        }
    }    
    return false;
}
    // checks if the event's kinematic variables pass the cuts
//bool CutSet::PassCut(double Q2, double y, double v, double w, double z) const {
//  return true ;    
//}

bool CutSet::PassCutsHadrons( const Particle& hadron)  {
    double z = hadron.Getz();
    double pt2 = hadron.Getpt2();
    double phih = hadron.Getphih();
    double h_pid= hadron.GetPID();
    //h_zpre->Fill(z);
    if (z >= cutZMin && z <= cutZMax) {
        //h_pt2pre->Fill(pt2);
        //h_pt4R->Fill(pt2);
        //h_z4R->Fill(hadron.Getz());
        //h_phihpre->Fill(phih);
        return true;
    }
    
    return false;
}

//do another  fct for hadron cutts (-count hadrons, call function passcuts de l'electron )
//returns  
  
//create funct to recover hadrons in event.cpp 
//fct how many hadrons in an evt 
//create passcutshadron w/ argument is particle hadron  

bool CutSet::PassCuts4R(const Event& event, const Particle& hadron)  {

    double pt2 = hadron.Getpt2();
    double v  = event.Getnu();
    if (pt2 >= 0 && pt2 <= 2) {
        if (v >= 4 && v <=9 ) {
            return true;
        }
    }
    return false;
}















/*
//====Writing and Drawing histos saved below=======//

void CutSet::Chop(const std::string filename) {
    //this function recreates a new rootfile everytime is called 
    //useful to have different rootfiles if different cuts were implemented
    //argument: title of the wanted output file 
    // save histograms to the specified ROOT file
    TFile file = TFile(filename.c_str(), "RECREATE");
    h_Q2pre->Write();
    h_Q2pos->Write();
    h_xbpre->Write();
    h_ypre->Write();
    h_nupre->Write();
    h_W2pre->Write();
    h_zpre->Write();
    h_pt2pre->Write();
    h_phihpre->Write();
    h_Vzpre->Write();
    h_Vxpre->Write();
    h_Vypre->Write();
    h_Vzpos->Write();
    h_Qvx->Write();
    h_Q4R->Write();
    h_nu4R->Write();
    h_z4R->Write();
    h_pt4R->Write();
    //file.Write();
    // Write other histograms here
    file.Close();
}







void CutSet::DrawChop(const std::string filename) {
    TCanvas ChopC("Chop canvas", "Chop Histograms");
    ChopC.Divide(3, 3);
    ChopC.cd(1);
    h_Q2pre->Draw("hist");
    h_Q2pre->SetTitle("Q2 Distribution");
    ChopC.Update();
    TLine *cutQ=new TLine(cutQMin,ChopC.cd(1)->GetUymin(),cutQMin,ChopC.cd(1)->GetUymax());
    cutQ->SetLineWidth(2);
    cutQ->SetLineStyle(kDashed);
    cutQ->Draw("same");


    ChopC.cd(2);
    h_xbpre->Draw("hist");
    h_xbpre->SetTitle("xb Distribution");
    ChopC.cd(3);
    h_ypre->Draw("hist");
    h_ypre->SetTitle("y Distribution");
    
    ChopC.Update();
    TLine *cutY=new TLine(cutYMax,ChopC.cd(3)->GetUymin(),cutYMax,ChopC.cd(3)->GetUymax());
    cutY->SetLineWidth(2);
    cutY->SetLineStyle(kDashed);
    cutY->Draw("");
    
    ChopC.cd(4);
    h_nupre->Draw("hist");
    h_nupre->SetTitle("nu Distribution");

    ChopC.cd(5);
    h_W2pre->Draw("hist");
    h_W2pre->SetTitle("W2 Distribution");
    
    ChopC.Update();
    TLine *cutW=new TLine(cutWMax,ChopC.cd(5)->GetUymin(),cutWMax,ChopC.cd(5)->GetUymax());
    cutW->SetLineWidth(2);
    cutW->SetLineStyle(kDashed);
    cutW->Draw("");

    ChopC.cd(6);
    h_zpre->Draw("hist");
    h_zpre->SetTitle("z Distribution");
    ChopC.Update();
    TLine *cutzlow=new TLine(cutZMin,ChopC.cd(6)->GetUymin(),cutZMin,ChopC.cd(6)->GetUymax());
    cutzlow->SetLineWidth(2);
    cutzlow->SetLineStyle(kDashed);
    cutzlow->Draw("");
    TLine *cutzhigh=new TLine(cutZMax,ChopC.cd(6)->GetUymin(),cutZMax,ChopC.cd(6)->GetUymax());
    cutzhigh->SetLineWidth(2);
    cutzhigh->SetLineStyle(kDashed);
    cutzhigh->Draw("");

    ChopC.cd(7);
    h_pt2pre->Draw("hist");
    h_pt2pre->SetTitle("pt2 Distribution");
    ChopC.cd(8);
    h_phihpre->Draw("hist");
    h_phihpre->SetTitle("phih Distribution");


    ChopC.Print((filename + ".pdf").c_str());
}

*/