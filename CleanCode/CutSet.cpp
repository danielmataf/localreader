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

void CutSet::SetCutGen4Rat(){

    SetCutQ(Constants::RcutminQ,Constants::RcutmaxQ );
    SetCutY(Constants::RcutminY,Constants::RcutmaxY );
    SetCutW(Constants::RcutminW,Constants::RcutmaxW );
    SetCutZ(Constants::RcutminZ,Constants::RcutmaxZ );
    SetCutPt2(Constants::RcutminPt2,Constants::RcutmaxPt2 );    
    //cut in vertex should be fixed here (?)
}

//both functions could be used? propagate to main if below is used.
void CutSet::SetCutGentest(const Event& event){
    int target = event.GetTargetType();
    SetCutQ(Constants::RcutminQ,Constants::RcutmaxQ );
    SetCutY(Constants::RcutminY,Constants::RcutmaxY );
    SetCutW(Constants::RcutminW,Constants::RcutmaxW );
    SetCutZ(Constants::RcutminZ,Constants::RcutmaxZ );
    SetCutPt2(Constants::RcutminPt2,Constants::RcutmaxPt2 );    
    // switch statement could be used, but almost same thing as if statement
    if (target == 0) {
        SetCutVz(Constants::RcutminVzLD2,Constants::RcutmaxVzLD2 );
    }
    if (target == 1) {
        SetCutVz(Constants::RcutminVzSn,Constants::RcutmaxVzSn );
    }
    if (target == 2) {
        SetCutVz(Constants::RcutminVzCu,Constants::RcutmaxVzCu );
    }
    if (target == 3) {
        SetCutVz(Constants::RcutminVzC,Constants::RcutmaxVzC );
    }
}


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
    if (Vz >= cutVzMin && Vz <= cutVzMax ){
        if (Q2 >= cutQMin && Q2 <= cutQMax ){
            if (y >= cutYMin && y <= cutYMax){
                if (v >= cutVMin && v <= cutVMax ){
                    if (w >= cutWMin && w <= cutWMax ){
                        return true;
                    }
                }
            }
        }
    }    
    return false;
}

bool CutSet::PassCutsHadrons( const Particle& hadron)  {
    double z = hadron.Getz();
    double pt2 = hadron.Getpt2();
    double phih = hadron.Getphih();
    double h_pid= hadron.GetPID();
    if (z >= cutZMin && z <= cutZMax) {
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

