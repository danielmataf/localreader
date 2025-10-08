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

    cutnuMin = 0.0;
    cutnuMax = 10.0;


    //PID defaults
    hadPID = Constants::PION_PLUS_PID;

    //Chi2 defaults
    elchi2min= -5.0;
    elchi2max= 5.0;

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
    SetCutnu(Constants::Rcutminnu,Constants::Rcutmaxnu );
    SetCutW(Constants::RcutminW,Constants::RcutmaxW );
    SetCutZ(Constants::RcutminZ,Constants::RcutmaxZ );
    SetCutPt2(Constants::RcutminPt2,Constants::RcutmaxPt2 );  
    SetCutnu(Constants::RcutnuMin,Constants::RcutnuMax); //avoid acceptance corrections defined with R to keep coherence  
    //cut in vertex should be fixed here (?)
}


//I think the function below is trash ( 27 nov 2024)
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

//PID
void CutSet::SetCutHadPID( int hadronPID){
    hadPID = hadronPID;
}

void CutSet::SetCutelPID(int elPID){
    elPID = Constants::ELECTRON_PID;
}
//


//CHI2
void CutSet::SetCutChi2el(double minChi2, double maxChi2){
    elchi2min = minChi2;
    elchi2max = maxChi2;
}


//u,v,w cuts
void CutSet::SetCutlu(double minlu, double maxlu){
    cutluMin = minlu;
    cutluMax = maxlu;
}
void CutSet::SetCutlv(double minlv, double maxlv){
    cutlvMin = minlv;
    cutlvMax = maxlv;
}
void CutSet::SetCutlw(double minlw, double maxlw){
    cutlwMin = minlw;
    cutlwMax = maxlw;
}

//Nphe cuts
void CutSet::SetCutNphe15(double minNphe_15){
    cutNphe15 = minNphe_15;
}
void CutSet::SetCutNphe16(double minNphe_16){
    cutNphe16 = minNphe_16;
}

void CutSet::SetCutnu(double minnu, double maxnu){
    cutnuMin = minnu;
    cutnuMax = maxnu;
}



//check if an event passes the cuts
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

bool CutSet::PassCutOnlyVz(const Event& event){     //assimilate to cutset class with event.GetTarget to generalize TBD!!
    //This is only for LD2 ????
    //need to generalalize ASAP. Check if this is being used/propagated
    double Vz = event.GetVz();
    if (Vz >= Constants::RcutminVzLD2 && Vz <= Constants::RcutmaxVzLD2 ){    
        return true;
    }
    return false;
}

bool CutSet::PassCutxb(const Event& event, double minxb)  {
    double xb = event.Getxb();
    if (xb > minxb ){
        return true;
    }
    return false;
}

bool CutSet::PassCutxbMC(const Event& event)  {
    double xb = event.GetxbMC();
    if (xb > 0.8 ){
        return true;
    }
    return false;
}

bool CutSet::PassCutsElectrons(const Event& event)  {
    // recover kinematic variables from the event
    // Pass cut electrons also filters Vz initially!! (cool)
    double Q2 = event.GetQ2();
    double y  = event.Gety();
    double v  = event.Getnu();
    double w  = event.GetW2();
    double Vz = event.GetVz();
    double Vx = event.GetVx();
    double Vy = event.GetVy();
    if (Vz >= cutVzMin && Vz <= cutVzMax ){
        //std::cout << "Vz passed" << std::endl;
        if (Q2 >= Constants::RcutminQ && Q2 <= Constants::RcutmaxQ ){
            if (y >= Constants::RcutminY && y <= Constants::RcutmaxY){
                if (v >= Constants::Rcutminnu && v <= Constants::RcutnuMax ){             //not specified... is broad enough not to be considered
                    if (w >= Constants::RcutminW && w <= Constants::RcutmaxW ){
                        //std::cout << "nu value " << event.Getnu()<< std::endl;
                        //std::cout << "y max " << Constants::RcutmaxY<< std::endl;
                        return true;
                    }
                }
            }
        }
    }    
    return false;
}

bool CutSet::PassCutsElectronsMC(const Event& event)  {
    // recover kinematic variables from the event
    // Pass cut electrons also filters Vz initially!! (cool)
    double Q2 = event.GetQ2MC();
    double y  = event.GetyMC();
    double v  = event.GetnuMC();
    double w  = event.GetW2MC();
    double Vz = event.GetVzMC();
    if (Vz >= cutVzMin && Vz <= cutVzMax ){
        //std::cout << "Vz passed" << std::endl;
        if (Q2 >= Constants::RcutminQ && Q2 <= Constants::RcutmaxQ ){
            if (y >= Constants::RcutminY && y <= Constants::RcutmaxY){
                if (v >= Constants::Rcutminnu && v <= Constants::RcutnuMax ){             //not specified... is broad enough not to be considered
                    if (w >= Constants::RcutminW && w <= Constants::RcutmaxW ){
                        //std::cout << "nu value " << event.Getnu()<< std::endl;
                        //std::cout << "y max " << Constants::RcutmaxY<< std::endl;
                        return true;
                    }
                }
            }
        }
    }    
    return false;
}


bool CutSet::PassCutVzselection(const Event& event){
    double Vz = event.GetVz();
    int targetType = event.GetTargetType();

    if (Vz >= cutVzMin && Vz <= cutVzMax ){
        return true;
    }
    return false;
}


bool CutSet::PassCutsHadrons( const Particle& hadron)  {
    double z = hadron.Getz();
    double pt2 = hadron.Getpt2();
    double phih = hadron.Getphih();
    double h_pid= hadron.GetPID();
    double hadVz = hadron.GetParticleVertexZ();
    //PID is not being included !!!
    if (z >= cutZMin && z <= cutZMax) {
        if (pt2 >= Constants::RcutminPt2 && pt2 <= Constants::RcutmaxPt2) {
            //std::cout << "vz value =" << hadVz <<  std::endl;
            return true;
        }
    }
    
    return false;
}

bool CutSet::PassCutsVariables(const Event& event){    
    //All cuts in variables, vz included 
    if (PassCutsElectrons(event) ){
        for (const Particle& hadron : event.GetHadrons()) {
            if (PassCutsHadrons(hadron)){
                return true;
            }
        }
    }
    return false;

}

bool CutSet::PassCutsCalo(const Event& event)  {
    // recover calo data from the event
    double lu = event.electron.Getlu();
    double lv = event.electron.Getlv();
    double lw = event.electron.Getlw();
    double epcal = event.electron.GetEpcal();
    double eecalin = event.electron.GetEcalin();
    double ecalout = event.electron.GetEcalout();

    if (lu >= Constants::cal_lumin  ){
        if (lv >= Constants::cutlv_min ){
            if (lw >= Constants::cutlw_min){
                if (epcal >= Constants::cutEpcal_min  ){
                    return true;
                }
            }
        }
    }
    return false;
}

bool CutSet::PassCutsCherenkov(const Event& event ) {
    // recover cherenkov data from the event
    double Nphe15 = event.electron.Getnphe15();
    double Nphe16 = event.electron.Getnphe16();
    if (Nphe15 >= Constants::Nphe15min   ){
        //if (Nphe16 >= Constants::Nphe16min   ){
            return true;
        //}
    }
    return false;
}

bool CutSet::PassCutsDetectors(const Event& event)  {
    if (PassCutsCalo(event) && PassCutsCherenkov(event) ){
        return true;
    }
    return false;
}


bool CutSet::PassCutsAll(const Event& event)  {
    if ( PassCutsDetectors(event) ){
        if ( PassCutsVariables(event)){
            return true;
        }        
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

bool CutSet::PassCutTarget2Vz(const Event& event)  {
    //this function can only be used when calling a carbon on carbon ratio. 
    //maybe we need to run the whole code just for the C on C ratio function. Annoying but necessary
    //maybe keep the 
    int targetType = event.GetTargetType();
    double Vz = event.GetVz();
    if (targetType == 3) {
        return true;
    }
    return false;
}