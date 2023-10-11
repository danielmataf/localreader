#include "CutSet.h"
#include "Event.h"






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

    //OTHER CUTS 
    // (theta_el*180/PI>6)	// /!\ implementing CUT on theta coordinate for electrons!!!!
    // (mmass>1.2)		    //To remove the proton.
    // (v_scapip->E()<3 && v_scapip->E()>1)
}

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


bool CutSet::PassCutsElectrons(const Event& event)  {
    // recover kinematic variables from the event
    double Q2 = event.GetQ2();
    double y  = event.Gety();
    double v  = event.Getnu();
    double w  = event.GetW2();
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

    if (Q2 >= cutQMin && Q2 <= cutQMax ){
        
        if (y >= cutYMin && y <= cutYMax){
            if (v >= cutVMin && v <= cutVMax ){
                if (w >= cutWMin && w <= cutWMax ){
                      return true;
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

bool CutSet::PassCutsHadrons(const Particle& hadron)  {
    double z = hadron.Getz();
    double pt2 = hadron.Getpt2();
    double phih = hadron.Getphih();
    
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