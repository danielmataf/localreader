#ifndef CUTSET_H
#define CUTSET_H

#include "Event.h"  

class CutSet {
public:
    CutSet(); //initialize the cut values

    void SetCutQ(double minQ, double maxQ);
    void SetCutY(double minY, double maxY);
    void SetCutV(double minV, double maxV);
    void SetCutW(double minW, double maxW);
    void SetCutZ(double minZ, double maxZ);

    
    //check if an event passes the cuts
    bool PassCuts(const Event& event) ;
    
    //bool PassCut(double Q2, double y, double v, double w, double z) ;
    //we 
private:
    double cutQMin, cutQMax;
    double cutYMin, cutYMax;
    double cutVMin, cutVMax;
    double cutWMin, cutWMax;
    double cutZMin, cutZMax;
};

#endif 