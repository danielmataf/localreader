#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <TLorentzVector.h>

namespace Constants {
    // MASSES
    const double MASS_ELECTRON = 0.000511;  // GeV/c^2
    const double MASS_NEUTRON = 0.939565;   // GeV/c^2
    const double MASS_PROTON = 0.938272;    // GeV/c^2
    const double MASS_PION = 0.134976;      // GeV/c^2
    static const double MASS_PION_PLUS  = 0.13957;
    static const double MASS_PION_MINUS = 0.13957;
    static const double MASS_PION_ZERO  = 0.13497;
    static const double default_mass  = 0.001;

    // COUNTERS

    // PARTICLE PDG ID 
    static const int ELECTRON_PID = 11;
    static const int POSITRON_PID = -11;
    static const int PION_PLUS_PID = 211;
    static const int PION_MINUS_PID = -211;
    static const int PION_ZERO_PID = 111;
    // VALUES 
    const double PI = 3.14159265;

    //to call these values use in declaration e.g.:     Constants::MASS_ELECTRON;

    // Vectors     
    const TLorentzVector elBeam(0.0, 0.0, 11.0, 11.0); // incident lepton


    //Cut values for Ratio 
    //SetCutQ(1.5,10);
    //SetCutY(0.25, 0.85);
    //SetCutW(6,30);
    //SetCutZ(0.3,0.7);
    //SetCutPt2(0.0,2.0);
    static const double RcutminQ = 1.5 ;
    static const double RcutminY = 0.25 ;
    static const double RcutminW = 6.0 ;
    static const double RcutminZ = 0.3 ;
    static const double RcutminPt2 = 0.0 ;
    static const double RcutmaxQ = 10.0 ;
    static const double RcutmaxY = 0.85 ;
    static const double RcutmaxW = 30.0 ;
    static const double RcutmaxZ = 0.7 ;
    static const double RcutmaxPt2 = 2.0 ;   
    static const double Rcutminnu = 4.0 ;
    static const double Rcutmaxnu = 9.0 ;


    //cuts for Target specific vertex 
    static const double RcutminVzLD2 = -7.5 ;   
    static const double RcutmaxVzLD2 = -2.5 ;   
    static const double RcutminVzSn = -3.5 ;
    static const double RcutmaxVzSn = -1.5 ;
        //not used yet
    static const double RcutminVzCu = -8.5 ;    //arbitrary +-.5
    static const double RcutmaxVzCu = -6.5 ;    //arbitrary +-.5
    static const double RcutminVzC = -7.5 ;     //default (for LD2)
    static const double RcutmaxVzC = -2.5 ;     //default (for LD2)
    

    
}

#endif 