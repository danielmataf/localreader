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

}

#endif 