#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <TLorentzVector.h>

namespace Constants {
    // Masses
    const double MASS_ELECTRON = 0.000511;  // GeV/c^2
    const double MASS_NEUTRON = 0.939565;   // GeV/c^2
    const double MASS_PROTON = 0.938272;    // GeV/c^2
    const double MASS_PION = 0.134976;      // GeV/c^2
    static const int ELECTRON_PID = 11;
    static const int POSITRON_PID = -11;


    //to call these values use in declaration e.g.:     Constants::MASS_ELECTRON;

    // Incident lepton
    const TLorentzVector elBeam(0.0, 0.0, 11.0, 11.0); // Define your values here

    // Other constants...
}

#endif // CONSTANTS_H
