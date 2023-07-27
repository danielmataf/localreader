#ifndef VARIABLES_H
#define VARIABLES_H
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TH1F.h"




struct MyVariables {
    double theta_e ;
    double Q2 ;
    double W;
    double y;
    double gamnu;
    double x_b;
    double t;
    double z;
    double P_t;
    double phih;
    double Epho;
    double theta_el;
    double theta_pip;
    double phi_pip;
    double phi_el;
    double angle_norm;
    double sign;
    double rawphih;
    double rawphih2;
    double p_el;
    double p_pho;
    double P_trX;   //variables 4 transverse momentum of pion and z, then transverse momentum for hadron X    //---[De]
    double pipX, pipY, pipZ;
    double scalX, scalY;
    double px_scal, py_scal, pz_scal, px_scapi, py_scapi, pz_scapi, px_scah, py_scah, pz_scah;
    double mmass;
    double cone_ang;

    
};

#endif 