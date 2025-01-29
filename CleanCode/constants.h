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
    const TLorentzVector elBeam(0.0, 0.0, 10.5, 10.5); // incident lepton


    //Cut values for Ratio 
    //SetCutQ(1.5,10);
    //SetCutY(0.25, 0.85);
    //SetCutW(6,30);
    //SetCutZ(0.3,0.7);
    //SetCutPt2(0.0,2.0);

            //==Binning & limits for Monitoring==//
        //=monitor histograms (bin & limits)
    //physics values
    static const int bin_default = 100;
    static const int phihbin_default = 10;
    static const double Qmin_default = 0.0;
    static const double Qmax_default = 10.0;
    static const double xmin_default = 0.0;
    static const double xmax_default = 0.4;
    static const double numin_default = 0.0;
    static const double numax_default = 10.0;
    static const double Ymin_default = 0.0;
    static const double Ymax_default = 1.0;
    static const double Wmin_default = 0.0;
    static const double Wmax_default = 10.0;
    static const double zmin_default = 0.0;
    static const double zmax_default = 1.0;
    static const double pt2min_default = 0.0;
    static const double pt2max_default = 3.0;
    static const double phihmin_default = 0.0;
    static const double phihmax_default = 360.0;

    //coord values
    static const double phielmin_default = -180.0;
    static const double phielmax_default = 180.0;
    static const double thetaelmin_default = 0.0;
    static const double thetaelmax_default = 50.0;
    static const double phipimin_default = -180.0;
    static const double phipimax_default = 180.0;
    static const double thetapimin_default = 0.0;
    static const double thetapimax_default = 50.0;


    //vertex values 
    static const double Vxmin_default = -10.0;
    static const double Vxmax_default = 10.0;
    static const double Vymin_default = -10.0;
    static const double Vymax_default = 10.0;
    static const double Vzmin_default = -40.0;
    static const double Vzmax_default = 40.0;


    //Calo Values 
    static const double cal_xmin = -400.0;
    static const double cal_xmax = 400.0;      //values assigned from RGA
    static const double cal_ymin = -400.0;     //values assigned from RGA
    static const double cal_ymax = 400.0;
    static const double cal_zmin = -400.0;     //not defined yet
    static const double cal_zmax = 400.0;      //not defined yet
    static const double cal_lumin = 0.0;       
    static const double cal_lumax = 50.0;
    static const double cal_lvmin = 0.0;
    static const double cal_lvmax = 50.0;
    static const double cal_lwmin = 0.0;
    static const double cal_lwmax = 50.0;   
    static const double cal_Epmin = 0.06;
    static const double cal_Epmax = 0.40;  

    //Cherenkov Values
    static const double Nphe15min = 2.0;
    static const double Nphe16min = 2.0;

    //Binning & limits for Ratios

    static const double RcutminQ = 1.0 ;
    static const double RcutminY = 0.25 ;
    static const double RcutminW = 4.0 ;
    static const double RcutminZ = 0.3 ;
    static const double RcutminPt2 = 0.0 ;
    static const double RcutmaxQ = 5.0 ;
    static const double RcutmaxY = 0.85 ;
    static const double RcutmaxW = 30.0 ;
    static const double RcutmaxZ = 0.7 ;
    static const double RcutmaxPt2 = 1.2 ;   
    static const double Rcutminnu = 4.0 ;
    static const double Rcutmaxnu = 9.0 ;
    static const int Rbin = 100;
    static const int Rbin_nu = 6;
    static const int Rbin_z  = 6;
    static const int Rbin_pt2= 6;

    static const int Dptbin_nu  = 5;
    static const int Dptbin_z   = 5;
    static const int Dptbin_Q = 5;
    static const int Dptbin_x = 5;

    //binning & others 4 ratios 
    static const int Cratiobin = 5;
    static const int Cratiobin_nu = 5;
    static const int Cratiobin_z = 5;
    static const int Cratiobin_Q = 5;
    static const int Cratiobin_x = 5;
    static const int Cratiobin_phih = 5;
    static const double numinCratio = 4;
    static const double numaxCratio = 9;
    
    
    
    static const double numinDpt = 4;
    static const double numaxDpt = 9;
    static const double xminDpt = 0;
    static const double xmaxDpt = 0;
    static const double phihminDpt = 0;
    static const double phihmaxDpt = 360;


    //cuts for Target specific vertex 
    static const double RcutminVzLD2 = -7.5 ;   
    static const double RcutmaxVzLD2 = -2.5 ;   
    static const double RcutminVzSn = -4 ;
    static const double RcutmaxVzSn = -1 ;
    static const double RcutminVzCu = -9 ;    //arbitrary +-.5
    static const double RcutmaxVzCu = -6 ;    //arbitrary +-.5
    static const double RcutminVzC = -7.5 ;     //default (for LD2)
    static const double RcutmaxVzC = -2.5 ;     //default (for LD2)
    static const double RcutminVzC1 = -4 ;    
    static const double RcutmaxVzC1 = -1 ;
    static const double RcutminVzC2 = -9 ;
    static const double RcutmaxVzC2 = -6 ;
    

    static const double RcutminVzLD2sim = -7.5 ;   
    static const double RcutmaxVzLD2sim = -2.5 ;   
    static const double RcutminVzSnsim = -4 ;
    static const double RcutmaxVzSnsim = -1 ;
    static const double RcutminVzCusim = -9 ;    //arbitrary +-.5
    static const double RcutmaxVzCusim = -6 ;    //arbitrary +-.5
    static const double RcutminVzCsim = -7.5 ;     //default (for LD2)
    static const double RcutmaxVzCsim = -2.5 ;     //default (for LD2)
    static const double RcutminVzC1sim = -4 ;    
    static const double RcutmaxVzC1sim = -1 ;
    static const double RcutminVzC2sim = -9 ;
    static const double RcutmaxVzC2sim = -6 ;

    static const double RcutminVzLD2data = -8.5 ;   
    static const double RcutmaxVzLD2data = -3.5 ;   
    static const double RcutminVzSndata = -5 ;
    static const double RcutmaxVzSndata = -2 ;
    static const double RcutminVzCudata = -10 ;    //arbitrary +-.5
    static const double RcutmaxVzCudata = -7 ;    //arbitrary +-.5
    static const double RcutminVzC1data = -5 ;    
    static const double RcutmaxVzC1data = -2 ;
    static const double RcutminVzC2data = -10 ;
    static const double RcutmaxVzC2data = -7 ;


    //RGA cuts 
    static const double cutlv_min = 14.0;    //in cm ?  cuts: loose >9cm ; medium >14cm; tight >19cm  
    static const double cutlv_max = 1000.0;
    static const double cutlw_min = 14.0;    //in cm ? 
    static const double cutlw_max = 1000.0;
    static const double cutlu_min = 0.0;    //doen't seem to have cuts in RGA
    static const double cutlu_max = 10000.0;  

    //Calo Cuts
    static const double cutEpcal_min = 0.06; //see 8.1.1 in RGA, in DST unit in GeV  (here is 60 MeV -> 0.06 GeV)
    
    static const double cutcalEp_max = 2.00; //see 8.1.1 in RGA, in DST unit in GeV  (here is 60 MeV)

    //Chi2 cuts
    static const double cutelchi2_min = 0.0;    // TBD !!
    static const double cutelchi2_max = 100.0;  // TBD !!        
}

#endif 