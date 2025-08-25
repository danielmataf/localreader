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
    //NewVertexValues passv11   (OB ?)  
    static const double v11cutminVzLD2data = -8.5 ;   
    static const double v11cutmaxVzLD2data = -4.5 ;   
    static const double v11cutminVzSndata = -4.5 ;
    static const double v11cutmaxVzSndata = -2 ;
    static const double v11cutminVzCudata = -8.5 ;    //arbitrary +-.5
    static const double v11cutmaxVzCudata = -7 ;    //arbitrary +-.5
    static const double v11cutminVzC1data = -4.5 ;    
    static const double v11cutmaxVzC1data = -2.5 ;
    static const double v11cutminVzC2data = -9.5 ;
    static const double v11cutmaxVzC2data = -7.5 ;    
    //NewVertexValues pass1   HADRONS!! (OB)  MAthieu values in pass1 
//C1    -9.67505	-5.38745
//C2    -4.50786	-0.68064
//Cu    -11.03832	-5.85492
//Sn    -4.74283	-0.42643
    
    static const double v11cutminVzLD2datahadron = -15 ;   //keeping these as the same...email below
    static const double v11cutmaxVzLD2datahadron = -5 ;   //keeping these as the same...email below
    static const double v11cutminVzSndatahadron = -4.74283 ;
    static const double v11cutmaxVzSndatahadron = -0.42643 ;
    static const double v11cutminVzCudatahadron = -11.03832 ;    //arbitrary +-.5
    static const double v11cutmaxVzCudatahadron = -5.85492 ;    //arbitrary +-.5
    static const double v11cutminVzC1datahadron = -4.50786 ;    
    static const double v11cutmaxVzC1datahadron = -0.68064 ;
    static const double v11cutminVzC2datahadron = -9.67505 ;
    static const double v11cutmaxVzC2datahadron = -5.38745 ;    

    //NewVertexValues pass1   ELECTRONS!! (OB)  MAthieu values in email
    static const double v11cutminVzLD2dataelectron = -15 ;   //keeping these as the same...
    static const double v11cutmaxVzLD2dataelectron = -5 ;   //keeping these as the same...
    static const double v11cutminVzSndataelectron = -6.21 ;
    static const double v11cutmaxVzSndataelectron = -1.2 ;
    static const double v11cutminVzCudataelectron = -10.55 ;    
    static const double v11cutmaxVzCudataelectron = -6.62 ;    
    static const double v11cutminVzC2dataelectron = -10.56 ;//merging both carbons, no separation.!!
    static const double v11cutmaxVzC2dataelectron = -5 ;    //merging both carbons, no separation.!!    



        //==NewBinings 
    //bining xb,Q with theta assist (and W)
    static const double lowA =  1.017;    
    static const double highA = 3.122;
    static const double lowB =  1;
    static const double highB = 2.38;
    static const double lowC =  1;
    static const double highC = 1.759;
    static const double lowD =  2.38;
    static const double highD = 5.073;
    static const double lowE =  1.759;
    static const double highE = 3.368;
    static const double lowF =  1;
    static const double highF = 2.246;
    static const double lowG =  3.368;
    static const double highG = 5;
    static const double lowH =  2.246;
    static const double highH = 3.941;
    static const double lowI =  1.467;
    static const double highI = 2.487;
    static const double lowJ =  3.941;
    static const double highJ = 5;
    static const double lowK =  2.449;
    static const double highK = 4.168;


        //LARGER    bining xb,Q with theta assist (and W)
    static const double LARGEQlowA = 1.128313968;
    static const double LARGEQhighA = 1.259832756;
    static const double LARGEQlowB = 1;
    static const double LARGEQhighB = 1.259832756;
    static const double LARGEQlowC = 1.549536766;
    static const double LARGEQhighC = 2.635206776;
    static const double LARGEQlowD = 1.259832756;
    static const double LARGEQhighD = 2.098986509;
    static const double LARGEQlowE = 1;
    static const double LARGEQhighE = 1.600454801;
    static const double LARGEQlowF = 2.098986509;
    static const double LARGEQhighF = 3.389715788;
    static const double LARGEQlowG = 1.600454801;
    static const double LARGEQhighG = 2.551322891;
    static const double LARGEQlowH = 1;
    static const double LARGEQhighH = 1.850633467;
    static const double LARGEQlowI = 2.551322891;
    static const double LARGEQhighI = 4.09282768;
    static const double LARGEQlowJ = 1.850633467;
    static const double LARGEQhighJ = 2.930202024;
    static const double LARGEQlowK = 1;
    static const double LARGEQhighK = 2.042169211;
    static const double LARGEQlowL = 2.930202024;
    static const double LARGEQhighL = 4.875821862;
    static const double LARGEQlowM = 2.042169211;
    static const double LARGEQhighM = 3.310850781;
    static const double LARGEQlowN = 1.039050613;
    static const double LARGEQhighN = 2.220055645;
    static const double LARGEQlowO = 3.310850781;
    static const double LARGEQhighO = 5;
    static const double LARGEQlowP = 2.220055645;
    static const double LARGEQhighP = 3.941119738;
    static const double LARGEQlowQ = 1.400459522;
    static const double LARGEQhighQ = 2.486714814;
    static const double LARGEQlowR = 3.941119738;
    static const double LARGEQhighR = 5;
    static const double LARGEQlowS = 2.449190731;
    static const double LARGEQhighS = 4.483907471;

    //large limits for xb 
    static const double LargexblowAB = 0.06;
    static const double LargexbhighAB = 0.1;
    static const double LargexbhighCDE = 0.15;
    static const double LargexbhighFGH = 0.2;
    static const double LargexbhighIJK = 0.25;
    static const double LargexbhighLMN = 0.31;
    static const double LargexbhighOPQ = 0.44;
    static const double LargexbhighRS = 0.6;
    


        //LARGER    bining xb,Q with theta assist (and W)
    static const double FEWQlowA = 1;
    static const double FEWQhighA = 2.291;
    static const double FEWQlowB = 1;
    static const double FEWQhighB = 1.549;
    static const double FEWQlowC = 1.549;
    static const double FEWQhighC = 2.908;
    static const double FEWQlowD = 1.197;
    static const double FEWQhighD = 1.808;
    static const double FEWQlowE = 1;
    static const double FEWQhighE = 1.349;
    static const double FEWQlowF = 1.808;
    static const double FEWQhighF = 4.176;
    static const double FEWQlowG = 1.346;
    static const double FEWQhighG = 2.229;
    static const double FEWQlowH = 1;
    static const double FEWQhighH = 1.566;
    static const double FEWQlowI = 2.229;
    static const double FEWQhighI = 4.937;
    static const double FEWQlowJ = 1.566;
    static const double FEWQhighJ = 2.429;
    static const double FEWQlowK = 1;
    static const double FEWQhighK = 1.663;
    static const double FEWQlowL = 2.429;
    //static const double LARGEQhighL = 4.875821862;

    //large limits for xb 
    static const double FewxblowAB = 0.06;
    static const double FewxbhighAB = 0.13;
    static const double FewxbhighCDE = 0.17;
    static const double FewxbhighFGH = 0.26;
    static const double FewxbhighIJK = 0.32;
    static const double FewxbhighLMN = 0.32;
    //Pass1 Vz mean values electron
    static const double meanC1ibel = -7.38522; 
    static const double meanC1obel = -3.66242;
    static const double meanC2ibel = -2.49464;
    static const double meanC2obel = -8.53079;
    static const double meanLD2ibel = -5.0;
    static const double meanLD2obel = -5.0;
    static const double meanSnibel = -2.0;
    static const double meanSnobel = -2.0;
    static const double meanCuibel = -8.49862;
    static const double meanCuobel = -3.63868;

    //Pass1 Vz mean values hadron pi+
    static const double meanC1ibhad = -3.36609;
    static const double meanC1obhad = -2.59425;
    static const double meanC2ibhad = -8.28039;
    static const double meanC2obhad = -7.53125;
    static const double meanLD2ibhad = -5.0;
    static const double meanLD2obhad = -5.0;
    static const double meanSnibhad = -2.0;
    static const double meanSnobhad = -2.0;
    static const double meanCuibhad = -8.0;
    static const double meanCuobhad = -8.0;





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
    static const double Rcutminx = 0.05 ;
    static const double RcutminY = 0.25 ;
    static const double RcutminW = 4.0 ;
    static const double RcutminZ = 0.3 ;
    static const double RcutminPt2 = 0.0 ;
    static const double RcutmaxQ = 5.0 ;
    static const double Rcutmaxx = 0.5 ;
    static const double RcutmaxY = 0.85 ;
    static const double RcutmaxW = 30.0 ;
    static const double RcutmaxZ = 0.7 ;
    static const double RcutmaxPt2 = 1.2 ;   
    static const double Rcutminnu = 2.5 ;
    static const double Rcutmaxnu = 9.0 ;
    //adding phih cuts
    static const double Rcutminphih = 0.0 ;
    static const double Rcutmaxphih = 360.0 ;

    //Adding below cut for nu to avoid acceptance corrections
    static const double RcutnuMin = 0.0; // avoid acceptance corrections 
    static const double RcutnuMax = 7.0; // avoid acceptance corrections
        //for the name I added the R to keep consistency with the other cuts used in general .  

    static const int Rbin = 100;
    static const int Rbin_Q = 6;
    static const int Rbin_x = 6;
    static const int Rbin_nu = 6;
    static const int Rbin_z  = 6;
    static const int Rbin_pt2= 6;

    static const int Dptbin_nu  = 6;
    static const int Dptbin_z   = 6;
    static const int Dptbin_Q = 6;
    static const int Dptbin_x = 6;

    //binning & others 4 ratios 
    static const int Cratiobin = 6;
    static const int Cratiobin_nu = 6;
    static const int Cratiobin_z = 6;
    static const int Cratiobin_Q = 6;
    static const int Cratiobin_x = 6;
    static const int Cratiobin_phih = 6;
    static const double numinCratio = 6;
    static const double numaxCratio = 6;
    
    
    
    //static const double numinDpt = 4;
    static const double numinDpt = 2.8;
    static const double numaxDpt = 9;
    static const double xminDpt = 0.07;
    static const double xmaxDpt = 0.5;
    static const double phihminDpt = 0;
    static const double phihmaxDpt = 360;


    //cuts for Target specific vertex 
    static const double RcutminVzLD2 = -7.5 ;   
    static const double RcutmaxVzLD2 = -2.5 ;   
    static const double RcutminVzSn = -4 ;
    static const double RcutmaxVzSn = -1 ;
    static const double RcutminVzCu = -10 ;    //arbitrary +-.5
    static const double RcutmaxVzCu = -7 ;    //arbitrary +-.5
    static const double RcutminVzC = -7.5 ;     //default (for LD2)
    static const double RcutmaxVzC = -2.5 ;     //default (for LD2)
    static const double RcutminVzC1 = -4 ;    
    static const double RcutmaxVzC1 = -1 ;
    static const double RcutminVzC2 = -9 ;
    static const double RcutmaxVzC2 = -8 ;
    

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