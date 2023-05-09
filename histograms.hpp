#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

#include "TH1F.h"

int nubin = 100;  //--------------------------------------------------------------------binning -->USE 10 FOR RATIO PLOTS
int phibin= 10;
int Rbin = 10;
int QminX = 0;
int QmaxX = 6;
int vminX = 0;
int vmaxX = 11;
int zminX = 0;
int zmaxX = 1;
int PtminX = 0;
int PtmaxX = 2;

struct VarHistograms {
    TH1F* hist_Q;
    TH1F* hist_v;
    TH1F* hist_y;
    TH1F* hist_x;
    TH1F* hist_W;
    TH1F* hist_phih;
    TH1F* hist_t;
    TH1F* hist_P_t;
    TH1F* hist_z;
    TH1F* hist_MX;
    TH1F* hist_MM_epi;
    TH1F* hist_MM_eX;
    TH1F* hist_MM_piX;
    TH1F* hist_MM_epiX;

    VarHistograms() {
        hist_Q = new TH1F("Q2", "Q2", nubin, QminX, QmaxX);
        hist_v = new TH1F("gamnu", "gamnu", nubin, vminX, vmaxX);
        hist_y = new TH1F("y", "y", nubin, 0.0, 1.0);
        hist_x = new TH1F("x_bj", "x_bj", nubin, 0.0, 0.6);
        hist_W = new TH1F("W", "W", nubin, 1.5, 5);
        hist_phih = new TH1F("phih", "phih", 12, 0,360);
        hist_t = new TH1F("t", "t", nubin, -15,0);
        hist_P_t = new TH1F("P_t", "P_t", nubin, 0,2);
        hist_z = new TH1F("z", "z", nubin, zminX,zmaxX);
        hist_MX = new TH1F("MX", "MX", nubin, 0,20);
        hist_MM_epi = new TH1F("MM_epi", "MM_epi", nubin, 0,20);
        hist_MM_eX = new TH1F("MM_eX", "MM_eX", nubin, 0,20);
        hist_MM_piX = new TH1F("MM_piX", "MM_piX", nubin, 0,20);
        hist_MM_epiX = new TH1F("MM_epiX", "MM_epiX", nubin, 0,20);
    }

    ~VarHistograms() {
        delete hist_Q;
        delete hist_v;
        delete hist_y;
        delete hist_x;
        delete hist_W;
        delete hist_phih;
        delete hist_t;
        delete hist_P_t;
        delete hist_z;
        delete hist_MM_epi;
        delete hist_MM_eX;
        delete hist_MM_piX;
        delete hist_MM_epiX;
    }
};



struct MCHistograms {
    TH1F* hist_MC_Q;
    TH1F* hist_MC_v;
    TH1F* hist_MC_y;
    TH1F* hist_MC_x;
    TH1F* hist_MC_W;
    TH1F* hist_MC_phih;
    TH1F* hist_MC_t;
    TH1F* hist_MC_P_t;
    TH1F* hist_MC_z;
    TH1F* hist_MC_MX;
    TH1F* hist_MC_MM_epi;
    TH1F* hist_MC_MM_eX;
    TH1F* hist_MC_MM_piX;
    TH1F* hist_MC_MM_epiX;

    MCHistograms() {
    hist_MC_Q = new TH1F("MCQ2", "MCQ2", nubin, 0.0, 6);
    hist_MC_v = new TH1F("MCgamnu", "MCgamnu", nubin, 0.0, 11);
    hist_MC_y = new TH1F("MCy", "MCy", nubin, 0.0, 1.0);
    hist_MC_x = new TH1F("MCx_bj", "MCx_bj", nubin, 0.0, 0.6);
    hist_MC_W = new TH1F("MCW", "MCW", nubin, 1.5, 5);
    hist_MC_t = new TH1F("MCt", "MCt", nubin, -15,0);
    hist_MC_phih = new TH1F("MCphih", "MCphih", 12, 0,360);
    hist_MC_P_t = new TH1F("MC_P_t", "MC_P_t", nubin, 0,1);
    hist_MC_z = new TH1F("MC_z", "MC_z", nubin, 0,1);
    hist_MC_MX = new TH1F("MC_MX", "MC_MX", nubin, 0,20);
    hist_MC_MM_epi = new TH1F("MC_MM_epi", "MC_MM_epi", nubin, 0,20);
    hist_MC_MM_eX = new TH1F("MC_MM_eX", "MC_MM_eX", nubin, 0,20);
    hist_MC_MM_piX = new TH1F("MC_MM_piX", "MC_MM_piX", nubin, 0,20);
    hist_MC_MM_epiX = new TH1F("MC_MM_epiX", "MC_MM_epiX", nubin, 0,20);

    }

    ~MCHistograms() {
        delete hist_MC_Q;
        delete hist_MC_v;
        delete hist_MC_y;
        delete hist_MC_x;
        delete hist_MC_W;
        delete hist_MC_phih;
        delete hist_MC_t;
        delete hist_MC_P_t;
        delete hist_MC_z;
        delete hist_MC_MM_epi;
        delete hist_MC_MM_eX;
        delete hist_MC_MM_piX;
        delete hist_MC_MM_epiX;
    }
};



    
    //=========Notes============//
	//create also momenta for MC values							TBD !!!
    
    //==========================//


struct hist_hadX {

    TH1F* hpx_X;
    TH1F* hpy_X;
    TH1F* hpz_X;
    TH1F* h_P_trX;
    TH2F* MC_coor2D_missX;
    TH2F* coor2D_missX;
    
    hist_hadX() {
    hpx_X = new TH1F("px_X", "px_X", nubin, -5,5);			//these three lines are hypothetical hadron X using missing vector w/ epi
    hpy_X = new TH1F("py_X", "py_X", nubin, -5,5);
    hpz_X = new TH1F("pz_X", "pz_X", nubin, 0,10);
    h_P_trX = new TH1F("P_trX", "P_trX", nubin, 0,2);		//transverse momentum for X-HADRON
    coor2D_missX = new TH2F ("thvsphi_missX", "thvsphi_missX", nubin, -180, 180, nubin, 0,60);
    MC_coor2D_missX= new TH2F ("MCthvsphi_missX", "MCthvsphi_missX", nubin, -180, 180, nubin, 0,60);
    
    }

    ~hist_hadX() {
        delete hpx_X;
        delete hpy_X;
        delete hpz_X;
        delete h_P_trX;
        delete MC_coor2D_missX;
        delete coor2D_missX;
    }
};

struct hist_pip {
    TH1F* h_MCtheta_pip; 
    TH1F* h_MCphi_pip;
    TH1F* h_theta_pip;
    TH1F* h_phi_pip;
    TH2F* coor2D_pip;
    TH2F* MC_coor2D_pip;
    TH1F* hpx_pi;
    TH1F* hpy_pi;
    TH1F* hpz_pi;
    TH1F* hID_pip; 

    hist_pip() {
    hID_pip = new TH1F("ID_pip", "ID_pip", nubin, -5,5);
    h_MCtheta_pip = new TH1F("MCtheta_pip", "MCtheta_pip", nubin, 0,90);
    h_MCphi_pip = new TH1F("MCphi_pip", "MCphi_pip", nubin, -180,180);
    h_theta_pip = new TH1F("theta_pip", "theta_pip", nubin, 0,90);
    h_phi_pip = new TH1F("phi_pip", "phi_pip", nubin, -180,180);
    coor2D_pip = new TH2F ("thvsphi_pho", "thvsphi_pho", nubin, 0, 180, nubin, 0, 60);
    MC_coor2D_pip = new TH2F ("MC_thvsphi_pip", "MC_thvsphi_pip", nubin, 0, 180, nubin, 0, 60);
    hpx_pi = new TH1F("px_pi", "px_pi", nubin, -5,5);
    hpy_pi = new TH1F("py_pi", "py_pi", nubin, -5,5);
    hpz_pi = new TH1F("pz_pi", "pz_pi", nubin, 0,10);
    
    }

    ~hist_pip() {
        delete hID_pip;
        delete h_MCtheta_pip;
        delete h_MCphi_pip;
        delete h_theta_pip;
        delete h_phi_pip;
        delete coor2D_pip;
        delete MC_coor2D_pip;
        delete hpx_pi;
        delete hpy_pi;
        delete hpz_pi;
      
    
    }
};

struct hist_electron {
    TH1F* hID_el; 
    TH1F* h_MCtheta_el;
    TH1F* h_MCphi_el;
    TH1F* h_theta_el;
    TH1F* h_phi_el;
    TH2F* coor2D_el;
    TH2F* MC_coor2D_el;
    TH1F* hpx_el;
    TH1F* hpy_el;
    TH1F* hpz_el;
    TH1F* h_scalE;
    TH1F* h_MCscalE;
    hist_electron() {
    hID_el = new TH1F("ID_el", "ID_el", 1000, -5,5);
    h_MCtheta_el = new TH1F("MCtheta_el", "MCtheta_el", nubin, 0,90);
    h_MCphi_el = new TH1F("MCphi_el", "MCphi_el", nubin, -180,180);
    h_theta_el = new TH1F("theta_el", "theta_el", nubin, 0,90);
    h_phi_el = new TH1F("phi_el", "phi_el", nubin, -180,180);
    coor2D_el = new TH2F ("thvsphi_EL", "thvsphi_EL", nubin, -180, 180, nubin, 0,60);
    MC_coor2D_el = new TH2F ("MC_thvsphi_EL", "MC_thvsphi_EL", nubin, -180, 180, nubin, 0,60);
    hpx_el = new TH1F("px_el", "px_el", nubin, -5,5);
    hpy_el = new TH1F("py_el", "py_el", nubin, -5,5);
    hpz_el = new TH1F("pz_el", "pz_el", nubin, 0,10);
    TH1F* h_scalE = new TH1F("nrg(el)", "nrg(el)", nubin, 0.0, 10);
    TH1F* h_MCscalE = new TH1F("MCnrg(el)", "MCnrg(el)", nubin, 0.0, 10);
    }

    ~hist_electron() {
        delete hID_el;
        delete h_MCtheta_el;
        delete h_MCphi_el;
        delete h_theta_el;
        delete h_phi_el;
        delete coor2D_el;
        delete MC_coor2D_el;
        delete hpx_el;
        delete hpy_el;
        delete hpz_el;
    
    }
};









#endif