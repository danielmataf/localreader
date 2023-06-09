#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

#include "TH1F.h"

int nubin = 10;  //--------------------------------------------------------------------binning -->USE 10 FOR RATIO PLOTS
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


struct canvaslist
{
    TCanvas* cnew; 
    canvaslist(){
        cnew=new TCanvas("new","new"); //create new canvas



    }
    ~canvaslist(){
        delete cnew;

    }
};

struct cutlist
{
    double cut1;
    double cut2; 
    double cutQ;
    double cutY;
    double cutV;
    double cutW;
    double cutZmin;
    double cutZmax;

    cutlist(){
        cut1 = 5;   //test 
        cut2 = 10;  //test
        cutQ = 1.5;
        cutY = 0.85;
        cutV = 8.5;
        cutW = 3.16; //sqrt(10)
        cutZmin = 0.3;
        cutZmax = 0.7;
    }
    ~cutlist(){

    }
};



canvaslist myCanvas;
cutlist myCuts;


/*
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
*/


//creation of a class for kinematical variables and others
//no need to differentiate between the nuclea target, main should do that distinction, and separate it when writting output files.

class VarHistograms {
public:
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
        hist_phih = new TH1F("phih", "phih", 12, 0, 360);
        hist_t = new TH1F("t", "t", nubin, -15, 0);
        hist_P_t = new TH1F("P_t", "P_t", nubin, 0, 2);
        hist_z = new TH1F("z", "z", nubin, zminX, zmaxX);
        hist_MX = new TH1F("MX", "MX", nubin, 0, 20);
        hist_MM_epi = new TH1F("MM_epi", "MM_epi", nubin, 0, 20);
        hist_MM_eX = new TH1F("MM_eX", "MM_eX", nubin, 0, 20);
        hist_MM_piX = new TH1F("MM_piX", "MM_piX", nubin, 0, 20);
        hist_MM_epiX = new TH1F("MM_epiX", "MM_epiX", nubin, 0, 20);
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

//myHists.hist_Q->Write();

//maybe no need to distinguish between REC and MC 

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

/*
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
*/  

class hist_hadX {
public:
    TH1F* hpx_X;
    TH1F* hpy_X;
    TH1F* hpz_X;
    TH1F* h_P_trX;
    TH2F* MC_coor2D_missX;
    TH2F* coor2D_missX;
    
    hist_hadX() {
        hpx_X = new TH1F("px_X", "px_X", nubin, -5, 5);
        hpy_X = new TH1F("py_X", "py_X", nubin, -5, 5);
        hpz_X = new TH1F("pz_X", "pz_X", nubin, 0, 10);
        h_P_trX = new TH1F("P_trX", "P_trX", nubin, 0, 2);
        coor2D_missX = new TH2F("thvsphi_missX", "thvsphi_missX", nubin, -180, 180, nubin, 0, 60);
        MC_coor2D_missX = new TH2F("MCthvsphi_missX", "MCthvsphi_missX", nubin, -180, 180, nubin, 0, 60);
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

/*
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
*/
class hist_pip {
public:
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
        hID_pip = new TH1F("ID_pip", "ID_pip", nubin, -5, 5);
        h_MCtheta_pip = new TH1F("MCtheta_pip", "MCtheta_pip", nubin, 0, 90);
        h_MCphi_pip = new TH1F("MCphi_pip", "MCphi_pip", nubin, -180, 180);
        h_theta_pip = new TH1F("theta_pip", "theta_pip", nubin, 0, 90);
        h_phi_pip = new TH1F("phi_pip", "phi_pip", nubin, -180, 180);
        coor2D_pip = new TH2F("thvsphi_pho", "thvsphi_pho", nubin, 0, 180, nubin, 0, 60);
        MC_coor2D_pip = new TH2F("MC_thvsphi_pip", "MC_thvsphi_pip", nubin, 0, 180, nubin, 0, 60);
        hpx_pi = new TH1F("px_pi", "px_pi", nubin, -5, 5);
        hpy_pi = new TH1F("py_pi", "py_pi", nubin, -5, 5);
        hpz_pi = new TH1F("pz_pi", "pz_pi", nubin, 0, 10);
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



/*
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
    h_scalE = new TH1F("nrg(el)", "nrg(el)", nubin, 0.0, 10);
    h_MCscalE = new TH1F("MCnrg(el)", "MCnrg(el)", nubin, 0.0, 10);
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
        delete h_scalE;
        delete h_MCscalE;
    
    }
};
*/
class hist_electron {
public:
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
        hID_el = new TH1F("ID_el", "ID_el", 1000, -5, 5);
        h_MCtheta_el = new TH1F("MCtheta_el", "MCtheta_el", nubin, 0, 90);
        h_MCphi_el = new TH1F("MCphi_el", "MCphi_el", nubin, -180, 180);
        h_theta_el = new TH1F("theta_el", "theta_el", nubin, 0, 90);
        h_phi_el = new TH1F("phi_el", "phi_el", nubin, -180, 180);
        coor2D_el = new TH2F("thvsphi_EL", "thvsphi_EL", nubin, -180, 180, nubin, 0, 60);
        MC_coor2D_el = new TH2F("MC_thvsphi_EL", "MC_thvsphi_EL", nubin, -180, 180, nubin, 0, 60);
        hpx_el = new TH1F("px_el", "px_el", nubin, -5, 5);
        hpy_el = new TH1F("py_el", "py_el", nubin, -5, 5);
        hpz_el = new TH1F("pz_el", "pz_el", nubin, 0, 10);
        h_scalE = new TH1F("nrg(el)", "nrg(el)", nubin, 0.0, 10);
        h_MCscalE = new TH1F("MCnrg(el)", "MCnrg(el)", nubin, 0.0, 10);
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
        delete h_scalE;
        delete h_MCscalE;
    }
};



struct linelist {

    TLine *cutMMX_Sn;
    TLine *line1; //test
    TLine *lineQ;
    TLine *lineY;
    TLine *lineV;
    TLine *lineW;
    TLine *lineZmin;
    TLine *lineZmax;

    linelist(){
        line1 = new TLine(myCuts.cut1, 0, myCuts.cut1, 0);
        lineQ = new TLine(myCuts.cutQ, 0, myCuts.cut1, 0);
        lineY = new TLine(myCuts.cutY, 0, myCuts.cut1, 0);
        lineV = new TLine(myCuts.cutV, 0, myCuts.cut1, 0);
        lineW = new TLine(myCuts.cutW, 0, myCuts.cut1, 0);
        lineZmin = new TLine(myCuts.cutZmin, 0, myCuts.cut1, 0);
        lineZmax = new TLine(myCuts.cutZmax, 0, myCuts.cut1, 0);
        //can I call a histogram minimum ?  in order to place the Cut. Or call any specific value for that matter   
        //the specific value of x-axis is given by the cut value. need to give a value of y-axis for the line to be drawn 
        //cutMMX_Sn = new TLine(1.2,cnew->cd(9)->GetUymin(),1.2,cnew->cd(9)->GetUymax());
    }


    ~linelist(){
        delete line1;
        delete lineQ;
        delete lineY;
        delete lineV;
        delete lineW;
        delete lineZmax;
        delete lineZmin;
    }

};

//following class multipratio uses all the histograms needed for variables in the multiplicity ratio
//Using one histogram per Nucleus.
//Also using one histogram per variable.... not ideal, maybe more iseful to stock histograms in a vector of variable_histogram for each Nucleus 

/*
    TH1F* hist_Q_Sn_e = new TH1F("Q2_Sn_e", "Q2_Sn_e", nubin, QminX, QmaxX);
    TH1F* hist_Q_De_e= new TH1F("Q2_De_e", "Q2_De_e", nubin, QminX, QmaxX);
    TH1F* hist_v_Sn_e = new TH1F("v_Sn_e", "v_Sn_e", nubin, vminX, vmaxX);
    TH1F* hist_v_De_e = new TH1F("v_De_e", "v_De_e", nubin, vminX, vmaxX);
    TH1F* hist_z_Sn_e = new TH1F("z_Sn_e", "z_Sn_e", nubin, zminX, zmaxX);
    TH1F* hist_z_De_e = new TH1F("z_De_e", "z_De_e", nubin, zminX, zmaxX);
    TH1F* hist_P_t_Sn_e = new TH1F("pt_Sn_e", "pt_Sn_e", nubin, PtminX, PtmaxX);
    TH1F* hist_P_t_De_e = new TH1F("pt_De_e", "pt_De_e", nubin, PtminX, PtmaxX);

    TH1F* h_meancos_Sn = new TH1F("meancos_Sn", "meancos_Sn", phibin, -1, 1);
    TH1F* h_meancos_De = new TH1F("meancos_De", "meancos_De", phibin, -1,1);
    TH1F* h_cosSn = new TH1F("h_cosSn", "h_cosSn", 10, -1,1);
    TH1F* h_cosDe = new TH1F("h_cosDe", "h_cosDe", 10, -1,1);
*/

class MultipRatio {     //Multipratio graphes, use them !!!!
public:
    TH1F* hist_Q_e;
    TH1F* hist_v_e ;
    TH1F* hist_z_e;
    TH1F* hist_P_t_e;
    

    TGraphErrors* R_Sn_Q;
    TGraphErrors* R_Sn_v;
    TGraphErrors* R_Sn_z;
    TGraphErrors* R_Sn_pt;
    TGraphErrors* R_D_Q;
    TGraphErrors* R_D_v;
    TGraphErrors* R_D_z;
    TGraphErrors* R_D_pt;
    TGraphErrors* R_Cu_Q;
    TGraphErrors* R_Cu_v;
    TGraphErrors* R_Cu_z;
    TGraphErrors* R_Cu_pt;

    MultipRatio() {
        hist_Q_e = new TH1F("onlye_Q", "onlye_Q", nubin, QminX ,QmaxX);
        hist_v_e = new TH1F("onlye_v", "onlye_v", nubin, vminX ,vmaxX);
        hist_z_e = new TH1F("onlye_z", "onlye_z", nubin, zminX ,zmaxX);
        hist_P_t_e = new TH1F("onlye_pt", "onlye_pt", nubin, PtminX ,PtmaxX);




        R_Sn_Q = new TGraphErrors();
        R_Sn_v = new TGraphErrors();
        R_Sn_z = new TGraphErrors();
        R_Sn_pt = new TGraphErrors();
        R_D_Q = new TGraphErrors();
        R_D_v = new TGraphErrors();
        R_D_z = new TGraphErrors();
        R_D_pt = new TGraphErrors();
        R_Cu_Q = new TGraphErrors();
        R_Cu_v = new TGraphErrors();
        R_Cu_z = new TGraphErrors();
        R_Cu_pt = new TGraphErrors();
    }

    ~MultipRatio() {
        delete R_Sn_Q ;
        delete R_Sn_v  ;
        delete R_Sn_z ;
        delete R_Sn_pt;
        delete R_D_Q  ;
        delete R_D_v  ;
        delete R_D_z ;
        delete R_D_pt ;
        delete R_Cu_Q ;
        delete R_Cu_v ;
        delete R_Cu_z ;
        delete R_Cu_pt ;

        delete hist_Q_e;
        delete hist_v_e ;
        delete hist_z_e ;
        delete hist_P_t_e;

    }
};


//Cratio determines the cosinus ratio of the three -eventual- nuclei but follows same principle as with multipratio
//can be improved
class Cratio {                  //Cratio graphes!!!! use them 
public:
    TGraphErrors* C_Sn_Q;
    TGraphErrors* C_Sn_v;
    TGraphErrors* C_Sn_z;
    TGraphErrors* C_Sn_pt;
    TGraphErrors* C_D_Q;
    TGraphErrors* C_D_v;
    TGraphErrors* C_D_z;
    TGraphErrors* C_D_pt;
    TGraphErrors* C_Cu_Q;
    TGraphErrors* C_Cu_v;
    TGraphErrors* C_Cu_z;
    TGraphErrors* C_Cu_pt;

    Cratio() {
        C_Sn_Q = new TGraphErrors();
        C_Sn_v = new TGraphErrors();
        C_Sn_z = new TGraphErrors();
        C_Sn_pt = new TGraphErrors();
        C_D_Q = new TGraphErrors();
        C_D_v = new TGraphErrors();
        C_D_z = new TGraphErrors();
        C_D_pt = new TGraphErrors();
        C_Cu_Q = new TGraphErrors();
        C_Cu_v = new TGraphErrors();
        C_Cu_z = new TGraphErrors();
        C_Cu_pt = new TGraphErrors();
    }

    ~Cratio() {
        delete C_Sn_Q ;
        delete C_Sn_v  ;
        delete C_Sn_z ;
        delete C_Sn_pt;
        delete C_D_Q  ;
        delete C_D_v  ;
        delete C_D_z ;
        delete C_D_pt ;
        delete C_Cu_Q ;
        delete C_Cu_v ;
        delete C_Cu_z ;
        delete C_Cu_pt ;
    }
};

//Dpt following the same principle as before 
//Cratio graphes!!!! use them 
class g_Dpt {
public:
    TGraphErrors* Dpt_Sn_Q;
    TGraphErrors* Dpt_Sn_v;
    TGraphErrors* Dpt_Sn_z;
    TGraphErrors* Dpt_D_Q;
    TGraphErrors* Dpt_D_v;
    TGraphErrors* Dpt_D_z;
    TGraphErrors* Dpt_Cu_Q;
    TGraphErrors* Dpt_Cu_v;
    TGraphErrors* Dpt_Cu_z;

    g_Dpt() {
        Dpt_Sn_Q = new TGraphErrors();
        Dpt_Sn_v = new TGraphErrors();
        Dpt_Sn_z = new TGraphErrors();
        Dpt_D_Q = new TGraphErrors();
        Dpt_D_v = new TGraphErrors();
        Dpt_D_z = new TGraphErrors();
        Dpt_Cu_Q = new TGraphErrors();
        Dpt_Cu_v = new TGraphErrors();
        Dpt_Cu_z = new TGraphErrors();
    }

    ~g_Dpt() {
        delete Dpt_Sn_Q ;
        delete Dpt_Sn_v  ;
        delete Dpt_Sn_z ;
        delete Dpt_D_Q  ;
        delete Dpt_D_v  ;
        delete Dpt_D_z ;
        delete Dpt_Cu_Q ;
        delete Dpt_Cu_v ;
        delete Dpt_Cu_z ;
    }
};



#include <TH1F.h>

class h_Dpt {
public:
    TH1F* w_ptQ_Sn;
    TH1F* w_ptQ_De;
    TH1F* w_ptv_Sn;
    TH1F* w_ptv_De;
    TH1F* w_ptz_Sn;
    TH1F* w_ptz_De;

    h_Dpt() {
        w_ptQ_Sn = new TH1F("weightedptQSn", "weightedptQSn", nubin, QminX, QmaxX);
        w_ptQ_De = new TH1F("weightedptQDe", "weightedptQDe", nubin, QminX, QmaxX);
        w_ptv_Sn = new TH1F("weightedptvSn", "weightedptvSn", nubin, vminX, vmaxX);
        w_ptv_De = new TH1F("weightedptvDe", "weightedptvDe", nubin, vminX, vmaxX);
        w_ptz_Sn = new TH1F("weightedptzSn", "weightedptzSn", nubin, zminX, zmaxX);
        w_ptz_De = new TH1F("weightedptzDe", "weightedptzDe", nubin, zminX, zmaxX);

        //w_ptQ_Sn->Sumw2();
        //w_ptQ_De->Sumw2();
        //w_ptv_Sn->Sumw2();
        //w_ptv_De->Sumw2();
        //w_ptz_Sn->Sumw2();
        //w_ptz_De->Sumw2();
    }

    ~h_Dpt() {
        delete w_ptQ_Sn;
        delete w_ptQ_De;
        delete w_ptv_Sn;
        delete w_ptv_De;
        delete w_ptz_Sn;
        delete w_ptz_De;
    }
};

class h_cratio
{
public:
    TH1F* w_cQ_Sn ;
    TH1F* w_cQ_De ;
    TH1F* w_cv_Sn ;
    TH1F* w_cv_De ;
    TH1F* w_cz_Sn ;
    TH1F* w_cz_De ;
    TH1F* w_cpt_Sn;
    TH1F* w_cpt_De;
    h_cratio(){
         w_cQ_Sn = new TH1F("weightedcQSn", "weightedcQSn", nubin, QminX ,QmaxX);
         w_cQ_De = new TH1F("weightedcQDe", "weightedcQDe", nubin, QminX ,QmaxX);
         w_cv_Sn = new TH1F("weightedcvSn", "weightedcvSn", nubin, vminX ,vmaxX);
         w_cv_De = new TH1F("weightedcvDe", "weightedcvDe", nubin, vminX ,vmaxX);
         w_cz_Sn = new TH1F("weightedczSn", "weightedczSn", nubin, zminX ,zmaxX);
         w_cz_De = new TH1F("weightedczDe", "weightedczDe", nubin, zminX ,zmaxX);
         w_cpt_Sn = new TH1F("weightedcptSn", "weightedcptSn", nubin, PtminX ,PtmaxX);
         w_cpt_De = new TH1F("weightedcptDe", "weightedcptDe", nubin, PtminX ,PtmaxX);
    }
    ~h_cratio(){
        delete w_cQ_Sn ;
        delete w_cQ_De ;
        delete w_cv_Sn ;
        delete w_cv_De ;
        delete w_cz_Sn ;
        delete w_cz_De ;
        delete w_cpt_Sn;
        delete w_cpt_De;
    }

};


class h_sratio
{
public:
    TH1F* w_sQ_Sn ;
    TH1F* w_sQ_De ;
    TH1F* w_sv_Sn ;
    TH1F* w_sv_De ;
    TH1F* w_sz_Sn ;
    TH1F* w_sz_De ;
    TH1F* w_spt_Sn;
    TH1F* w_spt_De;
    h_sratio(){
         w_sQ_Sn = new TH1F("weightedsQSn", "weightedsQSn", nubin, QminX ,QmaxX);
         w_sQ_De = new TH1F("weightedsQDe", "weightedsQDe", nubin, QminX ,QmaxX);
         w_sv_Sn = new TH1F("weightedsvSn", "weightedsvSn", nubin, vminX ,vmaxX);
         w_sv_De = new TH1F("weightedsvDe", "weightedsvDe", nubin, vminX ,vmaxX);
         w_sz_Sn = new TH1F("weightedszSn", "weightedszSn", nubin, zminX ,zmaxX);
         w_sz_De = new TH1F("weightedszDe", "weightedszDe", nubin, zminX ,zmaxX);
         w_spt_Sn = new TH1F("weightedsptSn", "weightedsptSn", nubin, PtminX ,PtmaxX);
         w_spt_De = new TH1F("weightedsptDe", "weightedsptDe", nubin, PtminX ,PtmaxX);
    }
    ~h_sratio(){
        delete w_sQ_Sn ;
        delete w_sQ_De ;
        delete w_sv_Sn ;
        delete w_sv_De ;
        delete w_sz_Sn ;
        delete w_sz_De ;
        delete w_spt_Sn;
        delete w_spt_De;
    }

};


class h_c2ratio
{
public:
    TH1F* w_c2Q_Sn ;
    TH1F* w_c2Q_De ;
    TH1F* w_c2v_Sn ;
    TH1F* w_c2v_De ;
    TH1F* w_c2z_Sn ;
    TH1F* w_c2z_De ;
    TH1F* w_c2pt_Sn;
    TH1F* w_c2pt_De;
    h_c2ratio(){
         w_c2Q_Sn = new TH1F("weightedc2QSn", "weightedc2QSn", nubin, QminX ,QmaxX);
         w_c2Q_De = new TH1F("weightedc2QDe", "weightedc2QDe", nubin, QminX ,QmaxX);
         w_c2v_Sn = new TH1F("weightedc2vSn", "weightedc2vSn", nubin, vminX ,vmaxX);
         w_c2v_De = new TH1F("weightedc2vDe", "weightedc2vDe", nubin, vminX ,vmaxX);
         w_c2z_Sn = new TH1F("weightedc2zSn", "weightedc2zSn", nubin, zminX ,zmaxX);
         w_c2z_De = new TH1F("weightedc2zDe", "weightedc2zDe", nubin, zminX ,zmaxX);
         w_c2pt_Sn = new TH1F("weightedc2ptSn", "weightedc2ptSn", nubin, PtminX ,PtmaxX);
         w_c2pt_De = new TH1F("weightedc2ptDe", "weightedc2ptDe", nubin, PtminX ,PtmaxX);
    }
    ~h_c2ratio(){
        delete w_c2Q_Sn ;
        delete w_c2Q_De ;
        delete w_c2v_Sn ;
        delete w_c2v_De ;
        delete w_c2z_Sn ;
        delete w_c2z_De ;
        delete w_c2pt_Sn;
        delete w_c2pt_De;
    }

};

#endif