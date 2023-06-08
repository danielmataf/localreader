
#include "reader.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TCanvas.h>
#include <TH2.h>
#include <TF1.h>
#include <TRatioPlot.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TLine.h>
#include <TArrow.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define PI 3.14159265
//#include "dcFiducialCuts.hpp"
#include "linspace.hpp"
#include "variablecalc.hpp"
#include "histograms.hpp"

//-------------------------------------------------------//
//recovering old analysis script from ALERT. "filelist" is the list of files to open.
//make sure to change the # of files in the for loop below to analysis wanted files
//Reconstruction done from mod_pythia (working code)
//files used here are properly reconstructed -> working with proper magnetic fields (no offset)
//
// code working, hold for classoftware version (intermitent - seems fixed now)
// pions found. using SIDIS in Pi+ (only)
// proceeding to determine other variables useful for TMDs in RGD
// using 2 types of files: Sn_files and De_files in order to compare nuclei to nucleons. /!\
// using data filtering in variables. DC cuts ( dcfiducialcuts.hpp -working). CAL cuts (RE CHECK). others
// TMD comparison study (R, Dpt, etc) were removed, find them in recovery_posCLAS.cpp
// Here we are nnot recovering MC databanks, remains TBD, has to be done for headers & all
//REFERENCE: Version not READY, on development. /!\
//DATE: 29/03/2023
//next TBD: improving code, fixing sloppy conditions.   HEADERS TBD, find a way to determine also MC 
// ctrl+F on 'TBD' to see what's pending
//---------------------------------------------------------//


using namespace std;
int main() {
    VarHistograms myHists;
    MCHistograms MCHists;
    hist_electron elHists;
    hist_pip pipHists;
    hist_hadX hadHists;
    TCanvas* cc=new TCanvas("first","first"); //create new canvas
    cc->Divide(3,3);
    TCanvas* c2=new TCanvas("three","three");
    c2->Divide(2,3);
    TCanvas* c3=new TCanvas("third","third");
    c3->Divide(3,4);
    TCanvas* c4=new TCanvas("fourth","fourth");
    c4->Divide(4,3);
    TCanvas* cN=new TCanvas("nathan","nathan");
    cN->Divide(2,2);
    TCanvas* pres=new TCanvas("presentation","presentation");
    pres->Divide(3,3);

    double MC_theta_e, MC_Q2, MC_W, MC_y, MC_gamnu,MC_x_b, MC_t, MC_phih, MC_Epho; //for MC  ---[De]
    double theta_e, Q2, W, y,gamnu,x_b, t , phih , Epho,  rawphih, rawphih2, MCrawphih, MCrawphih2;   //for REC  ---[De]

    double MC_angle_norm, MC_sign, angle_norm, sign;   //---[De]

    double t_He, phih_He, angle_norm_He,sign_He;		//useless
    double MC_ang_phph, ang_phph;				//useless?
    double theta_el, phi_el, theta_pho, phi_pho, p_el, p_pho , angphi_He , angth_He, p_MC_He,p_He,theta_pip, phi_pip;
    double MCtheta_el, MCphi_el, MCtheta_pho, MCphi_pho, MCp_el, MCp_pho;
    double  phi_He , theta_He;					//useless
    double MCtheta_He, MCphi_He,  MCp_He ;			//useless
    double MCtheta_pip, MCphi_pip,  MCp_pip ;  //---[De]
    double MC_P_t, P_t, MC_Epi, Epi, MC_z, z, MC_P_trX, P_trX;	//variables 4 transverse momentum of pion and z, then transverse momentum for hadron X    //---[De]

    //cart coord for Pion+
    double pipX, pipY, pipZ;
    double scalX, scalY;
    //
    double MC_mmass, mmass;
    //4 missing plots//
    double missE, missM, missP, Delta_phi, angle_copl, p_T, p_x_He, p_y_He, MM2_pho, cone_ang;
    double MM_epiX,MM_epi,MM_eX,MM_piX, MC_MM_epiX, MC_MM_epi, MC_MM_eX, MC_MM_piX;
    double MissP_epiX, MissP_epi, MissP_eX, MissP_piX, MC_MissP_epiX, MC_MissP_epi, MC_MissP_eX, MC_MissP_piX;
    double phi_missX, theta_missX, MC_phi_missX, MCtheta_missX;
    //momenta for all parts
    double px_scal, py_scal, pz_scal, px_scapi, py_scapi, pz_scapi, px_scah, py_scah, pz_scah;
    double MCimp_pip, imp_pip, imp_e, MCimp_e;
    double deltaphih;
    double phih_De;		//-----value 4 phi in De, create outside of the De loop to keep it's existence in Sn loop (?)
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
    //===== Histogram creation =====//
	//test (DELETE THIS)//
    

    TH1F* h_scapiE = new TH1F("nrg(pi)", "nrg(pi)", nubin, 0.0, 10);
    TH1F* h_MCscapiE = new TH1F("MCnrg(pi)", "MCnrg(pi)", nubin, 0.0, 10);
    // REC histograms  //
        //see histograms.hpp
    //simulations histograms from REAL MC//
        // see histograms.hpp

    TH1F* h_MCmiss_Pepi = new TH1F("MCmiss_epi", "MCmiss_epi", nubin, 0,20);
    TH1F* h_MCmiss_PeX = new TH1F("MCmiss_eX", "MCmiss_eX", nubin, 0,20);
    TH1F* h_MCmiss_PpiX = new TH1F("MCmiss_piX", "MCmiss_piX", nubin, 0,20);
    TH1F* h_MCmiss_PepiX = new TH1F("MCmiss_epiX", "MCmiss_epiX", nubin, 0,20);

    TH1F* h_miss_Pepi = new TH1F("miss_epi", "miss_epi", nubin, 0,20);
    TH1F* h_miss_PeX = new TH1F("miss_eX", "miss_eX", nubin, 0,20);
    TH1F* h_miss_PpiX = new TH1F("miss_piX", "miss_piX", nubin, 0,20);
    TH1F* h_miss_PepiX = new TH1F("miss_epiX", "miss_epiX", nubin, 0,20);
    // BSA //
    TH1* h_Phphi  = new TH1F("Ph_phi",  "Ph_phi", phibin ,0,360);
    TH1* h_Nhphi  = new TH1F("Nh_phi",  "Nh_phi", phibin ,0,360);
    TH1* h_BSA  = new TH1F("BSA",  "BSA", phibin ,0,360);
    TH1* h_MC_Phphi  = new TH1F("MC_Ph_phi",  "MC_Ph_phi", phibin ,0,360);
    TH1* h_MC_Nhphi  = new TH1F("MC_Nh_phi",  "MC_Nh_phi", phibin ,0,360);
    TH1* h_MC_BSA  = new TH1F("MC_BSA",  "MC_BSA", phibin ,0,360);
    //PHI Analysis//
    TH1* h_cphi  = new TH1F("cos_phi",  "cos_phi", phibin ,-1,1);
    TH1* h_sphi  = new TH1F("sin_phi",  "sin_phi", phibin ,-1,1);
    TH1* h_c2phi  = new TH1F("cos_2phi",  "cos_2phi", phibin ,-1,1);
    TH1* h_s2phi  = new TH1F("sin_2phi",  "sin_2phi", phibin ,-1,1);
    TH1* h_MCcphi  = new TH1F("MCcos_phi",  "MCcos_phi", phibin ,-1,1);
    TH1* h_MCsphi  = new TH1F("MCsin_phi",  "MCsin_phi", phibin ,-1,1);
    TH1* h_MCc2phi  = new TH1F("MCcos_2phi",  "MCcos_2phi", phibin ,-1,1);
    TH1* h_MCs2phi  = new TH1F("MCsin_2phi",  "MCsin_2phi", phibin ,-1,1);
    TH2F* h_cs = new TH2F("cos-sin","cos-sin",phibin , -1,1,phibin,-1,1);
    TH2F* h_2cs = new TH2F("cos2-sin2","cos2-sin2",phibin , -1,1,phibin,-1,1);
    TH2F* h_cc = new TH2F("cos-cos","cos-cos",100 , -1,1,100,-1,1);
    TH2F* h_ss = new TH2F("sin-sin","sin-sin",100 , -1,1,100,-1,1);
    TH2F* h_c2c2 = new TH2F("cos2-cos2","cos2-cos2",100 , -1,1,100,-1,1);
    TH2F* h_s2s2 = new TH2F("sin2-sin2","sin2-sin2",100 , -1,1,100,-1,1);
    TH1* h_instant  = new TH1F("instant phih",  "instant phih", nubin ,-3.2,3.2);
    TH1* h_MCinstant  = new TH1F("instant MCphih",  "instant MCphih", nubin ,-3.2,3.2);
    TH2F* h_InstvInst = new TH2F("instant versus instant","instant versus instant",nubin , 0,3.2,nubin,0,3.2);
    TH2F* h_RawvRaw = new TH2F("Raw versus Raw","Raw versus Raw",nubin , -360,360,nubin,-360,360);
    TH2F* h_trueRaw = new TH2F("true raw (angle)","true raw (angle)",nubin , -2,2,nubin,-2,2);

    TH1*  h_DPhih = new TH1F("Delta Phih",  "Delta Phih", nubin,-20,20);
    TH2F* h_DphihVimp = new TH2F("Imp_pi Vs Imp_pi","Imp_pi Vs Imp_pi",nubin , 0,10,nubin,0,10);
    TH2F* h_impeVimpe = new TH2F("Dphih Vs Imp_e","Dphih Vs Imp_e",nubin , 0,360,nubin,0,10);
    TH1*  h_Dimpe = new TH1F("Delta impe",  "Delta impe", nubin,-10,10);
    TH1*  h_Dimppi = new TH1F("Delta imppi",  "Delta imppi", nubin,-10,10);
    TH1*  h_DPhie = new TH1F("Delta phie",  "Delta phie", nubin,-60,60);
    TH1*  h_DPhipi = new TH1F("Delta phipi",  "Delta phipi", nubin,-60,60);

    TH2F* h_cordphiel = new TH2F("phi(el)","phi(el)",nubin , -180,180,nubin,-180,180);
    TH2F* h_cordphipi = new TH2F("phi(pi)","phi(pi)",nubin , -180,180,nubin,-180,180);

    TH1*  h_blank = new TH1F("blank",  "blank", nubin,0,100);
    //  2Ds  //
    TH2F* DC_allpre = new TH2F("DCcutspre","DCcutspre", nubin , -3,3,1000,-3,3);
    TH2F* DC_allpos = new TH2F("DCcutspos","DCcutspos", nubin , -3,3,1000,-3,3);

    TH2F* DC_scalpre = new TH2F("DCscal","DCscal", nubin , -3,3,1000,-3,3);
    TH2F* DC_scalpos = new TH2F("DCscalpos","DCscalpos", nubin , -3,3,1000,-3,3);
	//other//
    TH2F* h_QvX = new TH2F("QvsX","QvsX", nubin , 0,.8,100,0,8);	//part1
        TH2F* h_QvX2 = new TH2F("QvsX2","QvsX2", nubin , 0,.5,100,0,8); //part2
    TH2F* h_zvsPt = new TH2F("zvsp","zvsp", nubin , 0,1,100,0,.8);
        TH2F* h_zvsPt2 = new TH2F("zvsp2","zvsp2", nubin , 0,1,100,0,.8);
 	//comparisons//
     TH2F* h_QvQ = new TH2F("QvsQ","QvsQ", nubin ,  0.0, 6,nubin, 0.0, 6);
     TH2F* h_xvx = new TH2F("xvsx","xvsx", nubin ,0.0, 0.4,nubin,0.0, 0.4);
     TH2F* h_yvy = new TH2F("yvsy","yvsy", nubin , 0.0, 1.0,nubin,0.0, 1.0);
     TH2F* h_VvV = new TH2F("VvsV","VvsV", nubin , 0.0, 11,nubin,0.0, 11);
     TH2F* h_WvW = new TH2F("WvsW","WvsW", nubin ,1.5, 5,nubin,1.5, 5);
     TH2F* h_tvt = new TH2F("tvst","tvst", nubin , -15,0,nubin,-15,0);	// NOT DONE YET    -PENDING
     TH2F* h_phivphi = new TH2F("phivsphi","phivsphi", 108 ,0,360,108,0,360);
     TH2F* h_zvz = new TH2F("zvsz","zvsz", nubin , 0,1,nubin, 0,1);
     TH2F* h_ptpt = new TH2F("PtvsPt","PtvsPt", nubin , 0,2,nubin, 0,2);

    


	//====Comparisons2======//
    TH1F* Dcphi = new TH1F ("Delta c phi", "Delta c phi", nubin, -2, 2);
    TH1F* Dsphi = new TH1F ("Delta s phi", "Delta s phi", nubin, -2, 2);
    TH1F* Dc2phi = new TH1F ("Delta c 2phi", "Delta c 2phi", nubin, -2, 2);
    TH1F* Ds2phi = new TH1F ("Delta s 2phi", "Delta s 2phi", nubin, -2, 2);

	//==Particle Distributions==//

    TH2F* thVP_e = new TH2F ("thVP_e", "thVP_e", nubin, 1.5, 7, nubin, 0,40);
    TH2F* thVP_pi = new TH2F ("thVP_pi", "thVP_pi", nubin, 1, 7, nubin, 0, 40);
    TH2F* QVx = new TH2F ("QVx", "QVx", nubin, 0, 0.5, nubin, 0,8);
    TH2F*  zVpt = new TH2F ("zVpt", "zVpt", nubin, 0, 1.6, nubin, 0, 0.8);
    TH2F* xVz = new TH2F ("xVz", "xVz", nubin, 0, 1, nubin, 0,1);
    TH2F* PtrVz = new TH2F ("PtrVz", "PtrVz", nubin, 0, 1, nubin, 0, 2);
    TH2F* QVz = new TH2F ("QVz", "QVz", nubin, 0, 0.8, nubin, 0,8);
    TH2F*  PhVx = new TH2F ("PhVx", "PhVx", nubin, 0, 1, nubin, 0, 8);

	//=========Sn HISTOGRAMS=======//
    TH1F* hist_Q_Sn = new TH1F("Q2_Sn", "Q2_Sn", nubin, QminX, QmaxX);
    TH1F* hist_v_Sn = new TH1F("gamnu_Sn", "gamnu_Sn", nubin, vminX, vmaxX);
    TH1F* hist_y_Sn = new TH1F("y_Sn", "y_Sn", nubin, 0.0, 1.0);
    TH1F* hist_x_Sn = new TH1F("x_bj_Sn", "x_bj_Sn", nubin, 0.0, 0.6);
    TH1F* hist_W_Sn = new TH1F("W_Sn", "W_Sn", nubin, 1.5, 5);
    TH1F* hist_phih_Sn = new TH1F("phih_Sn", "phih_Sn", 12, 0,360);
    TH1F* hist_t_Sn = new TH1F("t_Sn", "t_Sn", nubin, -15,0);
    TH1F* hist_P_t_Sn = new TH1F("P_t_Sn", "P_t_Sn", nubin, 0,2);
    TH1F* hist_Pt2_Sn = new TH1F("Pt2_Sn", "Pt2_Sn", nubin, 0,2);
    TH1F* hist_Pt2_De = new TH1F("Pt2_De", "Pt2_De", nubin, 0,2);
    TH1F* hist_z_Sn = new TH1F("z_Sn", "z_Sn", nubin, zminX,zmaxX);
    TH1F* hist_MX_Sn = new TH1F("MX_Sn", "MX_Sn", nubin, 0,20);

    TH1* h_cphi_Sn  = new TH1F("cos_phi_Sn",  "cos_phi_Sn", phibin ,-1,1);
    TH1* h_sphi_Sn  = new TH1F("sin_phi_Sn",  "sin_phi_Sn", phibin ,-1,1);
    TH1* h_c2phi_Sn  = new TH1F("cos_2phi_Sn",  "cos_2phi_Sn", phibin ,-1,1);
    TH1* h_s2phi_Sn  = new TH1F("sin_2phi_Sn",  "sin_2phi_Sn", phibin ,-1,1);
    TH1* h_MCcphi_Sn  = new TH1F("MCcos_phi_Sn",  "MCcos_phi_Sn", phibin ,-1,1);
    TH1* h_MCsphi_Sn  = new TH1F("MCsin_phi_Sn",  "MCsin_phi_Sn", phibin ,-1,1);
    TH1* h_MCc2phi_Sn  = new TH1F("MCcos_2phi_Sn",  "MCcos_2phi_Sn", phibin ,-1,1);
    TH1* h_MCs2phi_Sn  = new TH1F("MCsin_2phi_Sn",  "MCsin_2phi_Sn", phibin ,-1,1);

    TH2F*  compcAN = new TH2F ("compcAN", "compcAN", nubin, -1, 1, nubin, -1, 1);
    TH2F*  compsAN = new TH2F ("compsAN", "compsAN", nubin, -1, 1, nubin, -1, 1);
    TH2F*  compc2AN = new TH2F ("compc2AN", "compc2AN", nubin, -1, 1, nubin, -1, 1);
    TH2F*  comps2AN = new TH2F ("comps2AN", "comps2AN", nubin, -1, 1, nubin, -1, 1);
    //=====histograms for all electrons /!\ ======//
	//(no coincidence) (for multip ratio)
    TH1F* hist_Q_Sn_e = new TH1F("Q2_Sn_e", "Q2_Sn_e", nubin, QminX, QmaxX);
    TH1F* hist_Q_De_e = new TH1F("Q2_De_e", "Q2_De_e", nubin, QminX, QmaxX);
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

	//====Histos for <cosphi>/<cosphi>====//
    TH2F* h_cQ_Sn = new TH2F("cosQSn", "cosQSn", nubin, -1,1, nubin,QminX ,QmaxX);
    TH2F* h_cQ_De = new TH2F("cosQDe", "cosQDe", nubin, -1,1, nubin,QminX ,QmaxX);
    TH2F* h_cv_Sn = new TH2F("cosvSn", "cosvSn", nubin, -1,1, nubin,vminX ,vmaxX);
    TH2F* h_cv_De = new TH2F("cosvDe", "cosvDe", nubin, -1,1, nubin,vminX ,vmaxX);
    TH2F* h_cz_Sn = new TH2F("coszSn", "coszSn", nubin, -1,1, nubin,zminX ,zmaxX);
    TH2F* h_cz_De = new TH2F("coszDe", "coszDe", nubin, -1,1, nubin,zminX ,zmaxX);
    TH2F* h_cpt_Sn = new TH2F("cosptSn", "cosptSn", nubin, -1,1, nubin,PtminX ,PtmaxX);
    TH2F* h_cpt_De = new TH2F("cosptDe", "cosptDe", nubin, -1,1, nubin,PtminX ,PtmaxX);

	//====Histos for <Deltapt>/====//
    TH2F* h_ptQ_Sn = new TH2F("ptQSn", "ptQSn", nubin, PtminX,PtmaxX, nubin,QminX ,QmaxX);
    TH2F* h_ptQ_De = new TH2F("ptQDe", "ptQDe", nubin, PtminX,PtmaxX, nubin,QminX ,QmaxX);
    TH2F* h_ptv_Sn = new TH2F("ptvSn", "ptvSn", nubin, PtminX,PtmaxX, nubin,vminX ,vmaxX);
    TH2F* h_ptv_De = new TH2F("ptvDe", "ptvDe", nubin, PtminX,PtmaxX, nubin,vminX ,vmaxX);
    TH2F* h_ptz_Sn = new TH2F("ptzSn", "ptzSn", nubin, PtminX,PtmaxX, nubin,zminX ,zmaxX);
    TH2F* h_ptz_De = new TH2F("ptzDe", "ptzDe", nubin, PtminX,PtmaxX, nubin,zminX ,zmaxX);
	
    TH1F* w_ptQ_Sn = new TH1F("weightedptQSn", "weightedptQSn", nubin, QminX ,QmaxX);
    TH1F* w_ptQ_De = new TH1F("weightedptQDe", "weightedptQDe", nubin, QminX ,QmaxX);
    TH1F* w_ptv_Sn = new TH1F("weightedptvSn", "weightedptvSn", nubin, vminX ,vmaxX);
    TH1F* w_ptv_De = new TH1F("weightedptvDe", "weightedptvDe", nubin, vminX ,vmaxX);
    TH1F* w_ptz_Sn = new TH1F("weightedptzSn", "weightedptzSn", nubin, zminX ,zmaxX);
    TH1F* w_ptz_De = new TH1F("weightedptzDe", "weightedptzDe", nubin, zminX ,zmaxX);
    w_ptQ_Sn->Sumw2();
    w_ptQ_De->Sumw2();
    w_ptv_Sn->Sumw2();
    w_ptv_De->Sumw2();
    w_ptz_Sn->Sumw2();
    w_ptz_De->Sumw2();
    TH1F* w_cQ_Sn = new TH1F("weightedcQSn", "weightedcQSn", nubin, QminX ,QmaxX);
    TH1F* w_cQ_De = new TH1F("weightedcQDe", "weightedcQDe", nubin, QminX ,QmaxX);
    TH1F* w_cv_Sn = new TH1F("weightedcvSn", "weightedcvSn", nubin, vminX ,vmaxX);
    TH1F* w_cv_De = new TH1F("weightedcvDe", "weightedcvDe", nubin, vminX ,vmaxX);
    TH1F* w_cz_Sn = new TH1F("weightedczSn", "weightedczSn", nubin, zminX ,zmaxX);
    TH1F* w_cz_De = new TH1F("weightedczDe", "weightedczDe", nubin, zminX ,zmaxX);
    TH1F* w_cpt_Sn = new TH1F("weightedcptSn", "weightedcptSn", nubin, PtminX ,PtmaxX);
    TH1F* w_cpt_De = new TH1F("weightedcptDe", "weightedcptDe", nubin, PtminX ,PtmaxX);
    w_cQ_Sn->Sumw2();
    w_cQ_De->Sumw2();
    w_cv_Sn->Sumw2();
    w_cv_De->Sumw2();
    w_cz_Sn->Sumw2();
    w_cz_De->Sumw2();
    w_cpt_Sn->Sumw2();
    w_cpt_De->Sumw2();
    TH1F* diffphiSn = new TH1F("diffphiSn", "diffphiSn", nubin, -20 ,20);
    TH1F* diffphiDe = new TH1F("diffphiDe", "diffphiDe", nubin, -20 ,20);
     TH2F* precutQx1 = new TH2F("precutQxDe", "precutQxDe", 100, 0.0,0.6, 100,QminX ,QmaxX);
    TH2F* poscutQx1 = new TH2F("poscutQxDe", "poscutQxDe", 100, 0.0,0.6, 100,QminX ,QmaxX);
     TH2F* precutQx2 = new TH2F("precutQxSn", "precutQxSn", 100, 0.0,0.6, 100,QminX ,QmaxX);
    TH2F* poscutQx2 = new TH2F("poscutQxSn", "poscutQxSn", 100, 0.0,0.6, 100,QminX ,QmaxX);



    //QHistogram q_hist;



	//=== Utilities ==== //
    double E_temp;

    //== counters ==//
    int counter_el = 0;   //electron counter
    int counter_pho= 0;     //photon counter
    int counter_pi = 0;   //pion (+) counter
    int counter_pim = 0;  //pion (-) counter
    int counter_pit = 0;//total pion counter
    int counter_pro= 0;     //proton counter
    int counter_neu = 0;   //neutron counter
    int counter_ev = 0;      //event counter
    int j ;   //photon saver
    int k;
    int counter_MC_el  = 0; //Simulated electron counter
    int counter_MC_pho = 0; //etc
    int counter_MC_He = 0;
    int counter_MC_pi = 0;   //pion (+) counter
    int counter_MC_pim = 0;  //pion (-) counter
    int counter_MC_pi0 = 0;  //pion (0) counter
    int counter_MC_pit = 0;//total pion counter
    int counter_MC_pro= 0;     //proton counter
    int counter_MC_neu = 0;   //neutron counter
    double e_sup_pho = 0;
    double e_MC_sup_pho =0;
    int counter_plus = 0;
    int counter_MC_plus = 0;
    int counter_minus = 0;
    int counter_MC_minus = 0;
    int counter_Ph = 0;
    int counter_Nh = 0;
    int counter_He = 0;
    int counterMC_He = 0;
    int coinci_He = 0;
    int ahit;
    int counter2_He  = 0;
    int counter2_pho = 0;
    int counter_end  = 0;
    int counter_elpi = 0;
    int counter_scalE =0;
    int counter_pi_De = 0;
    int counter_pi_Sn = 0;
    int counter_el_De = 0;
    int counter_el_Sn = 0;
    int MC_counter_pi_De = 0;
    int MC_counter_pi_Sn = 0;
    int MC_counter_el_De = 0;
    int MC_counter_el_Sn = 0;
    int MC_coinci_De = 0;
    int MC_coinci_Sn = 0;
    int coinci_De = 0;
    int coinci_Sn = 0;
    int counter_Deall_e = 0;
    int counter_Snall_e = 0;
    //ID	(chi2, eventID for each particle)
    double IDel, IDpip;

    //CALORIMETER variables
    double lu, lv, lw, du, dv, dw, layer, lay0, lay1, lay2, lay3;

    //== bools ==/
    bool val_MC_e, val_MC_ph, val_MC_He, val_e, val_ph, val_He, val_pip, val_pim, val_neu, val_pro,val_pit, val_MC_pit, val_MC_pip, val_MC_pim, val_MC_neu, val_MC_pro; //resets condition for electrons, and allow to count photons

    //== values  ==//
    double m_e = 0.000511, Mn = 0.939565, Mpro = 0.938272, M_He = 3.727379, m_pi = 0.134976, M_D2 = 1.875;
	//ELECTRON MASS ADJUST /!\ (???) //it was at 0.000000511 (not correct i guess)
    
	//static const char* const filelist[] = {"../../mywork/REC_CLAS22/r_eD-01.hipo","../../mywork/REC_CLAS22/r_eSn-01.hipo"};

	static const char* const filelist[] =	{"../../files2read/r_eD-01.hipo","../../files2read/r_eD-02.hipo","../../files2read/r_eSn-01.hipo","../../files2read/r_eSn-02.hipo"};
		
	//static const char* const filelist[] = {"../../mywork/REC_CLAS22/merged.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo","../../mywork/REC_CLAS22/out_out_rgd_d_1.hipo"};
//{"../../mywork/REC_CLAS22/readertest.hipo","../../mywork/REC_CLAS22/r_N1.hipo",
//,"../../mywork/REC_CLAS22/r_P9.hipo","../../mywork/REC_CLAS22/r_P10.hipo"};
    //number of files:1
   int files = 4;
	//int files is the total nber of files. change here when files are added on filelist[] --> filelist[files]


    for (int filenbr = 0; filenbr<files ; filenbr++){
        //filelist[filenbr]; is a const char....most probably...
        hipo::reader reader;
	    cout<<(filenbr*100)/files<<"% complete"<<endl;
        //const char *filename = "../../path/to/file/here";
        reader.open(filelist[filenbr]);
        hipo::dictionary factory;
        reader.readDictionary(factory);
        hipo::bank REC_Particle_bank(factory.getSchema("REC::Particle"));  //call REC Bank
        hipo::bank MC_Particle_bank(factory.getSchema("MC::Particle"));       //call MC Bank
        //***************//
        hipo::bank AHDC_Particle_bank(factory.getSchema("AHDCRec::Track"));
        hipo::bank AHDCMC_Particle_bank(factory.getSchema("AHDCRec::MC"));
        hipo::bank RECTraj(factory.getSchema("REC::Traj"));
	    hipo::bank RECCal(factory.getSchema("REC::Calorimeter"));
        //**************//
        hipo::event event;
        while (reader.next()) {

            counter_ev += 1;
            e_sup_pho = 0; //reset the photon max energy at the beginning of every new event.
            e_MC_sup_pho = 0;
            j = 0;   //reset the photon saver
            k = 0;
            ahit = 0;
	        counter_scalE =0;		//resetting the electron energy counter at every event (we want to recover the most energetic one)//
            val_e=false;
            val_ph=false;
            val_He=false;    //reset boolean
            val_MC_e  = false;
            val_MC_ph = false;
            val_MC_He = false;
	        val_pip = false;
	        val_pim = false;
	        val_pit = false;
	        val_neu = false;
	        val_pro = false;
	        val_MC_pit = false;
	        val_MC_pip = false;
	        val_MC_pim = false;
	        val_MC_neu = false;
	        val_MC_pro = false;
            reader.read(event);
            event.getStructure(REC_Particle_bank);
            event.getStructure(MC_Particle_bank);
            event.getStructure(AHDC_Particle_bank);
            event.getStructure(AHDCMC_Particle_bank);
	        event.getStructure(RECTraj);
	        event.getStructure(RECCal);
            int numRECrows = REC_Particle_bank.getRows();
            int numMCRows  = MC_Particle_bank.getRows();
	        int numRECcal  = RECCal.getRows();
	        //cout<<numRECcal<<endl;
	        for (int c = 0; c < numRECcal; ++c) {
	 	        lu = RECCal.getFloat("lu", c);
		        lv = RECCal.getFloat("lv", c);
		        lw = RECCal.getFloat("lw", c);
		        du = RECCal.getFloat("du", c);
		        dv = RECCal.getFloat("dv", c);
		        dw = RECCal.getFloat("dw", c);
       		    //layer = RECCal.getFloat("layer",c);
	        }
	        int numTRAJrows= RECTraj.getRows();
	        if (numTRAJrows==-1 || 0){continue;}
            //-----Q2//
            int numALERTRows=AHDC_Particle_bank.getRows();
	        int numALERTMCRows=AHDCMC_Particle_bank.getRows();
            //==   Vectors ==//
            //Vectors are defined HERE in order to be reset after every iteration
            TLorentzVector *v_MC_scal = new TLorentzVector; //scattered lepton
            TLorentzVector *v_incil = new TLorentzVector;   //incident lepton
            v_incil->SetPxPyPzE(0.0, 0.0, 11.0, 11.0); //has no MC version, since it is a predefined vector, REC and MC versions are the same
            TVector3 p_incil = v_incil->Vect();  //3D
            TLorentzVector *v_inciHe = new TLorentzVector;   //incident Helium
            v_inciHe->SetPxPyPzE(0.0, 0.0, 0.0, M_D2);
	        TLorentzVector *v_nucleontarg = new TLorentzVector;   //incident NUCLEONNNN
            v_nucleontarg->SetPxPyPzE(0.0, 0.0, 0.0, Mn);
            TLorentzVector *v_MC_scaph = new TLorentzVector;
            TLorentzVector *v_MC_vipho = new TLorentzVector;
            TLorentzVector *v_MC_diff  = new TLorentzVector;
            TLorentzVector *v_MC_diffmass  = new TLorentzVector;
	        TLorentzVector *v_MC_dif4epi = new TLorentzVector;
	        TLorentzVector *v_MC_dif4piX = new TLorentzVector;
	        TLorentzVector *v_MC_dif4eX = new TLorentzVector;
	        TLorentzVector *v_MC_dif4epiX = new TLorentzVector;
	        TLorentzVector *v_dif4epi = new TLorentzVector;
	        TLorentzVector *v_dif4piX = new TLorentzVector;
	        TLorentzVector *v_dif4eX = new TLorentzVector;
	        TLorentzVector *v_dif4epiX = new TLorentzVector;
            TLorentzVector *v_diffmass  = new TLorentzVector;	//re using this vector for X(hadron)
	        TLorentzVector *v_diffmiss  = new TLorentzVector;
	        TLorentzVector *v_misspho  = new TLorentzVector;		//for cone angle calc.
            TLorentzVector *v_MC_scaHe  = new TLorentzVector;
            TLorentzVector *v_scal  = new TLorentzVector; //scattered lepton
            TLorentzVector *v_scaph = new TLorentzVector; //scattered REAL photon
            TLorentzVector *v_vipho = new TLorentzVector; //virtual photon
            TLorentzVector *v_diff  = new TLorentzVector;
            TLorentzVector *v_scaHe  = new TLorentzVector; //ALERT!!!!!!!!
            TLorentzVector *v_diffHe  = new TLorentzVector;
	        TLorentzVector *v_HELIUM  = new TLorentzVector;
		//new vectors for SIDIS - TMDs (REC and MCs)
	        TLorentzVector *v_scapip  = new TLorentzVector;
	        TLorentzVector *v_scapim  = new TLorentzVector;
	        TLorentzVector *v_scapit  = new TLorentzVector;
	        TLorentzVector *v_scaneu  = new TLorentzVector;
	        TLorentzVector *v_scapro  = new TLorentzVector;
	        TLorentzVector *v_MC_scapip  = new TLorentzVector;
	        TLorentzVector *v_MC_scapim  = new TLorentzVector;
	        TLorentzVector *v_MC_scapit  = new TLorentzVector;
	        TLorentzVector *v_MC_scaneu  = new TLorentzVector;
	        TLorentzVector *v_MC_scapro  = new TLorentzVector;
   	        //TLorentzVector *v_null  = new TLorentzVector;

            //======GET ALL MC EVENTS======//
            for (int i = 0; i < numMCRows; ++i) {
        		if (MC_Particle_bank.getInt("pid", i)==211 || -211){
                    //Get simulated PION+
		            val_MC_pit = true;		//affecting boolean for every single pion detected
			//if distinction is needed it can be done later
			//pions are distinguished (+-) in this code so they can be separated in analysis
                    counter_MC_pit += 1;
		            v_MC_scapit->SetPx(MC_Particle_bank.getFloat("px", i));
                    v_MC_scapit->SetPy(MC_Particle_bank.getFloat("py", i));
                    v_MC_scapit->SetPz(MC_Particle_bank.getFloat("pz", i));
                    TVector3 p_MC_scapit = v_MC_scapit->Vect();
		    //creating a 3vect to define the 4vector's 4th coordinate w/ root
                    v_MC_scapit->SetVectM(p_MC_scapit, m_pi);

                }
		        if (MC_Particle_bank.getInt("pid", i)==211){
                    //Get simulated PION+
		            val_MC_pip = true;
                    counter_MC_pi += 1;
		            if (filenbr <= 2) {MC_counter_pi_De += 1;}
		            if (filenbr > 2) {MC_counter_pi_Sn += 1;}
		            v_MC_scapip->SetPx(MC_Particle_bank.getFloat("px", i));
                    v_MC_scapip->SetPy(MC_Particle_bank.getFloat("py", i));
                    v_MC_scapip->SetPz(MC_Particle_bank.getFloat("pz", i));
                    TVector3 p_MC_scapip = v_MC_scapip->Vect();
		            //creating a 3vect to define the 4vector's 4th coordinate w/ root
                    v_MC_scapip->SetVectM(p_MC_scapip, m_pi);
                }
		        if (MC_Particle_bank.getInt("pid", i)==-211 ){
                    //Get the simulated PION -     //detect PION MINUS
                    val_MC_pim = true;
		            counter_MC_pim += 1;
		            v_MC_scapim->SetPx(MC_Particle_bank.getFloat("px", i));
                    v_MC_scapim->SetPy(MC_Particle_bank.getFloat("py", i));
                    v_MC_scapim->SetPz(MC_Particle_bank.getFloat("pz", i));
                    TVector3 p_MC_scapim = v_MC_scapim->Vect();
                    v_MC_scapim->SetVectM(p_MC_scapim, m_pi);
                }
		        /*if (MC_Particle_bank.getInt("pid", i)==111 ){
                    //Get the simulated PION 0     //detect PION 0

		        counter_MC_pi0 += 1;
                }*/
		        if (MC_Particle_bank.getInt("pid", i)==2212 ){
                    //Get the simulated PROTON
		            val_MC_pro = true;
                    counter_MC_pro += 1;
		            v_MC_scapro->SetPx(MC_Particle_bank.getFloat("px", i));
                    v_MC_scapro->SetPy(MC_Particle_bank.getFloat("py", i));
                    v_MC_scapro->SetPz(MC_Particle_bank.getFloat("pz", i));
                    TVector3 p_MC_scapro = v_MC_scapro->Vect();
                    v_MC_scapro->SetVectM(p_MC_scapro, Mpro);
                }
		        if (MC_Particle_bank.getInt("pid", i)==2112 ){
                    //Get the simulated NEUTRON      //
                    val_MC_neu = true;
		            counter_MC_neu += 1;
		            v_MC_scaneu->SetPx(MC_Particle_bank.getFloat("px", i));
                    v_MC_scaneu->SetPy(MC_Particle_bank.getFloat("py", i));
                    v_MC_scaneu->SetPz(MC_Particle_bank.getFloat("pz", i));
                    TVector3 p_MC_scaneu = v_MC_scaneu->Vect();
                    v_MC_scaneu->SetVectM(p_MC_scaneu, Mn);
                }
                if (MC_Particle_bank.getInt("pid", i)==11 && val_MC_e==false){
                    //Get simulated ELECTRON
                    val_MC_e = true;
                    counter_MC_el += 1;
		            if (filenbr <= 2) {MC_counter_el_De += 1;}
		            if (filenbr > 2) {MC_counter_el_Sn += 1;}
                    v_MC_scal->SetPx(MC_Particle_bank.getFloat("px", i));
                    v_MC_scal->SetPy(MC_Particle_bank.getFloat("py", i));
                    v_MC_scal->SetPz(MC_Particle_bank.getFloat("pz", i));
                    TVector3 p_MC_scal = v_MC_scal->Vect();
                    v_MC_scal->SetVectM(p_MC_scal, m_e);
                }
                if (MC_Particle_bank.getInt("pid", i)==22 && val_MC_ph==false){
                    //Get simulated PHOTON
                    val_MC_ph = true;
                    v_MC_scaph->SetPx(MC_Particle_bank.getFloat("px", i));
                    v_MC_scaph->SetPy(MC_Particle_bank.getFloat("py", i));
                    v_MC_scaph->SetPz(MC_Particle_bank.getFloat("pz", i));
                    v_MC_scaph->SetE(sqrt(v_MC_scaph->Px()*v_MC_scaph->Px() + v_MC_scaph->Py()*v_MC_scaph->Py() + v_MC_scaph->Pz()*v_MC_scaph->Pz() ));
                }
                if (MC_Particle_bank.getInt("pid", i)==1000020040 && val_MC_He==false){
                    //Get simulated HELIUM
                    val_MC_He = true;
		            counter_MC_He += 1;
                    v_MC_scaHe->SetPx(MC_Particle_bank.getFloat("px", i));
                    v_MC_scaHe->SetPy(MC_Particle_bank.getFloat("py", i));
                    v_MC_scaHe->SetPz(MC_Particle_bank.getFloat("pz", i));
                    TVector3 p_MC_scaHe = v_MC_scaHe->Vect();
                    v_MC_scaHe->SetVectM(p_MC_scaHe, M_He);
                }
            }
            //======GET REC EVENTS=======//*
	        for (int i = 0; i < numRECrows; ++i) {
                /*
                  //test vector function !!!!
		        if (REC_Particle_bank.getInt("pid", i)==211 || -211){
                    TLorentzVector myVector;
                myVector= fillVector2(myVector,m_pi,  i);
                }
                
                    //Get the REC PION      //detect ANY pion /!\
		            val_pit = true;
                    counter_pit += 1;
		            v_scapit->SetPx(REC_Particle_bank.getFloat("px", i));
                    v_scapit->SetPy(REC_Particle_bank.getFloat("py", i));
                    v_scapit->SetPz(REC_Particle_bank.getFloat("pz", i));
                    TVector3 p_scapit = v_scapit->Vect();
                    v_scapit->SetVectM(p_scapit,m_pi);
                }
                */
		        if (REC_Particle_bank.getInt("pid", i)==211 ){
                    //Get the REC PION      //detect PION
		            val_pip = true;
                    counter_pi += 1;
		            if (filenbr <= 2) {counter_pi_De += 1;}
		            if (filenbr > 2) {counter_pi_Sn += 1;}
		            IDpip = REC_Particle_bank.getFloat("chi2pid", i);
		            v_scapip->SetPx(REC_Particle_bank.getFloat("px", i));
                    v_scapip->SetPy(REC_Particle_bank.getFloat("py", i));
                    v_scapip->SetPz(REC_Particle_bank.getFloat("pz", i));
                    TVector3 p_scapip = v_scapip->Vect();
                    v_scapip->SetVectM(p_scapip,m_pi);
                }
		        if (REC_Particle_bank.getInt("pid", i)==-211 ){
                    //Get the REC PION      //detect PION MINUS
                    val_pim = true;
		            counter_pim += 1;
		            v_scapim->SetPx(REC_Particle_bank.getFloat("px", i));
                    v_scapim->SetPy(REC_Particle_bank.getFloat("py", i));
                    v_scapim->SetPz(REC_Particle_bank.getFloat("pz", i));
                    TVector3 p_scapim = v_scapim->Vect();
                    v_scapim->SetVectM(p_scapim,m_pi);
                }
		        if (REC_Particle_bank.getInt("pid", i)==2212 ){
                    //Get the REC PROTON
                    val_pro = true;
		            counter_pro += 1;
		            v_scapro->SetPx(REC_Particle_bank.getFloat("px", i));
                    v_scapro->SetPy(REC_Particle_bank.getFloat("py", i));
                    v_scapro->SetPz(REC_Particle_bank.getFloat("pz", i));
                    TVector3 p_scapro = v_scapro->Vect();
                    v_scapro->SetVectM(p_scapro,Mpro);
                }
		        if (REC_Particle_bank.getInt("pid", i)==2112 ){
                    //Get the REC NEUTRON      //
                    val_neu = true;
		            counter_neu += 1;
		            v_scaneu->SetPx(REC_Particle_bank.getFloat("px", i));
                    v_scaneu->SetPy(REC_Particle_bank.getFloat("py", i));
                    v_scaneu->SetPz(REC_Particle_bank.getFloat("pz", i));
                    TVector3 p_scaneu = v_scaneu->Vect();
                    v_scaneu->SetVectM(p_scaneu,Mn);
                }
                if (REC_Particle_bank.getInt("pid", i)==11 && val_e==false ){
                    //Get the REC electron      //detect electron
                    val_e = true;
                    counter_el += 1;
		            if (filenbr <= 2) {counter_el_De += 1;}
		            if (filenbr > 2) {counter_el_Sn += 1;}
        		    IDel = REC_Particle_bank.getFloat("chi2pid", i);
                    v_scal->SetPx(REC_Particle_bank.getFloat("px", i));
                    v_scal->SetPy(REC_Particle_bank.getFloat("py", i));
                    v_scal->SetPz(REC_Particle_bank.getFloat("pz", i));
                    TVector3 p_scal = v_scal->Vect();
                    v_scal->SetVectM(p_scal,m_e);

                }
            }
            if (val_e==true && val_pip==true){		//---IGNORE THIS
		    //condition with coincidence on general pion not working
		    //we r using pim, can be switched to pip
		    //probably due to many pions in 1 event ......... TBD
		    // MOST PROBABLY A USELESS CONDITION [HERE]
                counter_elpi += 1;

		    //--ALERT used to be here--//

            }
				/////////////////////
                // ---BEGIN CALC---//
				/////////////////////
                //get the vqlues on qny possible vector contained in the event
                //suppose: ward against empty vectors due to the condition on val_e =true for ex
            //theta_el = v_scal->Theta();	//coord_el
            theta_el = cordtheta(v_scal);   
            phi_el = v_scal->Phi();		//coord_el
		    theta_pip = v_scapip->Theta()  ;	//bof
		    phi_pip = v_scapip->Phi() ;	//bof
            *v_vipho = *v_incil - *v_scal;          //-----keep
            TVector3 p_scal = v_scal->Vect();
		    p_el = p_scal.Mag();
            theta_e= acos( p_incil.Dot(p_scal) / ( p_incil.Mag() * p_scal.Mag() ) ); //useful
            //Q2  = 4 * v_incil->E() * v_scal->E() * pow( sin(theta_e/2), 2);
		    Q2 = calcQ(v_incil,v_scal);
            gamnu= calcgamnu(v_incil, v_scal);
            y = calcy(v_incil,gamnu);
            x_b = calcx(Q2, Mn, gamnu);
		    z = calcz(v_scapip,gamnu);
            W = calcW(Q2, Mn, gamnu);
            TVector3 p_vipho= v_vipho->Vect();
		    TVector3 p_scapip = v_scapip->Vect();
			TVector3 n1 = p_vipho.Cross(p_incil);
			TVector3 n2 = p_vipho.Cross(p_scapip);
            TVector3 n3 = n1.Cross(n2);
            *v_diffmass= *v_incil + *v_nucleontarg -*v_scal - *v_scapip;	//this vector corresponds to the 	X-HADRON
		    TVector3 p_diffmass= v_diffmass->Vect();
            MM_epiX= calcmissm(v_incil, v_nucleontarg,v_scal,v_scapip);

		    //theta_missX = v_dif4epi->Theta()  ;
		    //phi_missX = v_dif4epi->Phi() ;
            mmass = v_diffmass->M2(); //SQUARED		//bof
		    
 		    //p_T = sqrt( pow(p_x_He,2) + pow(p_y_He,2)	) ;
            px_scal=v_scal->Px();
 		    py_scal=v_scal->Py();
		    pz_scal=v_scal->Pz();
		    px_scapi=v_scapip->Px();
 		    py_scapi=v_scapip->Py();
		    pz_scapi=v_scapip->Pz();
            *v_diff  = *v_diffmass - *v_inciHe; 	//- *v_vipho;
            TVector3 p_scaph= v_scaph->Vect();
		    p_pho = p_scaph.Mag();
            Epho = v_scaph->E();
            t = v_diff->M2()  ;
		    //cone_ang =acos(( p_misspho.Dot(p_scaph) ) / ( p_misspho.Mag()*p_scaph.Mag() )) *180/PI;//?
		    //cout<<"test";
		    //bool dccut = DC_cuts_electron(RECTraj, 0);				//NEEDED (TBD)
		    //cout<<"t="<<t<<endl;
            //TVector3 n2 = p_vipho.Cross(autre_MC_scaHe);
            
            P_t =(n2.Mag() )/ (p_vipho.Mag() );
		    //TVector3 crossX = p_vipho.Cross(p_diffmass);
            angle_norm=( n1.Dot(n2) ) / ( n1.Mag()*n2.Mag() );
            sign = n3.Dot(p_vipho);
            P_trX = calcP_trX(v_vipho,v_diffmass );
		    pipX = v_scapip->X();
            pipY = v_scapip->Y();
	        scalX = v_scal->X();
		    scalY = v_scal->Y();
			if (sign<0 ){
	            phih =acos(angle_norm) *180.0/PI+180;
                if (filenbr%2 == 0) {h_Phphi->Fill(phih);}	
                if (filenbr%2 != 0) {h_Nhphi->Fill(phih);}	
                h_instant->Fill(acos(angle_norm));
                //phih = phih +180;
                rawphih =acos(angle_norm)*180.0/PI;
                rawphih2 = angle_norm;
	        }
            if (sign>0 ){					
	            phih =acos(angle_norm) *-180.0/PI+180;
                if (filenbr%2 == 0) {h_Phphi->Fill(phih);}	
                if (filenbr%2 != 0) {h_Nhphi->Fill(phih);}	
   	            h_instant->Fill(-acos(angle_norm));
                //phih = phih +180;
                rawphih =-acos(angle_norm)*180.0/PI;
                rawphih2 = angle_norm;
            } 
	        //======= START ANALYSIS 4 DEUTERIUM =======//
            if (filenbr <=2 ) {
		        //==== CALCS and FILLINGS for MC (Deuterium) ===//
                    //DELETED, if needed check previous versions//

                ////////////////////////////////////////////////////
                TVector3 autre_MC_scaph = v_MC_scaph->Vect();
                TVector3 autre_MC_scaHe = v_MC_scaHe->Vect();
                ////////////////////////////////////////////////////
	            //==== CALCS and FILLINGS for REC ===//		---[De]
	            if (val_e==true){				
	                if (Q2>1.5){
		            //hist_Q_De_e->Fill(Q2);
		                if (gamnu<8.5 && gamnu>2.5){
		                    //hist_v_De_e->Fill(gamnu);
			                //if (z>0.3 && z<0.7){
		                    //hist_z_De_e->Fill(z);
				            if (W>2) {		//!!!!!!!!//
					            counter_Deall_e += 1;
					            hist_Q_De_e->Fill(Q2);
					            hist_v_De_e->Fill(gamnu);					
				            }
				            //}
			            }
		            }
                    if (val_pip==true  ){	//this condition involves the condition in val_e
                        coinci_De += 1;
                        //hist_MX->Fill(mmass);
                        /*
                        
                        */
                        DC_allpre->Fill(pipX,pipY);
		                DC_scalpre->Fill(scalX,scalY);
                        elHists.h_theta_el->Fill(theta_el*180/PI);
		                elHists.h_phi_el->Fill(phi_el*180/PI);
		                elHists.coor2D_el->Fill(phi_el*180/PI,theta_el*180/PI);
		                pipHists.h_theta_pip->Fill(theta_pip*180/PI);
		                pipHists.h_phi_pip->Fill(phi_pip*180/PI);
		                pipHists.coor2D_pip->Fill(phi_pip*180/PI,theta_pip*180/PI);
		                //coor2D_missX->Fill(phi_missX*180/PI,theta_missX*180/PI);
            		    h_QvX->Fill(x_b,Q2);
		                h_zvsPt->Fill(z,P_t);
		                elHists.hID_el->Fill(IDel);
		                precutQx1->Fill(x_b,Q2);
                        if (abs(IDel)<3){		//first condition
		                    pipHists.hID_pip->Fill(IDpip);
		                    if (abs (IDpip)<2){
		                        //hist_Q->Fill(Q2);			//-------------------Q---//
		                        if (Q2>1.5){
                                    //myHists.hist1->Fill(y);
                                    myHists.hist_y->Fill(y);
		                            if (y<0.85 && y>0.25){
	    	    	                    //hist_z->Fill(z);		//-------------------z---//
			                            if (z>0.3 && z<0.7){
			    	                        //if (dccut == true){	// /!\ implementing DC CUT!!!!		(compatibility TBD)
				                                if (theta_el*180/PI>6){	// /!\ implementing CUT on theta coordinate for electrons!!!!
            			        		        //hist_v->Fill(gamnu);	//-------------------v---//
	    	    		   	                        myHists.hist_x->Fill(x_b);
            	    		    	                myHists.hist_W->Fill(W);
					                                if (W>2) {
	    	    		    	    	                myHists.hist_phih->Fill(phih);
	    	    		    	    	                //hist_P_t->Fill(P_t);	//-------------------pt--//
	    	    		    	    	                myHists.hist_MX->Fill(mmass);
						                                h_cphi->Fill(cos(phih));
						                                if (mmass>1.2){				//To remove the proton. Don't forget to draw the ct in the plot
				    	    	                            elHists.hpx_el->Fill(px_scal);
				    	    	                            elHists.hpy_el->Fill(py_scal);
				    	    	                            elHists.hpz_el->Fill(pz_scal);
					    	                                myHists.hist_t->Fill(t);
				    	    	                            pipHists.hpx_pi->Fill(px_scapi);
				    	    	                            pipHists.hpy_pi->Fill(py_scapi);
				    	    	                            pipHists.hpz_pi->Fill(pz_scapi);
				    	    	                            hadHists.hpx_X->Fill(v_diffmass->Px());
				    	    	                            hadHists.hpy_X->Fill(v_diffmass->Py());
				    	    	                            hadHists.hpz_X->Fill(v_diffmass->Pz());
  				    	    	                            hadHists.h_P_trX->Fill(P_trX);
				    	    	                            DC_allpos->Fill(pipX,pipY); //post fiducial filling pions
				    	    	                            DC_scalpos->Fill(scalX,scalY);//post fiducial filling electrons
					    	                                elHists.h_scalE->Fill(v_scal->E());
                            						        h_scapiE->Fill(v_scapip->E());
						                                    if (v_scapip->E()<3 && v_scapip->E()>1){            //+++++++++++++++++++++++++++//
				    	    	                                //if (lv>90 && lw>90) {				//calorimeter cut on 9cm (?)
						    	                                h_zvsPt2->Fill(z,P_t);
					    	                                    //h_cphi->Fill(cos(phih));			//Cos & Sin analysis
					    	                                    h_sphi->Fill(sin(phih));
					    	                                    h_c2phi->Fill(cos(2*phih));
					    	                                    h_s2phi->Fill(sin(2*phih));
						                                        h_cs->Fill(cos(phih),sin(phih));
						                                        h_2cs->Fill(cos(2*phih),sin(2*phih));
					            	                            h_QvX2->Fill(x_b,Q2);
							                                    h_cQ_De->Fill(cos(phih),Q2);
							                                    h_cv_De->Fill(cos(phih),gamnu);
							                                    h_cz_De->Fill(cos(phih),z);
							                                    h_cpt_De->Fill(cos(phih),P_t);
							                                    h_ptQ_De->Fill(P_t*P_t,Q2);
							                                    h_ptv_De->Fill(P_t*P_t,gamnu);
							                                    h_ptz_De->Fill(P_t*P_t,z);
							                                    w_ptQ_De->Fill(Q2,P_t*P_t);
							                                    w_ptv_De->Fill(gamnu, P_t*P_t);
							                                    w_ptz_De->Fill(z, P_t*P_t);
							                                    myHists.hist_Q->Fill(Q2);		//-----------------------Q--///
							                                    myHists.hist_v->Fill(gamnu);		//-----------------------v--///
							                                    myHists.hist_z->Fill(z);		//-----------------------z--///
							                                    myHists.hist_P_t->Fill(P_t);		//----------------------pt--//
							                                    hist_Pt2_De->Fill(P_t*P_t);
    							                                w_cQ_De->Fill(Q2,cos(phih));
    							                                w_cv_De->Fill(gamnu,cos(phih));
							                                    w_cz_De->Fill(z,cos(phih));
    							                                w_cpt_De->Fill(P_t*P_t,cos(phih));
							                                    poscutQx1->Fill(x_b,Q2);
					            	                            //q condition was here for the MC comparisons DELETED
							                                }
                       			       	                    // }					//end if on MC vs REC
 			 	            	                        }
					                                }
				                                }
			   	                            //}			//DC cut (compatibility TBD)
			                            }
		                            }
		                        }
		                    }
		                }
                    }
		        }       //condition for true e 
            }
		    //================================================================================================================//
		    //================================================================================================================//
		    //===============================================   Sn    ========================================================//
		    //================================================================================================================//
		    //================================================================================================================//
	        if (filenbr >2) {
	            if (val_e==true){		//=========for All electrons=========//
	    		    //counter inside the cut condition
    	    	    theta_el = v_scal->Theta();	//coord_el
                    phi_el = v_scal->Phi();		//coord_el
		            theta_pip = v_scapip->Theta()  ;	//bof
		            phi_pip = v_scapip->Phi() ;	//bof
                    *v_vipho = *v_incil - *v_scal;
                    TVector3 p_scal = v_scal->Vect();
		            p_el = p_scal.Mag();
                    theta_e= acos( p_incil.Dot(p_scal) / ( p_incil.Mag() * p_scal.Mag() ) ); //useful
                    Q2  = 4 * v_incil->E() * v_scal->E() * pow( sin(theta_e/2), 2);
		            gamnu= v_incil->E() - v_scal->E();
		            z = (v_scapip->E())/gamnu;
		            hist_z_Sn_e->Fill(z);
		            TVector3 p_vipho= v_vipho->Vect();
		            TVector3 p_scapip = v_scapip->Vect();
		            TVector3 n2 = p_vipho.Cross(p_scapip);
		            P_t =(n2.Mag() )/ (p_vipho.Mag() );
		            W = sqrt(Mn*Mn + 2*Mn*gamnu - Q2);
	                if (Q2>1.5){
		                //hist_Q_Sn_e->Fill(Q2);     	        
		                if (gamnu<8.5 && gamnu>2.5){
			                //hist_v_Sn_e->Fill(gamnu);
			                //if (z>0.3 && z<0.7){
				            //hist_z_Sn_e->Fill(z);
				            if (W>2) {
					            counter_Snall_e += 1;	
					            hist_Q_Sn_e->Fill(Q2);
					            hist_v_Sn_e->Fill(gamnu);					
				            }
				            //}
			            }
		            }
	            }
      	        
                ////////////////////////////////////////////////////
                TVector3 autre_MC_scaph = v_MC_scaph->Vect();
                TVector3 autre_MC_scaHe = v_MC_scaHe->Vect();
                ////////////////////////////////////////////////////

    	        //==== CALCS and FILLINGS for REC ===//		---[Sn]
                if (val_e==true  && val_pip==true  ){	//&& val_ph==true && v_scaph->E()>2
		            //other conditions can be added to the if-bracket
	                coinci_Sn += 1;
		            precutQx2->Fill(x_b,Q2);
		            if (abs(IDel)<3){		//first    condition
		                if (abs (IDpip)<2){     //this both cqn be merged in one condition TBDD
				            //hist_Q_Sn->Fill(Q2);				///-----------------------Q---//
		                    if (Q2>1.5){
                                hist_y_Sn->Fill(y);
		                        if (y<0.85 && y>0.25){
	    	    	                //hist_z_Sn->Fill(z);				//---------------------z--//	
			                        if (z>0.3 && z<0.7){
				                        //if (dccut == true){	// /!\ implementing DC CUT!!!!
				                            if (theta_el*180/PI>6){	// /!\ implementing CUT on theta coordinate for electrons!!!!
                        					    //hist_v_Sn->Fill(gamnu);		//---------------------v--//
	    	    		                   	    hist_x_Sn->Fill(x_b);
            	    		    	            hist_W_Sn->Fill(W);
					                            if (W>2) {
	    	    		    	    	            hist_phih_Sn->Fill(phih);
	    	    		    	    	            //hist_P_t_Sn->Fill(P_t);		//-------------------pt---//
	    	    		    	    	            hist_MX_Sn->Fill(mmass);
						                            h_cphi_Sn->Fill(cos(phih));
						                            if (mmass>1.2){				//To remove the proton. Don't forget to draw the ct in the plot
						                                if (v_scapip->E()<3 && v_scapip->E()>1){            ///+++++++++++++++++++++//
				    	    	                            //if (lv>90 && lw>90) {				//calorimeter cut on 9cm (?)
						    	                            //h_cphi_Sn->Fill(cos(phih));			//Cos & Sin analysis
					    	                                h_sphi_Sn->Fill(sin(phih));
					    	                                h_c2phi_Sn->Fill(cos(2*phih));
					    	                                h_s2phi_Sn->Fill(sin(2*phih));
							                                //h_delete->Fill(cos(phih),Q2 );
							                                h_cQ_Sn->Fill(cos(phih),Q2);
							                                h_cv_Sn->Fill(cos(phih),gamnu);
							                                h_cz_Sn->Fill(cos(phih),z);
							                                h_cpt_Sn->Fill(cos(phih),P_t);
							                                h_ptQ_Sn->Fill(P_t*P_t,Q2);
							                                h_ptv_Sn->Fill(P_t*P_t,gamnu);
							                                h_ptz_Sn->Fill(P_t*P_t,z);
							                                w_ptQ_Sn->Fill(Q2,P_t*P_t);
							                                w_ptv_Sn->Fill(gamnu, P_t*P_t);
							                                w_ptz_Sn->Fill(z, P_t*P_t);
							                                hist_Q_Sn->Fill(Q2);		//---------------------Q--//	
							                                hist_v_Sn->Fill(gamnu);		//---------------------v--//									
							                                hist_z_Sn->Fill(z);		//---------------------z--//	
							                                hist_P_t_Sn->Fill(P_t);		//--------------------pt--//
							                                hist_Pt2_Sn->Fill(P_t*P_t);
							                                w_cQ_Sn->Fill(Q2,cos(phih));
    							                            w_cv_Sn->Fill(gamnu,cos(phih));
							                                w_cz_Sn->Fill(z,cos(phih));
    							                            w_cpt_Sn->Fill(P_t*P_t,cos(phih));
							                                poscutQx2->Fill(x_b,Q2);
                       			       	                }					
 			 	            	                    }
					                            }
				                            }   
			   	                        //}		//4 DC cut (compatibility TBD)
			                        }
		                        }
		                    }
		                }
		            }
                }
	        }
        }
    }


    



	cout<<"tot: "<<counter_elpi<<endl;
	cout<<"D : "<<coinci_De<<endl;
	cout<<"Sn: "<<coinci_Sn<<endl;

	cc->cd(1);
    myHists.hist_Q->SetLineColor(kBlue);
	hist_Q_Sn->SetLineColor(kRed);
    myHists.hist_Q->Draw(" hist");
	hist_Q_Sn->Draw(" hist,same");
    cc->Update();
    TLine *cutQ=new TLine(1.5,cc->cd(1)->GetUymin(),1.5,cc->cd(1)->GetUymax());
    cutQ->SetLineWidth(2);
    cutQ->SetLineStyle(kDashed);
    cutQ->Draw("");
    cc->cd(2);
    myHists.hist_y->SetLineColor(kBlue);
	hist_y_Sn->SetLineColor(kRed);
    myHists.hist_y->Draw("hist");
	hist_y_Sn->Draw("hist,same");
    cc->Update();
    TLine *cutY=new TLine(0.85,cc->cd(2)->GetUymin(),0.85,cc->cd(2)->GetUymax());
    cutY->SetLineWidth(2);
    cutY->SetLineStyle(kDashed);
    cutY->Draw("");

    cc->cd(3);
    myHists.hist_v->SetLineColor(kBlue);
	hist_v_Sn->SetLineColor(kRed);
    myHists.hist_v->Draw("hist");
	hist_v_Sn->Draw("hist,same");
    cc->Update();
    TLine *cutV=new TLine(8.5,cc->cd(3)->GetUymin(),8.5,cc->cd(3)->GetUymax());
    cutV->SetLineWidth(2);
    cutV->SetLineStyle(kDashed);
    cc->cd(4);
    myHists.hist_x->SetLineColor(kBlue);
	hist_x_Sn->SetLineColor(kRed);
    myHists.hist_x->Draw("hist");
	hist_x_Sn->Draw("hist,same");
    cc->cd(5);
    myHists.hist_W->SetLineColor(kBlue);
    myHists.hist_W->Draw("hist");
	hist_W_Sn->SetLineColor(kRed);
    hist_W_Sn->Draw("hist,same");
    cc->Update();
    TLine *cutW=new TLine(3.16,cc->cd(5)->GetUymin(),3.16,cc->cd(5)->GetUymax());
    cutW->SetLineWidth(2);
    cutW->SetLineStyle(kDashed);
    cutW->Draw("");
    cc->cd(8);
    myHists.hist_phih->SetLineColor(kBlue);
    myHists.hist_phih->Draw("hist");
	hist_phih_Sn->SetLineColor(kRed);
    hist_phih_Sn->Draw("hist,same");
     cc->cd(7);
    myHists.hist_z->SetLineColor(kBlue);
    myHists.hist_z->Draw("hist");
	hist_z_Sn->SetLineColor(kRed);
    hist_z_Sn->Draw("hist,same");
    cc->Update();
    TLine *cutz_min=new TLine(0.3,cc->cd(7)->GetUymin(),0.3,cc->cd(7)->GetUymax());
    cutz_min->SetLineWidth(2);
    cutz_min->SetLineStyle(kDashed);
    cutz_min->Draw("");
    TLine *cutz_max=new TLine(0.7,cc->cd(7)->GetUymin(),0.7,cc->cd(7)->GetUymax());
    cutz_max->SetLineWidth(2);
    cutz_max->SetLineStyle(kDashed);
    cutz_max->Draw("");
    cc->cd(6);						
    myHists.hist_P_t->SetLineColor(kBlue);
    myHists.hist_P_t->Draw("hist");
	hist_P_t_Sn->SetLineColor(kRed);
    hist_P_t_Sn->Draw("hist,same");
/*
    cc->cd(9);		
    hist_MX->SetLineColor(kBlue);
    hist_MX->Draw("hist");
	hist_MX_Sn->SetLineColor(kRed);
    hist_MX_Sn->Draw("hist,same");	//WEIRD, CHECK THIS TBD
    cc->Update();
    TLine *cutMMX=new TLine(1.2,cc->cd(9)->GetUymin(),1.2,cc->cd(9)->GetUymax());
    cutMMX->SetLineWidth(2);
    cutMMX->SetLineStyle(kDashed);
    cutMMX->Draw("");
*/
    cc->cd(9);
    //myHists.hist1->Draw("hist,same");

    cc->SaveAs("../correctionvar.pdf");
    cc->SaveAs("../correctionvar.root");

    TFile* file1 = new TFile("output1.root", "RECREATE");
    //previous line is a test 4 multipratio analysis in external cpp
    //analysis keeps in the following 
    myHists.hist_Q->Write();
    file1->Close();
    
}
