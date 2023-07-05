
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

//REFERENCE: Version not READY yet, on development. /!\
//DATE: 15/06/2023 (git)
//Sn loop was deleted. histogramm header allows now to have general histograms  
// all drawing for Sn specific histograms were commented. unable now to compare two histograms from Sn and De in the same canvas, from this code
//if the comparison is required? One may call the Sn and D histograms from root files and compare them outside of main.cpp 
//fix list of fles so we can get the nbr of files w/o needing to manually insert the nbr TBD
//fix list of files so it can go either to a D list or a Sn list TBD
// 
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
    MultipRatio R_hists;
    TCanvas* cc=new TCanvas("first","first"); //create new canvas
    cc->Divide(3,3);

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
    // REC histograms  //
        //see histograms.hpp
    //simulations histograms from REAL MC//
        // see histograms.hpp


    //  2Ds  //
    TH2F* DC_allpre = new TH2F("DCcutspre","DCcutspre", nubin , -3,3,1000,-3,3);
    TH2F* DC_allpos = new TH2F("DCcutspos","DCcutspos", nubin , -3,3,1000,-3,3);

    TH2F* DC_scalpre = new TH2F("DCscal","DCscal", nubin , -3,3,1000,-3,3);
    TH2F* DC_scalpos = new TH2F("DCscalpos","DCscalpos", nubin , -3,3,1000,-3,3);
	//other//               //create a class for comparisons in 2D, or exploit the existing histograms for yada yada TBD
    TH2F* h_QvX = new TH2F("QvsX","QvsX", nubin , 0,.8,100,0,8);	//part1
    TH2F* h_QvX2 = new TH2F("QvsX2","QvsX2", nubin , 0,.5,100,0,8); //part2
    TH2F* h_zvsPt = new TH2F("zvsp","zvsp", nubin , 0,1,100,0,.8);
    TH2F* h_zvsPt2 = new TH2F("zvsp2","zvsp2", nubin , 0,1,100,0,.8);

	//=========Sn HISTOGRAMS=======// -----Replaced by  general histograms that can be used for any nuclear target (see header) 
    //=====histograms for all electrons /!\ ======//
	//(no coincidence) (for multip ratio)

	//=== Utilities ==== //
    double E_temp;

    //== counters ==//
    int counter_onlye =0;   //only electron counter, no coincidence
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


    int typeoffile; 
    cout<<"which type of file ? 1=D , 2 = Sn, 3=Cu"<<endl;
    cin>> typeoffile;




    //All files
	//static const char* const filelist[] =	{"../../files2read/r_eD-01.hipo","../../files2read/r_eD-02.hipo","../../files2read/r_eSn-01.hipo","../../files2read/r_eSn-02.hipo"};
	
    // D files  
    //static const char* const filelist[] =	{"../../files2read/r_eD-01.hipo","../../files2read/r_eD-02.hipo","../../files2read/r_eD-01.hipo","../../files2read/r_eD-02.hipo"};

    //Sn files 
    //static const char* const filelist[] =	{"../../files2read/r_eSn-01.hipo","../../files2read/r_eSn-02.hipo","../../files2read/r_eSn-01.hipo","../../files2read/r_eSn-02.hipo"};

    //number of files:4
   //const int files = 4;       //rechange to 4 if error in local
	//int files is the total nber of files. change here when files are added on filelist[] --> filelist[files]

    //imrpoving file reader: 
    const int files = 100;
    std::string basePath = "/volatile/clas12/dmat/gen/Deut/";       //path of files in ifarm
    //std::string basePath = "/volatile/clas12/dmat/gen/Sn/";
    std::string EndPath = "/sidis_mc-master/r_out_rgd10k_d.hipo";
    //std::string EndPath = "/sidis_mc-master/r_out_rgd10k_Sn.hipo";
    std::vector<std::string> filelist;
    for (int filenb = 1; filenb <= files; filenb++) {
        std::string filePath = basePath + std::to_string(filenb) + EndPath;
        filelist.push_back(filePath);
    }
    
    for (int filenbr = 0; filenbr<files ; filenbr++){
        //filelist[filenbr]; is a const char....most probably...
        hipo::reader reader;
	    cout<<(filenbr*100)/files<<"% complete"<<endl;
        //const char *filename = "../../path/to/file/here";
        reader.open(filelist[filenbr].c_str());
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
            				/////////////////////
                // ---BEGIN CALC---//
				/////////////////////
                //get the vqlues on qny possible vector contained in the event
                //suppose: ward against empty vectors due to the condition on val_e =true for ex
            //theta_el = v_scal->Theta();	//coord_el
            
			
	        //======= START ANALYSIS 4 DEUTERIUM =======//
            //======= START ANALYSIS 4 ANY TYPE OF FILE =======//
            if (filenbr <=2 ) {                                         
		        //==== CALCS and FILLINGS for MC (Deuterium) ===//
                    //add an option in cin at the beginning of the code in order to have an option for MC analysis.... TBD !
                    //DELETED, if needed check previous versions//

	            //==== CALCS and FILLINGS for REC ===//		---[De]
	            if (val_e==true){	
                    counter_onlye +=1;
                    theta_el = cordtheta(v_scal);   
                    phi_el = v_scal->Phi();		//coord_el
		            
                    *v_vipho = *v_incil - *v_scal;          //-----keep
                    TVector3 p_scal = v_scal->Vect();
		            p_el = p_scal.Mag();
                    theta_e= acos( p_incil.Dot(p_scal) / ( p_incil.Mag() * p_scal.Mag() ) ); //useful
		            Q2 = calcQ(v_incil,v_scal);
                    gamnu= calcgamnu(v_incil, v_scal);
                    y = calcy(v_incil,gamnu);
                    x_b = calcx(Q2, Mn, gamnu);
		            W = calcW(Q2, Mn, gamnu);
                    TVector3 p_vipho= v_vipho->Vect();
		            
			        TVector3 n1 = p_vipho.Cross(p_incil);
			        
                    px_scal=v_scal->Px();
 		            py_scal=v_scal->Py();
		            pz_scal=v_scal->Pz();
		            
                    TVector3 p_scaph= v_scaph->Vect();
		            p_pho = p_scaph.Mag();
                    Epho = v_scaph->E();
                    
		            //cone_ang =acos(( p_misspho.Dot(p_scaph) ) / ( p_misspho.Mag()*p_scaph.Mag() )) *180/PI;//?
		                //bool dccut = DC_cuts_electron(RECTraj, 0);				//NEEDED (TBD)
	                scalX = v_scal->X();
		            scalY = v_scal->Y();			
	                if (Q2>1.5){
		                R_hists.hist_Q_e->Fill(Q2);
		                if (gamnu<8.5 && gamnu>2.5){
                            R_hists.hist_v_e->Fill(gamnu);
			                //if (z>0.3 && z<0.7){
		                    //hist_z_De_e->Fill(z);
				            if (W>2) {		//!!!!!!!!//
					            counter_Deall_e += 1;
					            //hist_Q_De_e->Fill(Q2);
					            //hist_v_De_e->Fill(gamnu);					
				            }
				            //}
			            }
		            }
                    if (val_pip==true  ){	//this condition involves the condition in val_e
                        coinci_De += 1;
                        theta_pip = v_scapip->Theta()  ;	//bof
                        phi_pip = v_scapip->Phi() ;	//bof
                        z = calcz(v_scapip,gamnu);
                        t = v_diff->M2()  ;
                        
                        TVector3 p_scapip = v_scapip->Vect();
                        TVector3 n2 = p_vipho.Cross(p_scapip);
                        TVector3 n3 = n1.Cross(n2);
                        *v_diffmass= *v_incil + *v_nucleontarg -*v_scal - *v_scapip;	//this vector corresponds to the 	X-HADRON
		                TVector3 p_diffmass= v_diffmass->Vect();
                        MM_epiX= calcmissm(v_incil, v_nucleontarg,v_scal,v_scapip);
                        mmass = v_diffmass->M2(); //SQUARED		//bof
                        px_scapi=v_scapip->Px();
 		                py_scapi=v_scapip->Py();
		                pz_scapi=v_scapip->Pz();
                        P_t =(n2.Mag() )/ (p_vipho.Mag() );
		                TVector3 crossX = p_vipho.Cross(p_diffmass);
                        angle_norm=( n1.Dot(n2) ) / ( n1.Mag()*n2.Mag() );
                        sign = n3.Dot(p_vipho);
                        P_trX = calcP_trX(v_vipho,v_diffmass );
		                pipX = v_scapip->X();
                        pipY = v_scapip->Y();
                        
                        if (sign<0 ){
	                        phih =acos(angle_norm) *180.0/PI+180;
                            //if (filenbr%2 == 0) {h_Phphi->Fill(phih);}	//histogral filling for eventual BSA TBD
                            //if (filenbr%2 != 0) {h_Nhphi->Fill(phih);}	//keeping them, useful
                            //phih = phih +180;
                            rawphih =acos(angle_norm)*180.0/PI;
                            rawphih2 = angle_norm;
	                    }
                        if (sign>0 ){					
	                        phih =acos(angle_norm) *-180.0/PI+180;
                            //if (filenbr%2 == 0) {h_Phphi->Fill(phih);}	//TBD
                            //if (filenbr%2 != 0) {h_Nhphi->Fill(phih);}	//
                            //phih = phih +180;
                            rawphih =-acos(angle_norm)*180.0/PI;
                            rawphih2 = angle_norm;
                        } 
                        
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
		                //precutQx1->Fill(x_b,Q2);
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
						                                //h_cphi->Fill(cos(phih));
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
                            						        h_scapiE->Fill(v_scapip->E());      //replace this one with a class histogram in hpp TBD
						                                    if (v_scapip->E()<3 && v_scapip->E()>1){            //+++++++++++++++++++++++++++//
				    	    	                                //if (lv>90 && lw>90) {				//calorimeter cut on 9cm (?)
					    	                                    //h_cphi->Fill(cos(phih));			//Cos & Sin analysis
							                                    myHists.hist_Q->Fill(Q2);		//-----------------------Q--///
							                                    myHists.hist_v->Fill(gamnu);		//-----------------------v--///
							                                    myHists.hist_z->Fill(z);		//-----------------------z--///
							                                    myHists.hist_P_t->Fill(P_t);		//----------------------pt--//
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
		    //===============================================   Sn loop erased   =============================================//
        }
    }


    



	cout<<"tot: "<<counter_elpi<<endl;
	cout<<"D : "<<coinci_De<<endl;
	cout<<"Sn: "<<coinci_Sn<<endl;
    //the couts for Sn and D are useless now

	cc->cd(1);
    myHists.hist_Q->SetLineColor(kBlue);
	//hist_Q_Sn->SetLineColor(kRed);
    myHists.hist_Q->Draw(" hist");
	//hist_Q_Sn->Draw(" hist,same");
    cc->Update();
    TLine *cutQ=new TLine(1.5,cc->cd(1)->GetUymin(),1.5,cc->cd(1)->GetUymax());
    cutQ->SetLineWidth(2);
    cutQ->SetLineStyle(kDashed);
    cutQ->Draw("");
    cc->cd(2);
    myHists.hist_y->SetLineColor(kBlue);
	//hist_y_Sn->SetLineColor(kRed);
    myHists.hist_y->Draw("hist");
	//hist_y_Sn->Draw("hist,same");
    cc->Update();
    TLine *cutY=new TLine(0.85,cc->cd(2)->GetUymin(),0.85,cc->cd(2)->GetUymax());
    cutY->SetLineWidth(2);
    cutY->SetLineStyle(kDashed);
    cutY->Draw("");

    cc->cd(3);
    myHists.hist_v->SetLineColor(kBlue);
	//hist_v_Sn->SetLineColor(kRed);
    myHists.hist_v->Draw("hist");
	//hist_v_Sn->Draw("hist,same");
    cc->Update();
    TLine *cutV=new TLine(8.5,cc->cd(3)->GetUymin(),8.5,cc->cd(3)->GetUymax());
    cutV->SetLineWidth(2);
    cutV->SetLineStyle(kDashed);
    cc->cd(4);
    myHists.hist_x->SetLineColor(kBlue);
	//hist_x_Sn->SetLineColor(kRed);
    myHists.hist_x->Draw("hist");
	//hist_x_Sn->Draw("hist,same");
    cc->cd(5);
    myHists.hist_W->SetLineColor(kBlue);
    myHists.hist_W->Draw("hist");
	//hist_W_Sn->SetLineColor(kRed);
    //hist_W_Sn->Draw("hist,same");
    cc->Update();
    TLine *cutW=new TLine(3.16,cc->cd(5)->GetUymin(),3.16,cc->cd(5)->GetUymax());
    cutW->SetLineWidth(2);
    cutW->SetLineStyle(kDashed);
    cutW->Draw("");
    cc->cd(8);
    myHists.hist_phih->SetLineColor(kBlue);
    myHists.hist_phih->Draw("hist");
	//hist_phih_Sn->SetLineColor(kRed);
    //hist_phih_Sn->Draw("hist,same");
     cc->cd(7);
    myHists.hist_z->SetLineColor(kBlue);
    myHists.hist_z->Draw("hist");
	//hist_z_Sn->SetLineColor(kRed);
    //hist_z_Sn->Draw("hist,same");
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
	//hist_P_t_Sn->SetLineColor(kRed);
    //hist_P_t_Sn->Draw("hist,same");
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
    
    
    cc->cd(1);
    R_hists.hist_Q_e->Draw("hist");
    cc->cd(2);
    R_hists.hist_v_e->Draw("hist");
    cc->cd(3);
    R_hists.hist_z_e->Draw("hist");
    cc->cd(4);
    R_hists.hist_P_t_e->Draw("hist");
    cc->SaveAs("../testelectron.pdf");
    cc->SaveAs("../testelectron.root");

cout<<"testcounter e"<<counter_onlye<<endl;

/*
    TFile* file1 = new TFile("output1.root", "RECREATE");
    //previous line is a test 4 multipratio analysis in external cpp
    //analysis keeps in the following 
    myHists.hist_Q->Write();
    file1->Close();
*/


    std::map<int, std::string> fileNames = {
        {1, "../output_D.root"},
        {2, "../output_Sn.root"},
        {3, "../output_Cu.root"}
    };

if (fileNames.find(typeoffile) == fileNames.end()) {
    cout << "Invalid option" << endl;
    // Handle invalid option error or provide appropriate action
}
else {
    std::string filename = fileNames[typeoffile];
    TFile* file = new TFile(filename.c_str(), "RECREATE");

    myHists.hist_Q->Write();
    myHists.hist_v->Write();
    myHists.hist_y->Write();
    myHists.hist_x->Write();
    myHists.hist_W->Write();
    myHists.hist_phih->Write();
    myHists.hist_t->Write();
    myHists.hist_P_t->Write();
    myHists.hist_z->Write();
    myHists.hist_MX->Write();
    R_hists.hist_Q_e->Write();
    R_hists.hist_v_e->Write();
    //R_hists.hist_z_e->Write();
    //R_hists.hist_P_t_e->Write();
    TTree tree("treecounter","treecounter");
    tree.Branch("counter_onlye", &counter_onlye, "counter_onlye/I");
    tree.Fill();
    file->Write();

    file->Close();
}
    // we have one only loop that could be used to calc all variables for any nuclear target 

    
}
