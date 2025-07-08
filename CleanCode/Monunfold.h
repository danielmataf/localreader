#ifndef MONUNFOLD_H
#define MONUNFOLD_H

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include "Event.h" // 
#include "CutSet.h"
#include "constants.h"

//This class is made initially made ton avoid overppulating the class monitoring but it will implement the same function but mostly to fill True data and simulated REC and MC data. 
//We will be filling then histograms for this databanks in the usual variables, Q², x, y, nu, W², z, pt², phi_h, vertexZ
//we will additionally compare eventually variable values for different databanks in the same TH2F histogram. for example a 2D histo comparing Q² for MC and REC data.
//first priority is to plot a histogram that plots values of vz for only MC.

class Monunfold {
public:
    Monunfold(CutSet, const std::string& targetName); //adding name to constructor to add unique names to histograms
    //Monitoring(CutSet);
    ~Monunfold();          

    void WriteHistogramsToFile(const std::string ); 
    void FillOnlyVz(const Event& ) ;
    void DrawOnlyVz(Monunfold& , const std::string ) ;

    int SetResponseIndex(const Event&); //returns int index for constructing response matrix per evt.

    void Fill(const Event& );
    void FillHistograms(const Event& );
    void FillHistogramsNoCuts( const Event&);
    void FillHistogramswCuts(const Event& );
    void FillHistogramswCutsMC(const Event& );
    void FillHistogramsNoCutsMC(const Event& ); 
    void FillHistComp(const Event&, const Event&);
    void FillHistCompwCuts(const Event&, const Event&);

    void FillMomentumHistograms(const Event& );

    void DrawHistograms(const std::string filename);
    void DrawHistoRec(const std::string);   //this 2 in order to separate MC from REC and monitor separetely   
    void DrawHistoMC(const std::string);    //this 2 in order to separate MC from REC and monitor separetely       
    //void DrawCompRECMC(const std::string);
    void DrawCompRECMC(const std::string& filename); 
    void DrawHistTrueandSIM(Monunfold& , const std::string ) ;  //this for plot comparisons in true"s REC and sim's REC 

    void DrawMomentumHistograms(const std::string);
    void DrawMomentumElectronHistograms(const std::string );
    void DrawMomentumHadronHistograms(const std::string );
    void DrawHistogramsNoCuts(const std::string);
    void DrawVertexHistograms(const std::string);
    void SaveHistRoot(const std::string& ) ;
    void SaveHistMCRoot(const std::string& ) ;
    void DrawMomentainSim(const std::string& ) ;

    




private:
    CutSet cut1;
    std::string targetName;     //adding target name to add to hist names
    
    int QminX = 0;
    int QmaxX = 6;
    int nubin = 100;
    int phibin= 10;
    int Rbin = 5;
    int xminX = 0;
    int xmaxX = 0.4;
    int yminX = 0;
    int ymaxX = 1;
    int numinX = 0;
    int numaxX = 10;
    int WminX = 0;
    int WmaxX = 20;
    int zminX = 0;
    int zmaxX = 1;
    int pt2minX = 0;
    int pt2maxX = 3;
    int phihminX = 0;
    int phihmaxX = 360;

    int counterel_R = 0;

    //PROPAGATE BORDERS AND LIMITS TO CONSTANTS.H AND USE THEM HERE!!!!! urgent TBD 

    //borders 4 ratio
    //the values here should take the values of the implemented cuts
    int numinR = 4;
    int numaxR = 9;
    int zminR = 0.2;
    int zmaxR = 0.7;    
    int pt2minR = 0;
    int pt2maxR = 2;



    //create pointer and init them here
    TH1F *h_Q2; // = new TH1F(("Q2_" + targetName).c_str(), "Q2", nubin, QminX , QmaxX);
    TH1F *h_xb; //= new TH1F("xb", "xb", nubin, xminX, xmaxX) ;
    TH1F *h_y;  //= new TH1F ("y" , "y" , nubin, yminX, ymaxX) ;
    TH1F *h_nu; //= new TH1F("nu", "nu", nubin,numinX,numaxX) ;
    TH1F *h_W2; //= new TH1F("W2", "W2", nubin, WminX, 30) ;
    TH1F *h_z;  //= new TH1F("z", "z", nubin, zminX, zmaxX) ;
    TH1F *h_pt2;    //= new TH1F("pt2", "pt2", nubin, pt2minX, pt2maxX) ;
    TH1F *h_phih;   //= new TH1F("phih", "phih", nubin, phihminX, phihmaxX) ;
    TH1F *h_vertexZ;    //= new TH1F("targetVz", "vertex4target", 100, -40, 40) ;
    TH1F *h_vertexZ_pi;    //= new TH1F("targetVx", "vertex4target", 100, -40, 40) ;
    TH1F *h_vertexX;    //= new TH1F("targetVy", "vertex4target", 100, -40, 40) ;
    TH1F *h_vertexY;    //= new TH1F("targetVy", "vertex4target", 100, -40, 40) ;
    

    TH1F *h_thetaelectron;
    TH1F *h_rapport;
    TH1F *h_thetaelectronMC;
    TH1F *h_rapportMC;

    //Adding simulation MC histos for comparison (unfolding )
    TH1F *h_Q2MC; // = new TH1F(("Q2_" + targetName).c_str(), "Q2", nubin, QminX , QmaxX);
    TH1F *h_xbMC; //= new TH1F("xb", "xb", nubin, xminX, xmaxX) ;
    TH1F *h_yMC;  //= new TH1F ("y" , "y" , nubin, yminX, ymaxX) ;
    TH1F *h_nuMC; //= new TH1F("nu", "nu", nubin,numinX,numaxX) ;
    TH1F *h_W2MC; //= new TH1F("W2", "W2", nubin, WminX, 30) ;
    TH1F *h_zMC;  //= new TH1F("z", "z", nubin, zminX, zmaxX) ;
    TH1F *h_pt2MC;    //= new TH1F("pt2", "pt2", nubin, pt2minX, pt2maxX) ;
    TH1F *h_phihMC;   //= new TH1F("phih", "phih", nubin, phihminX, phihmaxX) ;
    TH1F *h_vertexZMC;    //= new TH1F("targetVz", "vertex4target", 100, -40, 40) ;
    TH1F *h_vertexXMC;    //= new TH1F("targetVx", "vertex4target", 100, -40, 40) ;
    TH1F *h_vertexYMC;    //= new TH1F("targetVy", "vertex4target", 100, -40, 40) ;

    
    TH2F *h_Q2comp;    //= new TH2F("xQ2", "xQ2", nubin, xminX, xmaxX, nubin, QminX, QmaxX) ;
    TH2F *h_xbcomp; //= new TH1F("xb", "xb", nubin, xminX, xmaxX) ;
    TH2F *h_ycomp;  //= new TH1F ("y" , "y" , nubin, yminX, ymaxX) ;
    TH2F *h_nucomp; //= new TH1F("nu", "nu", nubin,numinX,numaxX) ;
    TH2F *h_W2comp; //= new TH1F("W2", "W2", nubin, WminX, 30) ;
    TH2F *h_zcomp;  //= new TH1F("z", "z", nubin, zminX, zmaxX) ;
    TH2F *h_pt2comp;    //= new TH1F("pt2", "pt2", nubin, pt2minX, pt2maxX) ;
    TH2F *h_phihcomp;   //= new TH1F("phih", "phih", nubin, phihminX, phihmaxX) ;
    TH2F *h_vertexZcomp;    //= new TH1F("targetVz", "vertex4target", 100, -40, 40) ;
    TH2F *h_vertexZ_picomp;    //= new TH1F("targetVx", "vertex4target", 100, -40, 40) ;
    TH1F *h_DeltaVz;    //= new TH1F("targetVz", "vertex4target", 100, -40, 40) ;
    TH1F *h_pid;    //= new TH1F("pid", "pid", 100, -250, 250) ;
    TH2F *h_xQ2;    //= new TH2F("xQ2", "xQ2", nubin, xminX, xmaxX, nubin, QminX, QmaxX) ;
    TH2F *h_xQ2pos; //= new TH2F("xQ2pos", "xQ2pos", nubin, xminX, xmaxX, nubin, QminX, QmaxX) ;
    TH2F *h_pt2zMC;  //= new TH2F("pt2z", "pt2z", nubin, pt2minX, pt2maxX, nubin, zminX, zmaxX) ;

    TH1F *h_px_el;
    TH1F *h_py_el;
    TH1F *h_pz_el;
    TH1F *h_ptot_el;
    TH1F *h_px_pi;
    TH1F *h_py_pi;
    TH1F *h_pz_pi;
    TH1F *h_ptot_pi;
    TH1F *h_theta_el;
    TH1F *h_phi_el;
    TH1F *h_theta_pi;
    TH1F *h_phi_pi;


    TH1F *h_px_elMC;
    TH1F *h_py_elMC;
    TH1F *h_pz_elMC;
    TH1F *h_ptot_elMC;
    TH1F *h_px_piMC;
    TH1F *h_py_piMC;
    TH1F *h_pz_piMC;
    TH1F *h_ptot_piMC;
    TH1F *h_theta_elMC;
    TH1F *h_phi_elMC;
    TH1F *h_theta_piMC;
    TH1F *h_phi_piMC;


    TH1F *h_E_el;
    TH1F *h_E_pi;
    TH1F *h_E_elMC;
    TH1F *h_E_piMC;

    TH1F *h_deltaphipi;
    TH1F *h_deltathetapi;



    TH2F *h_px_elcomp;
    TH2F *h_py_elcomp;
    TH2F *h_pz_elcomp;
    TH2F *h_ptot_elcomp;
    TH2F *h_px_picomp;
    TH2F *h_py_picomp;
    TH2F *h_pz_picomp;
    TH2F *h_ptot_picomp;
    TH2F *h_theta_elcomp;
    TH2F *h_phi_elcomp;
    TH2F *h_theta_picomp;
    TH2F *h_phi_picomp;
    TH2F *h_E_elcomp;
    TH2F *h_E_picomp;
    TH2F *h_vxcomp;
    TH2F *h_vycomp;


    //Monitor generated and sim with Deltas dor Pi+ reconstruction 
    TH1F *h_px_piDelta;
    TH1F *h_py_piDelta;
    TH1F *h_pz_piDelta;
    TH1F *h_ptot_piDelta;
    TH1F *h_theta_piDelta;
    TH1F *h_phi_piDelta;
    TH1F *h_E_piDelta;

    TH1F *h_Delta_px_el;
    TH1F *h_Delta_py_el;
    TH1F *h_Delta_pz_el;
    TH1F *h_Delta_ptot_el;
    TH1F *h_Delta_theta_el;
    TH1F *h_Delta_phi_el;



    TH1F *h_chi2_el; // = new TH1F("chi2", "chi2", 100, 0, 10) ; //there is no chi2 in generated data obviously
    TH1F *h_chi2_pi; // = new TH1F("chi2", "chi2", 100, 0, 10) ; //there is no chi2 in generated data obviously



    TH1F *h_evtnbrdiff;   //evtnbr_sim - evtnbr_mc should be 0 iw we are reading the same event.    


    TFile outputFile;
    //creating a root file 


};

#endif 
