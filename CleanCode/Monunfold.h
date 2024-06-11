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

    void Fill(const Event& );
    void FillHistograms(const Event& );
    void FillHistogramsNoCuts( const Event&);
    void FillHistogramswCuts(const Event& );
    void FillHistogramsNoCutsMC(const Event& ); 
    void Fill2DHistogramsRECMC(const Event& ); 
    void CompareHistograms(Monunfold& , Monunfold& , const std::string& ) ;

    void FillMomentumHistograms(const Event& );

    void DrawHistograms(const std::string);
    void DrawHistoRec(const std::string);   //this 2 in order to separate MC from REC and monitor separetely   
    void DrawHistoMC(const std::string);    //this 2 in order to separate MC from REC and monitor separetely       
    void DrawCompRECMC(const std::string);
    void DrawHistTrueandSIM(Monunfold& , const std::string ) ;  //this for plot comparisons in true"s REC and sim's REC 

    void DrawMomentumHistograms(const std::string);
    void DrawMomentumElectronHistograms(const std::string );
    void DrawMomentumHadronHistograms(const std::string );
    void DrawHistogramsNoCuts(const std::string);
    void DrawVertexHistograms(const std::string);




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
    int WmaxX = 10;
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





    TFile outputFile;
    //creating a root file 


};

#endif 
