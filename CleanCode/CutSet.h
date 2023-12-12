#ifndef CUTSET_H
#define CUTSET_H
#include <TFile.h>

#include "Event.h"  
#include <TH1F.h>
#include <TH2F.h>



class CutSet {
public:
    CutSet(); //initialize the cut values
    ~CutSet();

    void SetCutQ(double minQ, double maxQ);
    void SetCutY(double minY, double maxY);
    void SetCutV(double minV, double maxV);
    void SetCutW(double minW, double maxW);
    void SetCutZ(double minZ, double maxZ);
    void SetCutPt2(double minPt2, double maxPt2);

    void SetCutVx(double, double);
    void SetCutVy(double, double);
    void SetCutVz(double, double);

    void SetCutHadPID( int );


    
    //check if an event passes the cuts
    bool PassCutSnTarget (const Event& ) ;
    bool PassCutLD2Target (const Event& ) ;
    bool PassCutCuTarget (const Event& ) ;
    bool PassCutsElectrons(const Event& ) ;
    bool PassCutsHadrons(const Particle& ) ;
    bool PassCuts4R(const Event& , const Particle& );
    
    
    void Chop(const std::string );
    void DrawChop(const std::string );
    
    //bool PassCut(double Q2, double y, double v, double w, double z) ;
    //we 
private:
    //for cuts
    double cutQMin, cutQMax;
    double cutYMin, cutYMax;
    double cutVMin, cutVMax;
    double cutWMin, cutWMax;
    double cutZMin, cutZMax;
    double cutPt2Min, cutPt2Max;

    double cutVxMin,cutVxMax;
    double cutVyMin,cutVyMax;
    double cutVzMin,cutVzMax;

    //for target cuts (specific)
    double VzminSn, VzmaxSn;
    double VzminLD2, VzmaxLD2;
    double VzminCu, VzmaxCu;
    double VzminCxC, VzmaxCxC;

    double hadPID;

    //for histograms
    int QminX = 0;
    int QmaxX = 6;
    int nubin = 100;
    int phibin= 10;
    int Rbin = 10;
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

    //TH1F *h_Q2pre = new TH1F("Q2pre", "Q2pre", nubin, QminX , QmaxX);
    //TH1F *h_Q2pos = new TH1F("Q2pos", "Q2pos", nubin, QminX , QmaxX);
    //TH1F *h_xbpre= new TH1F("xbpre", "xbpre", nubin, xminX, xmaxX) ;
    //TH1F *h_ypre= new TH1F ("ypre" , "ypre" , nubin, yminX, ymaxX) ;
    //TH1F *h_nupre= new TH1F("nupre", "nupre", nubin,numinX,numaxX) ;
    //TH1F *h_W2pre= new TH1F("W2pre", "W2pre", nubin, WminX, 30) ;
    //TH1F *h_zpre= new TH1F("zpre", "zpre", nubin, zminX, zmaxX) ;
    //TH1F *h_pt2pre= new TH1F("pt2pre", "pt2pre", nubin, pt2minX, pt2maxX) ;
    //TH1F *h_phihpre= new TH1F("phihpre", "phihpre", 10, phihminX, phihmaxX) ;
    //TH1F *h_Vzpre= new TH1F("Vzpre", "Vzpre", nubin, -20, 30) ;
    //TH1F *h_Vxpre= new TH1F("Vxpre", "Vxpre", nubin, -5, 5) ;
    //TH1F *h_Vypre= new TH1F("Vypre", "Vypre", nubin, -5, 5) ;
//  
    //TH1F *h_Vzpos= new TH1F("Vzpos", "Vzpos", nubin, -20, 30) ;
    //TH2F *h_Qvx= new TH2F("Qvsx", "Qvsx", nubin, xminX, xmaxX, nubin, QminX , QmaxX);
//
    //TH1F *h_Q4R= new TH1F("Q2R", "A2R", 10, QminX , QmaxX) ;
    //TH1F *h_nu4R= new TH1F("nuR", "nuR", 10, numinX,numaxX) ;
    //TH1F *h_z4R= new TH1F("zR", "zR", 10, zminX, zmaxX) ;
    //TH1F *h_pt4R= new TH1F("pt2", "pt2R", 10, pt2minX, pt2maxX) ;


    //  TH1F *h_vertexZpre= new TH1F("targetVzpre", "vertex4targetpre", 100, -40, 40) ;
    //TH1F *h_vertexZpos= new TH1F("targetVzpre", "vertex4targetpre", 100, -40, 40) ;
    

};

#endif 