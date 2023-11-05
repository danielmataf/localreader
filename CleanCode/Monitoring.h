#ifndef MONITORING_H
#define MONITORING_H

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include "Event.h" // 
#include "CutSet.h"

class Monitoring {
public:
    Monitoring(CutSet);

    void FillHistograms(const Event& );
    void WriteHistogramsToFile(const std::string ); 
    void DrawHistograms(const std::string);
    //void FillQ2pre(const Event&  );
    //void Fillypre(const Event&  );
    //void Fillnupre(const Event&  );
    //void Fillnupre(const Event&  );
    //void FillW2pre(const Event&  );



private:
    CutSet cut1;
    
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
    int pt2maxX = 10;
    int phihminX = 0;
    int phihmaxX = 360;

    int counterel_R = 0;

    //create pointer and init them here
    TH1F *h_Q2 = new TH1F("Q2", "Q2", nubin, QminX , QmaxX);
    TH1F *h_xb= new TH1F("xb", "xb", nubin, xminX, xmaxX) ;
    TH1F *h_y= new TH1F ("y" , "y" , nubin, yminX, ymaxX) ;
    TH1F *h_nu= new TH1F("nu", "nu", nubin,numinX,numaxX) ;
    TH1F *h_W2= new TH1F("W2", "W2", nubin, WminX, 30) ;
    TH1F *h_z= new TH1F("z", "z", nubin, zminX, zmaxX) ;
    TH1F *h_pt2= new TH1F("pt2", "pt2", nubin, pt2minX, pt2maxX) ;
    TH1F *h_phih= new TH1F("phih", "phih", nubin, phihminX, phihmaxX) ;
    TH1F *h_vertexZ= new TH1F("targetVz", "vertex4target", 100, -40, 40) ;

    //add more histograms for other variables here
    
    //TH1F *h_Q2pre = new TH1F("Q2pre", "Q2pre", nubin, QminX , QmaxX);
    //TH1F *h_xbpre= new TH1F("xbpre", "xbpre", nubin, xminX, xmaxX) ;
    //TH1F *h_ypre= new TH1F ("ypre" , "ypre" , nubin, yminX, ymaxX) ;
    //TH1F *h_nupre= new TH1F("nupre", "nupre", nubin,numinX,numaxX) ;
    //TH1F *h_W2pre= new TH1F("W2pre", "W2pre", nubin, WminX, WmaxX) ;
    //TH1F *h_zpre= new TH1F("zpre", "zpre", nubin, zminX, zmaxX) ;
    //TH1F *h_pt2pre= new TH1F("pt2pre", "pt2pre", nubin, pt2minX, pt2maxX) ;
    //TH1F *h_phihpre= new TH1F("phihpre", "phihpre", nubin, phihminX, phihmaxX) ;




    TFile outputFile;
    //creating a root file 
};

#endif 