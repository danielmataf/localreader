#ifndef RATIO_H
#define RATIO_H

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TGraphErrors.h>

#include <vector>
#include "Event.h"
#include "CutSet.h"


class Ratio {
public:
    //Ratio(CutSet);
    Ratio(CutSet cutsD, CutSet cutsA, const std::string& targetName);  //cutsA is for the Nucleus considered 
                                        //propaate different nucleus to the eventreader mass as a n argument 
                                        // :!! TBD !
    
    
    //void FillHistograms1(const Event );
    
    //void WriteHistogramsToFile(const std::string ); 
    //void DrawHistograms(const std::string);
    void FillHistograms(const Event& );
    void WriteHistos(const std::string );
    void calcR();
    void calcRcarbon() ;        //new
    void PlotRatio(const std::string );
    void writeMatrixToFile(const std::string& );
    void calculateMRat(int , TH1F* , TH1F* , int , TH1F* , TH1F* , int ,  TGraphErrors* , TGraphErrors* ); 
    void multiplotR();
    void multiplotR( Ratio& , Ratio&);
    void multiplotR( Ratio& , Ratio& , Ratio& ,   Ratio& );
    
    std::vector<std::vector<std::vector<double>>> getRatMatrix() const{
        return ratMatrix;
    }
    std::vector<std::vector<std::vector<double>>> getErrorMatrix() const{
        return errorMatrix;
    }
    




private:
    //CutSet cut1;
    CutSet cutd;
    CutSet cuta;
    std::string targetName;
    //BINNING 4 MULTIBINNING
    int nubin = 100;
    int phibin= 10;
    int Rbin = 5   ;
    int Rbin_nu  = 5   ;
    int Rbin_z   = 5   ;
    int Rbin_pt2 = 5   ;
    
    int Qbin  ;  //=  
    int vbin  ;  //=  
    int xbin  ;  //=  
    int zbin  ;  //= 
    int pt2bin;  //=  

    //counter electrons only for pt2 & z (for nu we use a new histo see below)
    int counter_elSn = 0;
    int counter_elLD2 = 0;


    //Constants::RcutminQ
    //Constants::RcutmaxQ
    //Constants::RcutminY
    //Constants::RcutmaxY
    //Constants::RcutminW
    //Constants::RcutmaxW
    //
    //
    //
    //

    int QminX = 0;
    int QmaxX = 6;
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



    double numinR = 4;
    double numaxR = 9;
    double zminR = Constants::RcutminZ;
    double zmaxR = Constants::RcutmaxZ;    
    double pt2minR = Constants::RcutminPt2;
    double pt2maxR = Constants::RcutmaxPt2;

    //HISTOS FOR DEUTERIUM
    //create pointer and init them here
    //TH1F *h_Q2_D = new TH1F("Q2_D", "Q2_D", Rbin, QminX , QmaxX);
    //TH1F *h_xb_D= new TH1F("xb_D", "xb_D", Rbin, xminX, xmaxX) ;
    //TH1F *h_y_D= new TH1F ("y_D" , "y_D" , Rbin, yminX, ymaxX) ;
    //TH1F *h_nu_D= new TH1F("nu_D", "nu_D", Rbin,numinX,numaxX) ;
    //TH1F *h_nu_D_had= new TH1F("nu_D", "nu_D", Rbin,numinX,numaxX) ;
    //TH1F *h_W2_D= new TH1F("W2_D", "W2_D", Rbin, WminX, 30) ;
    //TH1F *h_z_D= new TH1F("z_D", "z_D", Rbin, zminX, zmaxX) ;
    //TH1F *h_pt2_D= new TH1F("pt2_D", "pt2_D", Rbin, pt2minX, pt2maxX) ;
    //TH1F *h_phih_D= new TH1F("phih_D", "phih_D", phibin, phihminX, phihmaxX) ;
    //add more histograms for other variables here

    //HISTOS FOR TARGET 
    //TH1F *h_xb_A= new TH1F("xb_A", "xb_A", Rbin, xminX, xmaxX) ;
    //TH1F *h_y_A= new TH1F ("y_A" ,"y_A" , Rbin, yminX, ymaxX) ;
    //TH1F *h_nu_A= new TH1F("nu_A", "nu_A", Rbin,numinX,numaxX) ;
    //TH1F *h_nu_A_had= new TH1F("nu_A", "nu_A", Rbin,numinX,numaxX) ;
    //TH1F *h_W2_A= new TH1F("W2_A", "W2_A", Rbin, WminX, 30) ;
    //TH1F *h_z_A= new TH1F("z_A", "z_A", Rbin, zminX, zmaxX) ;
    //TH1F *h_pt2_A= new TH1F("pt2_A", "pt2_A", Rbin, pt2minX, pt2maxX) ;
    //TH1F *h_phih_A= new TH1F("phih_A", "phih_A", phibin, phihminX, phihmaxX) ;

    //histos after passcuthadrons 
    TH3F *h_nu_z_pt2D; // = new TH3F("nu,z,pt2,D", "histo nu,z,pt2 for D", Rbin_nu,numinR,numaxR,Rbin_z,zminR, zmaxR,Rbin_pt2, pt2minR, pt2maxR  );
    TH3F *h_nu_z_pt2A; // = new TH3F("nu,z,pt2,A", "histo nu,z,pt2 for A", Rbin_nu,numinR,numaxR,Rbin_z,zminR, zmaxR,Rbin_pt2, pt2minR, pt2maxR  );
    //histo after passcutelectrons only e for nu unse only
    TH3F *h_nu_z_pt2A_onlye; // = new TH3F("nu,z,pt2,A onlye", "histo_e nu,z,pt2 for A", Rbin_nu,numinR,numaxR,Rbin_z,zminR, zmaxR,Rbin_pt2, pt2minR, pt2maxR  );
    TH3F *h_nu_z_pt2D_onlye; // = new TH3F("nu,z,pt2,D onlye", "histo_e nu,z,pt2 for A", Rbin_nu,numinR,numaxR,Rbin_z,zminR, zmaxR,Rbin_pt2, pt2minR, pt2maxR  );
    TH1F *h_nuA; // = new TH1F("nu_A", "nu_A", Rbin_nu,numinR,numaxR) ;
    TH1F *h_nuD; // = new TH1F("nu_D", "nu_D", Rbin_nu,numinR,numaxR) ;

    TH1F *h_nuC1; 
    TH1F *h_nuC2; 
    TH3F *h_nu_z_pt2C1; 
    TH3F *h_nu_z_pt2C2; 


    //Graphs
    TGraphErrors* graph_rat= new TGraphErrors();

    //Storage of points (and errors)
    std::vector<std::vector<std::vector<double>>> ratMatrix;    //three vectors for 3D matrix
    std::vector<std::vector<std::vector<double>>> errorMatrix;

    //TFile* outputFile;
    CutSet cutsD;
    CutSet cutsSn;

    //TFile outputFile;
    //creating a root file 
};

#endif 
