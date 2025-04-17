#ifndef CRATIO_H
#define CRATIO_H

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TGraphErrors.h>
#include <math.h>
#include <vector>
#include "Event.h"
#include "CutSet.h"

class cratio {

public: 
    cratio(CutSet cutsD, CutSet cutsA, const std::string& targetName);  //cutsA is for the Nucleus considered 
                                            //correction (TBD)in ratio needs to be propagated here too
    void FillDebug(const Event& event);
    void WriteDebugHistos(const std::string& filename);
    void FillHistograms(const Event& );
    void WriteHistos(const std::string );
    void calcCratio();
    void calcCratioCarbon( cratio& ) ; //this fct for self carbon ratio. 
    void multiplotCratioBis();  // this 4 plots self ratio CC
    void saveCratioHistos();


    //void PlotCratio(const std::string );
    void writeMatrixToFile(const std::string& );
    //void calculateCratio(int , TH1F* , TH1F* , int , TH1F* , TH1F* , int ,  TGraphErrors* , TGraphErrors* ); 
    void multiplotCratio();
    void multiplotCratio( cratio& , cratio&);

    
    void multiCratsimus( cratio&,  cratio& , cratio& );
    void multiCrattrue(  cratio&,  cratio& , cratio& );
    void multiCratall(   cratio&,  cratio& , cratio& , cratio&, cratio& , cratio&, cratio& );
    void multiCratall2 ( cratio& , cratio&, cratio& , cratio&,cratio& ); 

    //histograms for selfRatio, for proof ig 

    TH3F* getHwA() ;
    TH3F* getHA3D() ;
    TH3* getHwA2() ;
    //Debugging Ratio inspired, and also used in Dpt
    void ValidateHistograms();
    void LogBinContent();

    //getters from Ratio for MAtrixes 
        std::vector<std::vector<std::vector<double>>> getRatMatrix() const{
        return CratioMatrix;
    }
    std::vector<std::vector<std::vector<double>>> getErrorMatrix() const{
        return errorCratioMatrix;
    }


    //following elements are used for self ratio in CxC target
    std::vector<std::vector<std::vector<double>>> getRatMatrixbis() const{
        return CratioMatrixbis;
    }
    std::vector<std::vector<std::vector<double>>> getErrorMatrixbis() const{
        return errorCratioMatrixbis;
    }


~cratio();         //maybe I should put the destructor in the .cpp 




private:
    CutSet cutd;
    CutSet cuta;
    std::string targetName;
    //BINNING 4 MULTIBINNING
    int nubin = 100;
    int phibin= 10;
    int Cratiobin = 6   ;
    int Cratiobin_nu  = 6   ;
    int Cratiobin_z   = 6   ;
    int Cratiobin_Q = 6   ;
    int Cratiobin_x = 6   ;
    int Cratiobin_phih = 6   ;

    
    int Qbin  ;    
    int vbin  ;    
    int xbin  ;    
    int zbin  ;
    int phihbin;   
    int pt2bin;    

    //counter electrons only for pt2 & z (for nu we use a new histo see below)
    int counter_elSn = 0;
    int counter_elLD2 = 0;

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


    double numinCratio = 4;
    double numaxCratio = 9;
    double xminCratio = 0;
    double xmaxCratio = 0.4;
    double QminCratio = Constants::RcutminQ;       //correct these  4 Dpt. To be changed in constants  TBD   
    double QmaxCratio = Constants::RcutmaxQ;       //correct these  4 Dpt. To be changed in constants.h
    double zminCratio = Constants::RcutminZ;       //correct these  4 Dpt. To be changed in constants.h
    double zmaxCratio = Constants::RcutmaxZ;       //correct these  4 Dpt. To be changed in constants.h    
    double pt2minCratio = Constants::RcutminPt2;       //correct these  4 Dpt. To be changed in constants.h
    double pt2maxCratio = Constants::RcutmaxPt2;       //correct these  4 Dpt. To be changed in constants.h
    double phihminCratio = Constants::phihmin_default;   
    double phihmaxCratio = Constants::phihmax_default;  



//CORRECT HISTOS HERE. DONT USE PT2
    //histos after passcuthadrons 
    TH3F *hCratio_Q_nu_zD ; // = new TH3F("Cratio:nu,z,pt2,D", "Cratio:histo nu,z,pt2 for D",Cratiobin_Q, QminCratio, QmaxCratio, Cratiobin_nu,numinCratio,numaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );
    TH3F *hCratio_Q_nu_zA ; // = new TH3F("Cratio:nu,z,pt2,A", "Cratio:histo nu,z,pt2 for A",Cratiobin_Q, QminCratio, QmaxCratio, Cratiobin_nu,numinCratio,numaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );
    //useless now ? 
    TH3F *hCratio_nu_z_pt2A_onlye ; // = new TH3F("Cratio:Q2, nu,z,A onlye", "Cratio:histo_e Q2,nu,z for A", Cratiobin_Q, QminCratio, QmaxCratio,Cratiobin_nu,numinCratio,numaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );
    TH3F *hCratio_nu_z_pt2D_onlye ; // = new TH3F("Cratio:Q2, nu,z,D onlye", "Cratio:histo_e Q2,nu,z for A", Cratiobin_Q, QminCratio, QmaxCratio,Cratiobin_nu,numinCratio,numaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );
    TH1F *hCratio_nuA ; // = new TH1F("Cratio:nu_A", "Cratio:nu_A", Cratiobin_nu,numinCratio,numaxCratio) ;
    TH1F *hCratio_nuD ; // = new TH1F("Cratio:nu_D", "Cratio:nu_D", Cratiobin_nu,numinCratio,numaxCratio) ;




    //Weighted histograms for Cratio
    TH3F *h_wD_Cratio ; //= new TH3F("wD_Cratio", "wD_Cratio",Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,QminCratio,QmaxCratio,Cratiobin_z,zminCratio, zmaxCratio);
    TH3F *h_wA_Cratio ; //= new TH3F("wA_Cratio", "wA_Cratio",Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,QminCratio,QmaxCratio,Cratiobin_z,zminCratio, zmaxCratio);
    
    //count histograms for Cratio (3D)
    TH3F *h_D_Cratio3D ; //= new TH3F("countCratio:D", "count:wD_Cratio", Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,QminCratio,QmaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );
    TH3F *h_A_Cratio3D ; //= new TH3F("countCratio:A", "count:wD_Cratio", Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,QminCratio,QmaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );
    

    //histograms (3D, x,Q,z) for weighted  in squared cos  useful for VARIANCE
    TH3 *h_wD_sqCratio ; //= new TH3F("wD_sqCratio", "wD_sqCratio",Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,QminCratio,QmaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );//definition w/ 3 args
    TH3 *h_wA_sqCratio ; //= new TH3F("wA_sqCratio", "wA_sqCratio",Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,QminCratio,QmaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );



    TH3F *h_wD_CratioRE; 
    TH3F *h_wD_sqCratioRE;
    TH3 *h_D_Cratio3DRE; 

    TH3F *h_wA_CratioRE; 
    TH3F *h_wA_sqCratioRE;
    TH3 *h_A_Cratio3DRE; 

    //Debugging
    TH1F* h_debug_xB;
    std::vector<TH1F*> h_debug_Q2;
    std::vector<std::vector<TH1F*>> h_debug_z;
    bool debugInitialized = false;


    //Graphs
    TGraphErrors* graph_crat= new TGraphErrors();

    //Storage of points (and errors)
    std::vector<std::vector<std::vector<double>>> CratioMatrix;    //three vectors for 3D matrix
    std::vector<std::vector<std::vector<double>>> errorCratioMatrix;
        //bis are for  self Carbon
    std::vector<std::vector<std::vector<double>>> CratioMatrixbis;    //three vectors for 3D matrix
    std::vector<std::vector<std::vector<double>>> errorCratioMatrixbis;

    //this is for intermediary steps 
    std::vector<std::vector<std::vector<double>>> CosPhiA_Matrix;
    std::vector<std::vector<std::vector<double>>> CosPhiD_Matrix;
    std::vector<std::vector<std::vector<double>>> errorCosPhiA_Matrix;
    std::vector<std::vector<std::vector<double>>> errorCosPhiD_Matrix;

    //TFile* outputFile;
    CutSet cutsD;
    CutSet cutsSn;




};

#endif 