#ifndef DPT_H
#define DPT_H

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TGraphErrors.h>

#include <vector>
#include "Event.h"
#include "CutSet.h"

class deltaptsq {

public: 
    deltaptsq(CutSet , CutSet , const std::string& targetName);  //cutsA is for the Nucleus considered 
    
                                            //correction (TBD)in ratio needs to be propagated here too
    //~deltaptsq();          

    void FillHistograms(const Event& );
    void WriteHistos(const std::string );
    void calcDpt();
    void calcDptCarbon(deltaptsq& ) ;      
    void multiplotDptbis();

    void FillOnlyptandz(const Event& );
    void DrawOnlyptandz(const std::string& );


    //void PlotDpt(const std::string );
    void writeMatrixToFile(const std::string& );
    void calculateDpt(int , TH1F* , TH1F* , int , TH1F* , TH1F* , int ,  TGraphErrors* , TGraphErrors* ); 
    void multiplotDpt();
    void multiplotDpt( deltaptsq& , deltaptsq&);
    std::vector<std::vector<std::vector<double>>> getDptMatrix() const{
        return DptMatrix;
    }
    std::vector<std::vector<std::vector<double>>> getErrorMatrix() const{
        return errorDptMatrix;
    }


    void Dpttargetsimcomp(deltaptsq&); //this for comparison with simulationnnn
    void multiDptsimus( deltaptsq&,  deltaptsq& , deltaptsq& );
    void multiDpttrue(  deltaptsq&,  deltaptsq& , deltaptsq& );
    void multiDptall(   deltaptsq&,  deltaptsq& , deltaptsq& , deltaptsq&, deltaptsq& , deltaptsq&, deltaptsq& );
        //this fct plots every target in true and sil. C1 and C2 included
    void multiDptall2 ( deltaptsq& , deltaptsq&, deltaptsq& , deltaptsq&,deltaptsq& ); 
        //this same as b4 to plot only with the carbon on CxC


    TH3F* getHwA();
    TH3F* getHwA_sqpt2() ;
    TH3F* getHA3D() ;
    



private:
    CutSet cut1;
    CutSet cut2;
    std::string targetName;
    //BINNING 4 MULTIBINNING
    int nubin = 100;
    int phibin= 10;
    int Dptbin = 5   ;
    int Dptbin_nu  = 5   ;
    int Dptbin_z   = 5   ;
    int Dptbin_Q = 5   ;
    int Dptbin_x = 5   ;
    
    int Qbin  ;    
    int vbin  ;    
    int xbin  ;    
    int zbin  ;   
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



    double numinDpt = 4;
    double numaxDpt = 9;
    double xminDpt = 0.2;   
    double xmaxDpt = 0.4;
    double QminDpt = Constants::RcutminQ;       //correct these  4 Dpt. To be changed in constants  TBD   
    double QmaxDpt = Constants::RcutmaxQ;       //correct these  4 Dpt. To be changed in constants.h
    double zminDpt = Constants::RcutminZ;       //correct these  4 Dpt. To be changed in constants.h
    double zmaxDpt = Constants::RcutmaxZ;       //correct these  4 Dpt. To be changed in constants.h    
    double pt2minDpt = Constants::RcutminPt2;       //correct these  4 Dpt. To be changed in constants.h
    double pt2maxDpt = Constants::RcutmaxPt2;       //correct these  4 Dpt. To be changed in constants.h
    double phihminDpt = 0;
    double phihmaxDpt = 360;
    



//CORRECT HISTOS HERE. DONT USE PT2
    //histos after passcuthadrons 
    TH3F *hDpt_Q_nu_zD; // = new TH3F("Dpt:nu,z,pt2,D", "Dpt:histo nu,z,pt2 for D",Dptbin_Q, QminDpt, QmaxDpt, Dptbin_nu,numinDpt,numaxDpt,Dptbin_z,zminDpt, zmaxDpt  );
    TH3F *hDpt_Q_nu_zA; // = new TH3F("Dpt:nu,z,pt2,A", "Dpt:histo nu,z,pt2 for A",Dptbin_Q, QminDpt, QmaxDpt, Dptbin_nu,numinDpt,numaxDpt,Dptbin_z,zminDpt, zmaxDpt  );
    //useless now ? 
    TH3F *hDpt_nu_z_pt2A_onlye; // = new TH3F("Dpt:Q2, nu,z,A onlye", "Dpt:histo_e Q2,nu,z for A", Dptbin_Q, QminDpt, QmaxDpt,Dptbin_nu,numinDpt,numaxDpt,Dptbin_z,zminDpt, zmaxDpt  );
    TH3F *hDpt_nu_z_pt2D_onlye; // = new TH3F("Dpt:Q2, nu,z,D onlye", "Dpt:histo_e Q2,nu,z for A", Dptbin_Q, QminDpt, QmaxDpt,Dptbin_nu,numinDpt,numaxDpt,Dptbin_z,zminDpt, zmaxDpt  );
    TH1F *hDpt_nuA; // = new TH1F("Dpt:nu_A", "Dpt:nu_A", Dptbin_nu,numinDpt,numaxDpt) ;
    TH1F *hDpt_nuD; // = new TH1F("Dpt:nu_D", "Dpt:nu_D", Dptbin_nu,numinDpt,numaxDpt) ;




    //Weighted histograms for Delta pt
    TH3F *h_wD_pt; // = new TH3F("wD_pt", "wD_pt",Dptbin_x  , xminDpt, xmaxDpt,Dptbin_Q,QminDpt,QmaxDpt,Dptbin_z,zminDpt, zmaxDpt  );
    TH3F *h_wA_pt; // = new TH3F("wA_pt", "wA_pt",Dptbin_x  , xminDpt, xmaxDpt,Dptbin_Q,QminDpt,QmaxDpt,Dptbin_z,zminDpt, zmaxDpt  );
    
    //count histograms for Delta pt (3D)
    TH3F *h_D_pt3D; // = new TH3F("count:D", "count:wD_pt", Dptbin_x  , xminDpt, xmaxDpt,Dptbin_Q,QminDpt,QmaxDpt,Dptbin_z,zminDpt, zmaxDpt  );
    TH3F *h_A_pt3D; // = new TH3F("count:A", "count:wD_pt", Dptbin_x  , xminDpt, xmaxDpt,Dptbin_Q,QminDpt,QmaxDpt,Dptbin_z,zminDpt, zmaxDpt  );
    
    //histograms only to count pt squared values    USELESS  ???
    TH1F *h_D_onlypt2; // = new TH1F("count:D_pt2", "count:D_pt2", 5  , pt2minDpt, pt2maxDpt  );    //binning in 5 for coherence. Create a constant TBD
    TH1F *h_A_onlypt2; // = new TH1F("count:A_pt2", "count:A_pt2", 5  , pt2minDpt, pt2maxDpt  );    

    //histograms (3D, x,Q,z) for weighted  in squared pt2 (pt4) useful for VARIANCE
    TH3F *h_wD_sqpt2; // = new TH3F("wD_sqpt2", "wD_sqpt2",Dptbin_x  , xminDpt, xmaxDpt,Dptbin_Q,QminDpt,QmaxDpt,Dptbin_z,zminDpt, zmaxDpt  );//definition w/ 3 args
    TH3F *h_wA_sqpt2; // = new TH3F("wA_sqpt2", "wA_sqpt2",Dptbin_x  , xminDpt, xmaxDpt,Dptbin_Q,QminDpt,QmaxDpt,Dptbin_z,zminDpt, zmaxDpt  );

    TH1F *h_xb_A_had;
    TH1F *h_Q_A_had;
    TH1F *h_z_A_had;
    TH1F *h_z_D_had;
    TH1F *h_pt2_A_had;
    TH1F *h_pt2_D_had;




    //Graphs
    TGraphErrors* graph_rat= new TGraphErrors();

    std::vector<std::vector<std::vector<double>>> DptMatrix;    //three vectors for 3D matrix
    std::vector<std::vector<std::vector<double>>> errorDptMatrix;
    std::vector<std::vector<std::vector<double>>> DptMatrixbis;    //three vectors for 3D matrix
    std::vector<std::vector<std::vector<double>>> errorDptMatrixbis;

    //CutSet cutsD;
    //CutSet cutsSn;




};

#endif 