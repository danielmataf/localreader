#ifndef SRATIO_H
#define Sratio_H

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TGraphErrors.h>
#include <math.h>
#include <vector>
#include "Event.h"
#include "CutSet.h"

class sratio {

public: 
    sratio(CutSet cutsD, CutSet cutsA, const std::string& targetName);  //cutsA is for the Nucleus considered 
                                            //correction (TBD)in ratio needs to be propagated here too

    void FillHistograms(const Event& );
    void WriteHistos(const std::string );
    void calcSratio();

    void writeMatrixToFile(const std::string& );

    void DrawMonSinrat(const std::string&);
    
    void multiplotSratio(); //for only one target
    void multiplotSratio( sratio& , sratio&);   //for 3 targets (Sn, Cu, CxC)
    void multiplotSratio( sratio& , sratio& , sratio&);   //for 4 targets (Sn, Cu, C1, C2)

    void multiSratsimus(sratio&,sratio&,sratio&); //for 4 targets in simulation
    void multiSratRGD(sratio&,sratio&,sratio&); //for 4 targets in RGD data (real)

    void multiSratall(sratio&,sratio&,sratio&,sratio&,sratio&,sratio&,sratio&); //for 4 targets in sim and RGD (separation of C1 & C2)
    void multiSratall2(sratio&,sratio&,sratio&,sratio&,sratio&); // three ragets in sim and RGD ( no separation of CxC)


private:
    CutSet cutd;
    CutSet cuta;
    std::string targetName;
    //BINNING 4 MULTIBINNING
    int nubin = 100;
    int phibin= 10;
    int Cratiobin = 5   ;
    int Cratiobin_nu  = 5   ;
    int Cratiobin_z   = 5   ;
    int Cratiobin_Q = 5   ;
    int Cratiobin_x = 5   ;
    int Cratiobin_phih = 5   ;

    
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
    double xminCratio = 0.1;
    double xmaxCratio = 0.6;
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
    TH3F *hSratio_Q_nu_zD ; // = new TH3F("Cratio:nu,z,pt2,D", "Cratio:histo nu,z,pt2 for D",Cratiobin_Q, QminCratio, QmaxCratio, Cratiobin_nu,numinCratio,numaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );
    TH3F *hSratio_Q_nu_zA ; // = new TH3F("Cratio:nu,z,pt2,A", "Cratio:histo nu,z,pt2 for A",Cratiobin_Q, QminCratio, QmaxCratio, Cratiobin_nu,numinCratio,numaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );
    //useless now ? 
    TH3F *hSratio_nu_z_pt2A_onlye ; // = new TH3F("Cratio:Q2, nu,z,A onlye", "Cratio:histo_e Q2,nu,z for A", Cratiobin_Q, QminCratio, QmaxCratio,Cratiobin_nu,numinCratio,numaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );
    TH3F *hSratio_nu_z_pt2D_onlye ; // = new TH3F("Cratio:Q2, nu,z,D onlye", "Cratio:histo_e Q2,nu,z for A", Cratiobin_Q, QminCratio, QmaxCratio,Cratiobin_nu,numinCratio,numaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );
    TH1F *hSratio_nuA ; // = new TH1F("Cratio:nu_A", "Cratio:nu_A", Cratiobin_nu,numinCratio,numaxCratio) ;
    TH1F *hSratio_nuD ; // = new TH1F("Cratio:nu_D", "Cratio:nu_D", Cratiobin_nu,numinCratio,numaxCratio) ;




    //Weighted histograms for Cratio
    TH3F *h_wD_Sratio ; //= new TH3F("wD_Cratio", "wD_Cratio",Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,QminCratio,QmaxCratio,Cratiobin_z,zminCratio, zmaxCratio);
    TH3F *h_wA_Sratio ; //= new TH3F("wA_Cratio", "wA_Cratio",Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,QminCratio,QmaxCratio,Cratiobin_z,zminCratio, zmaxCratio);
    
    //count histograms for Cratio (3D)
    TH3F *h_D_Sratio3D ; //= new TH3F("countCratio:D", "count:wD_Cratio", Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,QminCratio,QmaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );
    TH3F *h_A_Sratio3D ; //= new TH3F("countCratio:A", "count:wD_Cratio", Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,QminCratio,QmaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );
    

    //histograms (3D, x,Q,z) for weighted  in squared cos  useful for VARIANCE
    TH3 *h_wD_sqSratio ; //= new TH3F("wD_sqCratio", "wD_sqCratio",Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,QminCratio,QmaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );//definition w/ 3 args
    TH3 *h_wA_sqSratio ; //= new TH3F("wA_sqCratio", "wA_sqCratio",Cratiobin_x  , xminCratio, xmaxCratio,Cratiobin_Q,QminCratio,QmaxCratio,Cratiobin_z,zminCratio, zmaxCratio  );


    //2D histograms for phi monitoring (useful?)
    TH2 *h_xphiD; //= new TH2F("xphiD", "xphiD", Cratiobin_x, xminCratio, xmaxCratio, Cratiobin_phih, phihminCratio, phihmaxCratio);
    TH2 *h_xphiA; //= new TH2F("xphiA", "xphiA", Cratiobin_x, xminCratio, xmaxCratio, Cratiobin_phih, phihminCratio, phihmaxCratio);
    TH2 *h_QphiD; //= new TH2F("QphiD", "QphiD", Cratiobin_Q, QminCratio, QmaxCratio, Cratiobin_phih, phihminCratio, phihmaxCratio);
    TH2 *h_QphiA; //= new TH2F("QphiA", "QphiA", Cratiobin_Q, QminCratio, QmaxCratio, Cratiobin_phih, phihminCratio, phihmaxCratio);
    TH2 *h_zphiD; //= new TH2F("zphiD", "zphiD", Cratiobin_z, zminCratio, zmaxCratio, Cratiobin_phih, phihminCratio, phihmaxCratio);
    TH2 *h_zphiA; //= new TH2F("zphiA", "zphiA", Cratiobin_z, zminCratio, zmaxCratio, Cratiobin_phih, phihminCratio, phihmaxCratio);
    TH2 *h_xQA; //= new TH2F("xQA", "xQA", Cratiobin_x, xminCratio, xmaxCratio, Cratiobin_Q, QminCratio, QmaxCratio);
    TH2 *h_xQD; //= new TH2F("xzA", "xzA", Cratiobin_x, xminCratio, xmaxCratio, Cratiobin_z, zminCratio, zmaxCratio);
    TH1 *h_phiMonA;
    TH1 *h_phiMonD;



    //Graphs
    TGraphErrors* graph_Srat= new TGraphErrors();

    //Storage of points (and errors)
    std::vector<std::vector<std::vector<double>>> SratioMatrix;    //three vectors for 3D matrix
    std::vector<std::vector<std::vector<double>>> errorSratioMatrix;

    //TFile* outputFile;
    CutSet cutsD;
    CutSet cutsSn;




};

#endif 
