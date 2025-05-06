#ifndef RATIO_H
#define RATIO_H

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TGraphErrors.h>
#include "THnSparse.h"
#include <array>
#include <map>

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
    void FillDebug(const Event& ev, bool isDeuterium);
    void WriteDebugHistos(const std::string& filename);
    void DrawHistos(Ratio& );
    void calcR();   //has been pretty much unchanged for a while now 
    void calcRin5D();   //Dont forget to make an output in txt ig 
    void calcR_xB_Q2_z(); //RE variables for the ratio
    void saveRhistos();
    //void calcRcarbon() ;        //new
    void DrawSelfHistos(Ratio& ); //new

    void calcRcarbon(Ratio&) ;        //new
    void calcRwCC( Ratio& );
    void PlotRatio(const std::string );
    void writeMatrixToFile(const std::string& );
    void calculateMRat(int , TH1F* , TH1F* , int , TH1F* , TH1F* , int ,  TGraphErrors* , TGraphErrors* ); 
    void multiplotR();
    void multiplotR(Ratio& );
    void multiplotR( Ratio& , Ratio&);
    void multiplotR( Ratio& , Ratio& , Ratio& ,   Ratio& );
    void multiplotRbis();

    void multiRsimus( Ratio&, Ratio& , Ratio& );
    void multiRtrue( Ratio&, Ratio& , Ratio& );
    void multiRall( Ratio& , Ratio&, Ratio& , Ratio&,Ratio& , Ratio&, Ratio& );
        //this fct plots every target in true and sil. C1 and C2 included
    void multiRall2 ( Ratio& , Ratio&, Ratio& , Ratio&,Ratio& );
    void Rtargetsimcomp (Ratio& ); 
        //this same as b4 to plot only with the carbon on CxC
    TH1F* getHNuA() ;
    TH3F* getHNuzptA() ;


    //debugging? 
    void ValidateHistograms();
    void LogBinContent();
    
    std::vector<std::vector<std::vector<double>>> getRatMatrix() const{
        return ratMatrix;
    }
    std::vector<std::vector<std::vector<double>>> getErrorMatrix() const{
        return errorMatrix;
    }

    //following elements are used for self ratio in CxC target
    std::vector<std::vector<std::vector<double>>> getRatMatrixbis() const{
        return ratMatrixbis;
    }
    std::vector<std::vector<std::vector<double>>> getErrorMatrixbis() const{
        return errorMatrixbis;
    }



    ~Ratio();

private:
    //CutSet cut1;
    CutSet cutd;
    CutSet cuta;
    std::string targetName;
    //BINNING 4 MULTIBINNING
    int nubin = 100;
    int phibin= 10;
    int Rbin = Constants::Rbin_nu   ;
    int Rbin_nu  = 6   ;
    int Rbin_z   = 6   ;
    int Rbin_pt2 = 6   ;
    
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

    //histos to monitor self ration in CxC target
    TH1F *h_nu_A_had;
    TH1F *h_z_A_had;
    TH1F *h_pt2_A_had;

    TH1F *h_z_A;
    TH1F *h_z_D;
    TH1F *h_nu_D_had;
    TH1F *h_pt2_D_had;

    //histos for new analysis with new variables 
    TH2F *h_xB_Q2_D;
    TH3F *h_xB_Q2_z_D;
    TH2F *h_xB_Q2_A;
    TH3F *h_xB_Q2_z_A;





    //histos for 5D calc 
    private:
    //3D histos for electron count
    TH3F* h_3D_A_e;  // (Q2, xB, nu) 
    TH3F* h_3D_D_e;  // (Q2, xB, nu) 

    //5D histos for hadron count 
    THnSparseD* h_5D_A_had;  // (Q2, xB, nu, z, pt2) 
    THnSparseD* h_5D_D_had;  // (Q2, xB, nu, z, pt2) 
    //edges for THn
    static const int Rdim = 5;  //using 5 dimensions 
    double* binEdges[Rdim]; //need to make 5  arrays one per dimension

    // Optionally, store number of bins and ranges if reused
    int bins[Rdim] = {Constants::Rbin_nu, Constants::Rbin_nu, Constants::Rbin_nu, Constants::Rbin_nu, Constants::Rbin_nu}; //kinda useless but should work all bins are always 6 
    double binMins[Rdim] = {Constants::RcutminQ, Constants::Rcutminx, Constants::Rcutminnu, Constants::RcutminZ, Constants::RcutminPt2};
    double binMaxs[Rdim] = {Constants::RcutmaxQ, Constants::Rcutmaxx, Constants::Rcutmaxnu, Constants::RcutmaxZ, Constants::RcutmaxPt2};    //using cut lows and highs per variable (should work) but we may be exagerating some of the cuts. If issue come here for fixing
    //if need to edit cuts, look at tendencies on graphs and either replace here for values or directly in constants.h



    //Graphs
    TGraphErrors* graph_rat= new TGraphErrors();

    //storage of points (and errors)
    std::vector<std::vector<std::vector<double>>> ratMatrix;    //three vectors for 3D matrix
    std::vector<std::vector<std::vector<double>>> errorMatrix;
    
    std::vector<std::vector<std::vector<double>>> ratMatrixbis;    //three vectors for 3D matrix
    std::vector<std::vector<std::vector<double>>> errorMatrixbis;

    //new matrixes for variables x,Q,z
    std::vector<std::vector<std::vector<double>>> ratMatrix_xB;    //three vectors for 3D matrix
    std::vector<std::vector<std::vector<double>>> errorMatrix_xB;

    //storage for the 5D ratio points 
    //(Q2, xb, nu, z, pt2)
    //std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> ratMatrix5D;
    //std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> errorMatrix5D;
    //trying using arrays instead of vectors for fixed size matrixes. since they binning should be always the same eg: 5 or 6
    double ratMatrix5D[Constants::Rbin_nu][Constants::Rbin_nu][Constants::Rbin_nu][Constants::Rbin_nu][Constants::Rbin_nu] = {{{{{0.0}}}}};
    double errorMatrix5D[Constants::Rbin_nu][Constants::Rbin_nu][Constants::Rbin_nu][Constants::Rbin_nu][Constants::Rbin_nu] = {{{{{0.0}}}}};
    //this to aavoid resizing



//==new attempt at 5D, using two methods, 1st one uses THnsparse? Seconfd one uses array (map? )
    //implemetning arrays here
    //5D bin array: [Q2, xB, nu, pt2, z]  
    using Bin5D = std::array<int, 5>;   //uselss since all bins are the same but lets keep it for consistency
    //3D bin array: [Q2, xB, nu] 
    using Bin3D = std::array<int, 3>;   //same here
    //5D maps (hadron counts)
    std::map<Bin5D, double> hadronCountsD_5D;
    std::map<Bin5D, double> hadronCountsA_5D;
    
    //3D maps (electron normalization)
    std::map<Bin3D, double> electronCountsD_3D;
    std::map<Bin3D, double> electronCountsA_3D;
    //implementing THnSparse here
    //THnSparseD* h_5D_D = nullptr; //this is for all variables 
    //THnSparseD* h_5D_A = nullptr; 
    //TH3D* h_3D_D = nullptr; //this should be for electron counts: in q xb and nu
    //TH3D* h_3D_A = nullptr;



    CutSet cutsD;
    CutSet cutsSn;

};

#endif 
