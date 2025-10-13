 #ifndef MONITORING_H
#define MONITORING_H

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include "Event.h" // 
#include "CutSet.h"
#include "constants.h"


class Monitoring {
public:
    Monitoring(CutSet, const std::string& targetName); //adding name to constructor to add unique names to histograms
    //Monitoring(CutSet);
    ~Monitoring();          

    void WriteHistogramsToFile(const std::string ); 
    void FillHistograms(const Event& );
    void Fill_sR_Histograms(const Event& event) ;
    void Fill_R_Histograms(const Event& event) ;
    void FillHistogramsNoCuts( const Event&);
    void FillHistogramswCuts(const Event& );
    void FillHistogramsNoCutsMC(const Event& ); 
    void FillResolutionFrom(Monitoring& other);

    void CheckLargeBins(const Event& ); //check large bins in order to fill them with the correct values imported from monunfold that checks simus
    //I just need to properly check actual data, then modifiy it
    void PrintRegionCounters(); //print large all bin event number in all regions
    void CheckFewBins(const Event& ); //check large bins in order to fill them with the correct values imported from monunfold that checks simus
    void PrintFewRegionCounters(); //print large all bin event number in all regions
    void CalcElectronRatio(Monitoring& monDeu);
    void CalcLARGEElectronRatio(Monitoring& monDeu);

    void FillPrePostCutHistograms(const Event& event);


    void initThetaBinning();            // called once in ctor
    void thetabinning(const Event& ev); // call from main event loop


    void FillMomentumHistograms(const Event& );

    void SaveKeyHistograms();
    void DrawHistograms(const std::string);
    void DrawHistTrueandSIM(Monitoring& , const std::string ) ;
    //void DrawHistTrueRECandMC(const std::string);
    void DrawR_Histograms(const std::string filename);
    void DrawHistogramsPos(const std::string ,const std::string );
    void DrawQmonitoring(Monitoring& , const std::string);

    void DrawMomentumHistograms(const std::string);
    void DrawMomentumElectronHistograms(const std::string );
    void DrawMomentumHadronHistograms(const std::string );
    void DrawEnergyHistograms(const std::string );
    void DrawCaloHistograms(const std::string );
    void DrawCherenkovHistograms(const std::string ) ;
    void DrawHistogramsNoCuts(const std::string);
    void DrawVertexHistograms(const std::string);
    void DrawHelicityHistograms(const std::string);
    void DrawCompRECMC(const std::string);
    void DrawMonSratio(const std::string);
    
    void SaveFINRoot(const std::string& ) ;

    void SaveHistRoot(const std::string& ) ;


    void FillDISforUnfoldDAT(const Event&);
    void saveDISforUnfoldRoot(const std::string& ) ;

    //void FillQ2pre(const Event&  );
    //void Fillypre(const Event&  );
    //void Fillnupre(const Event&  );
    //void Fillnupre(const Event&  );
    //void FillW2pre(const Event&  );



private:
    CutSet cut1;
    std::string targetName;     //adding target name to add to hist names
    
    int QminX = 0;
    int QmaxX = 6;
    int nubin = 100;
    int phibin= 10;
    int Rbin = 5;
    int xminX = 0;
    int xmaxX = 0.65;
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
    int counterprecuts = 0;
    int counterpassVz = 0;
    int counterpassQ2 = 0;
    int counterpassy = 0;
    int counterpassnu = 0;
    int counterpassw = 0;
    int counterpassz = 0;
    int counterpasspt2 = 0;
    int counterpasshadid = 0;


    //PROPAGATE BORDERS AND LIMITS TO CONSTANTS.H AND USE THEM HERE!!!!! urgent TBD 

    //borders 4 ratio
    //the values here should take the values of the implemented cuts
    int numinR = 4;
    int numaxR = 9;
    int zminR = 0.2;
    int zmaxR = 0.7;    
    int pt2minR = 0;
    int pt2maxR = 2;


    //Bininng with theta 
    TH2F* h_theta_xB;  // 
    TH2D* h_xB_Q2_quant; // histo with theta bin 
    //TH2D* h_xB_Q2_quant = nullptr;      // new histogram
    //Re binning using y and W for some reason 
    TH2D* h_y_W2_quant; // 

  std::vector<double> xB_edges;       // variable edges
  std::vector<double> Q2_edges;



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
    TH1F *h_vertexX;    //= new TH1F("targetVz", "vertex4target", 100, -40, 40) ;
    TH1F *h_vertexY;    //= new TH1F("targetVz", "vertex4target", 100, -40, 40) ;
    TH1F *h_vertexZ_pi;    //= new TH1F("targetVx", "vertex4target", 100, -40, 40) ;
    TH1F *h_DeltaVz;    //= new TH1F("targetVz", "vertex4target", 100, -40, 40) ;
    TH1F *h_pid;    //= new TH1F("pid", "pid", 100, -250, 250) ;
    TH2F *h_xQ2;    //= new TH2F("xQ2", "xQ2", nubin, xminX, xmaxX, nubin, QminX, QmaxX) ;
    TH2F *h_xQ2pos; //= new TH2F("xQ2pos", "xQ2pos", nubin, xminX, xmaxX, nubin, QminX, QmaxX) ;

    TH2F *h_pt2z; //= new TH2F("pt2z", "pt2z", nubin, pt2minX, pt2maxX, nubin, zminX, zmaxX) ;
    

    //define 2D histogrrams that compare coordinates theta and phi versus electron and hadron momenta (in X axis)
    TH2F *h_thetaP_el;
    TH2F *h_phiP_el;
    TH2F *h_thetaP_had;
    TH2F *h_phiP_had;

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
    TH1F *h_vertexXMC;    //= new TH1F("targetVz", "vertex4target", 100, -40, 40) ;
    TH1F *h_vertexYMC;    //= new TH1F("targetVz", "vertex4target", 100, -40, 40) ;

    TH2F *h_Q2comp;    //= new TH2F("xQ2", "xQ2", nubin, xminX, xmaxX, nubin, QminX, QmaxX) ;

    //momentum histograms
    TH1F *h_px_el;  //= new TH1F("px_ele", "px_ele", 100, 0, 10) ;
    TH1F *h_py_el;  //= new TH1F("py_ele", "py_ele", 100, 0, 10) ;
    TH1F *h_pz_el;  //= new TH1F("pz_ele", "pz_ele", 100, 0, 10) ;
    TH1F *h_ptot_el;    //= new TH1F("ptotal_ele", "ptotal_ele", 100, 0, 10) ;
    TH1F *h_px_pi;  //= new TH1F("px_pro", "px_pro", 100, 0, 10) ;
    TH1F *h_py_pi;  //= new TH1F("py_pro", "py_pro", 100, 0, 10) ;
    TH1F *h_pz_pi;  //= new TH1F("pz_pro", "pz_pro", 100, 0, 10) ;
    TH1F *h_ptot_pi;    //= new TH1F("ptotal_pro", "ptotal_pro", 100, 0, 10) ;
    TH1F *h_theta_el;   //= new TH1F("theta_ele", "theta_ele", 100, 0, 3.5) ;
    TH1F *h_phi_el; //= new TH1F("phi_ele", "phi_ele", 100, 0, 6.5) ;
    TH1F *h_theta_pi;   //= new TH1F("theta_pro", "theta_pro", 100, 0, 3.5) ;
    TH1F *h_phi_pi; //= new TH1F("phi_pro", "phi_pro", 100, 0, 6.5) ;
    TH2F *h_polcoord_el;
    TH2F *h_polcoord_pi;


    //Energy plots 
    TH1F *h_E_el;  //= new TH1F("E_ele", "E_ele", 100, 0, 10) ;
    TH1F *h_E_pi;  //= new TH1F("E_pro", "E_pro", 100, 0, 10) ;
    TH2F *h_E_el_theta; //= new TH2F("E_ele_theta", "E_ele_theta", 100, 0, 50, 100, 0, 10) ;
    TH2F *h_E_pi_theta; //= new TH2F("E_pro_theta", "E_pro_theta", 100, 0, 50, 100, 0, 10) ;
    TH2F *h_E_el_phi;   //= new TH2F("E_ele_phi", "E_ele_phi", 100, 0, 360, 100, 0, 10) ;
    TH2F *h_E_pi_phi;   //= new TH2F("E_pro_phi", "E_pro_phi", 100, 0, 360, 100, 0, 10) ;

    //pre and pos histos 
    TH1F *h_Vz_pre;
    TH1F *h_Vz_post;
    TH1F *h_Q2_pre;
    TH1F *h_Q2_post;
    TH1F *h_y_pre;
    TH1F *h_y_post;
    TH1F *h_nu_pre;
    TH1F *h_nu_post;
    TH1F *h_W2_pre;
    TH1F *h_W2_post;
    TH1F *h_z_pre;
    TH1F *h_z_post;
    TH1F *h_pt2_pre;
    TH1F *h_pt2_post;
    TH2F *h_sampl_el_pre;
    TH2F *h_sampl_el_post;


    //Calo histograms
    TH1F *h_lu; // = new TH1F("lu", "lu", 100, 0, 400) ;    //units cm can be zoomed in to 40 
    TH1F *h_lv; // = new TH1F("lv", "lv", 100, 0, 400) ;    //units cm can be zoomed in to 40
    TH1F *h_lw; // = new TH1F("lw", "lw", 100, 0, 400) ;    //units cm can be zoomed in to 40
    TH1F *h_epcal ; // = new TH1F("epcal", "epcal", 100, 0, 1) ;    
    TH1F *h_eecalin ; // = new TH1F("eecalin", "eecalin", 100, 0, 1) ;
    TH1F *h_epcalout ; // = new TH1F("epcalout", "epcalout", 100, 0, 1) ;
    TH2F *h_calXY; // = new TH2F("calxy", "calxy", 100, -400, 400, 100, -400, 400) ;
    TH2F *h_calEall; // = new TH2F("calEall", "calEall", 100, 0, 1, 100, 0, 1) ;
    TH1F *h_calSector; // = new TH1F("calSector", "calSector", 100, 0, 8) ;

    //Cherenkov histograms
    TH1F *h_Nphe15; // = new TH1F("Nphe15", "Nphe15", 100, 0, 60) ;
    TH1F *h_Nphe16; // = new TH1F("Nphe16", "Nphe16", 100, 0, 60) ;
    TH1F *h_Nphe15pos; // = new TH1F("Nphe15pos", "Nphe15pos", 100, 0, 60) ;
    TH1F *h_Nphe16pos; // = new TH1F("Nphe16pos", "Nphe16pos", 100, 0, 60) ;

    //Helicity histograms
    TH1F *h_helicity; // = new TH1F("helicity", "helicity", 100, -1, 1) ;
    TH1F *h_helicity_raw; // = new TH1F("helicity_pos", "helicity_pos", 100, -1, 1) ;
    //4 Asymmetries 
    TH1F *h_phih_plus; // = new TH1F("phih_plus", "phih_plus", 10, 0, 360) ;
    TH1F *h_phih_minus; // = new TH1F("phih_minus", "phih_minus", 10, 0, 360) ;
    TH1 *h_BSA; // = new TH1F("BSA", "BSA", 10, 0, 360) ; 

    //Chi2 histograms
    TH1F *h_chi2_el; // = new TH1F("chi2", "chi2", 100, 0, 10) ;
    TH1F *h_chi2_pi; // = new TH1F("chi2", "chi2", 100, 0, 10) ;
    TH2F *h_chi2_pid_pi;

    //pion chi2 fixing
    TH2F * h_luthetael; // = new TH2F("luthetapi", "luthetapi", 100, 0, 400, 100, 0, 80) ; 


    //plot for theta binning 
    TH2F *h_xQ2poscuts;


    //sampling fraction 
    TH2F *h_sampl_el; //



            //===Lines for Cuts(TBD)===// 
    //==variable lines==//

    //==calorimeter lines==//

    //==cherenkov lines==//
    //pending TBD



    //add more histograms for other variables here
    

    TH1F *R_nu_el;
    TH1F *R_nu_had;
    TH1F *R_z;
    TH1F *R_pt2;

    // histograms for monitoring resolution on theta and phi coordinates for at leats electron PER SECTOR only 6 sectors (check calo info en fait;) )
    TH1F *h_theta_el_sec1;
    TH1F *h_theta_el_sec2;
    TH1F *h_theta_el_sec3;
    TH1F *h_theta_el_sec4;
    TH1F *h_theta_el_sec5;
    TH1F *h_theta_el_sec6;
    TH1F *h_phi_el_sec1;
    TH1F *h_phi_el_sec2;
    TH1F *h_phi_el_sec3;
    TH1F *h_phi_el_sec4;
    TH1F *h_phi_el_sec5;
    TH1F *h_phi_el_sec6;
std::vector<float> theta_sector[6];
std::vector<float> phi_sector[6];

TH1F* h_phi_res[6];
TH1F* h_theta_res[6];


    //this part is for the unfolding application
    //enum {  NX_DAT = 5, NY_DAT = 3 };
    //enum {  NX_DAT = 7, NY_DAT = 3 };
    enum { NX_BIS = 5, NY_BIS = 3, NX_DAT = 7, NY_DAT = 3 }; // doubleing th nbr of histos...
    // we will apply unfolding on the initial DAT histo. BIS is used as decoy, for comparison with unfolded results (témoin?)

    // IMPORTANT here I believe the DAT histo should be filled with the meas binning (so large, coarse(?))
    //this for the correction to make sense I believe. so Im adding below the comment the proper values of NX_DAT and NY_DAT
    // and xEdgesDAT; thEdgesDAT are not hte same EITHER, need to be changed as well !!!
      //declare bin edges in fct to be called i another fct 
     struct HistoDefsBIS {
      //definition of decoy histogram edges 
      static const double* XEdgesBIS() {
        static const double xEdgesBIS[NX_BIS + 1] = {0.1, 0.11, 0.15, 0.19, 0.29,1.0};
        return xEdgesBIS;
      }
      static const double* ThEdgesBIS() {
        static const double thEdgesBIS[NY_BIS + 1] = {7.6, 8.8, 11.0, 23.0};
        return thEdgesBIS;
      }
    };
    static TH2F* MakeHistoBIS(const std::string& name, const std::string& title) {
      TH2F* h = new TH2F(name.c_str(), title.c_str(), NX_BIS, HistoDefsBIS::XEdgesBIS(), NY_BIS, HistoDefsBIS::ThEdgesBIS());
      h->GetYaxis()->SetRangeUser(0.0, 30.0);
      h->SetDirectory(nullptr); //this should avoid collisions ig MANDATORY ?
      return h;
    }




    struct HistoDefs {
      static const double* XEdgesDAT() {
//        static const double xEdgesREC[NX_REC + 1] = {0.075, 0.105, 0.13, 0.16, 0.20, 0.25, 0.36,1.0};
        static const double xEdgesDAT[NX_DAT + 1] = {0.1, 0.11, 0.13, 0.16, 0.20, 0.25, 0.36,1.0}; // values come from meas histo (REC) in monunfolding class
        //static const double xEdgesDAT[NX_DAT + 1] = {0.075, 0.11, 0.15, 0.19, 0.29, 1.0};
        return xEdgesDAT;
      }
      static const double* ThEdgesDAT() {
        static const double thEdgesDAT[NY_DAT + 1] = {7.6, 8.4, 10.3, 23.0};  //from REC or "meas" values of edges  
        //static const double thEdgesDAT[NY_DAT + 1] = {5.0, 8.8, 11.0, 27.0};
        return thEdgesDAT;
      }
    };
    //stating fct to create histos with previoulsy stablished edges 4 custom binning
    static TH2F* MakeHistoDAT(const std::string& name, const std::string& title) {
        //declaraiton here allows to call for targetname.str() in the constructor in cpp 
        //this method of declaration should be propagated also to the simulated data ig 
      TH2F* h = new TH2F(name.c_str(), title.c_str(),
                         NX_DAT, HistoDefs::XEdgesDAT(),
                         NY_DAT, HistoDefs::ThEdgesDAT());
      h->GetYaxis()->SetRangeUser(0.0, 30.0);
      h->SetDirectory(nullptr); //this should avoid collisions ig MANDATORY ?
      return h;
    }
    TH2F *h_xB_thetaelDAT;
    TH2F *h_xB_thetaelBIS;
    TH1F *h_thetaelDAT_1D;




        //==default histograms (no cuts)==//
    //physics histos
   // TH1F *h_Q2_defaultD = new TH1F("Q2_defaultD", "Q2_defaultD", Constants::bin_default ,Constants::Qmin_default , Constants::Qmax_default);
   // TH1F *h_xb_defaultD= new TH1F("xb_defaultD", "xb_defaultD", Constants::bin_default, Constants::xmin_default, Constants::xmax_default) ;
   // TH1F *h_y_defaultD= new TH1F ("y_defaultD" , "y_defaultD" , Constants::bin_default, Constants::Ymin_default, Constants::Ymax_default) ;
   // TH1F *h_nu_defaultD= new TH1F("nu_defaultD", "nu_defaultD", Constants::bin_default, Constants::numin_default, Constants::numax_default) ;
   // TH1F *h_nu_defaultD_had= new TH1F("nu_defaultDhad", "nu_defaultDhad ", Constants::bin_default, Constants::numin_default, Constants::numax_default) ; //(useful?? TBD)
   // TH1F *h_W2_defaultD= new TH1F("W2_defaultD", "W2_defaultD", Constants::bin_default, Constants::Wmin_default, 30) ;
   // TH1F *h_z_defaultD= new TH1F("z_defaultD", "z_defaultD", Constants::bin_default, Constants::zmin_default, Constants::zmax_default) ;
   // TH1F *h_pt2_defaultD= new TH1F("pt2_defaultD", "pt2_defaultD", Constants::bin_default, Constants::pt2min_default, Constants::pt2max_default) ;
   // TH1F *h_phih_defaultD= new TH1F("phih_defaultD", "phih_defaultD", Constants::phihbin_default, Constants::phihmin_default, Constants::phihmax_default) ;
//
   // //coord histos
   // TH1F *h_phiel_defaultD = new TH1F("phiel_defaultD", "phiel_defaultD", Constants::bin_default, Constants::phielmin_default, Constants::phielmax_default) ;
   // TH1F *h_thetael_defaultD = new TH1F("thetael_defaultD", "thetael_defaultD", Constants::bin_default, Constants::thetaelmin_default, Constants::thetaelmax_default) ;
   // TH1F *h_phipi_defaultD = new TH1F("phipi_defaultD", "phipi_defaultD", Constants::bin_default, Constants::phipimin_default, Constants::phipimax_default) ;
   // TH1F *h_thetapi_defaultD = new TH1F("thetapi_defaultD", "thetapi_defaultD", Constants::bin_default, Constants::thetapimin_default, Constants::thetapimax_default) ;
   // 
//
//
   // //vertex histos
   // TH1F *h_Vx_defaultD= new TH1F("Vx_defaultD", "Vx_defaultD", Constants::bin_default, Constants::Vxmin_default, Constants::Vxmax_default) ;
   // TH1F *h_Vy_defaultD= new TH1F("Vy_defaultD", "Vy_defaultD", Constants::bin_default, Constants::Vymin_default, Constants::Vymax_default) ;
   // TH1F *h_Vz_defaultD= new TH1F("Vz_defaultD", "Vz_defaultD", Constants::bin_default, Constants::Vzmin_default, Constants::Vzmax_default) ;
//
   // //calo histos
   // TH2F *h_calxy_defaultD= new TH2F("calxy_defaultD", "calxy_defaultD", Constants::bin_default, Constants::cal_xmin, Constants::cal_xmax, Constants::bin_default, Constants::cal_ymin, Constants::cal_ymax) ;
    //define a 2D histo with lu (lv or lw) in X axis anf E/p in Y axis (check what E/p is? TBD)


    //commenting histos that dont belong here. DELETE
    //TH1F *h_Q2posel = new TH1F("Q2truePOS", "Q2truePOS", 10, QminX , QmaxX);
    //TH1F *h_nuposel= new TH1F("nuposel", "nuposel", 10,numinX,numaxX) ;
    //TH1F *h_Q2_had = new TH1F("Q2_had", "Q2_had",10, QminX , QmaxX);
    //TH1F *h_nu_had= new TH1F("nu_had", "nu_had", 10, xminX, xmaxX) ;
    //TH1F *h_z_had= new TH1F("z_had", "z_had", 10, zminX, zmaxX) ;
    //TH1F *h_pt2_had= new TH1F("pt2_had", "pt2_had", 10, pt2minX, pt2maxX) ;


    //TH1F *h_xbpre= new TH1F("xbpre", "xbpre", nubin, xminX, xmaxX) ;
    //TH1F *h_ypre= new TH1F ("ypre" , "ypre" , nubin, yminX, ymaxX) ;
    //TH1F *h_nupre= new TH1F("nupre", "nupre", nubin,numinX,numaxX) ;
    //TH1F *h_W2pre= new TH1F("W2pre", "W2pre", nubin, WminX, WmaxX) ;
    //TH1F *h_zpre= new TH1F("zpre", "zpre", nubin, zminX, zmaxX) ;
    //TH1F *h_pt2pre= new TH1F("pt2pre", "pt2pre", nubin, pt2minX, pt2maxX) ;
    //TH1F *h_phihpre= new TH1F("phihpre", "phihpre", nubin, phihminX, phihmaxX) ;




    TFile outputFile;
    //creating a root file 




    //counters for large regions 
    int counterLargeD_regA = 0;
    int counterLargeD_regB = 0;
    int counterLargeD_regC = 0;
    int counterLargeD_regD = 0;
    int counterLargeD_regE = 0;
    int counterLargeD_regF = 0;
    int counterLargeD_regG = 0;
    int counterLargeD_regH = 0;
    int counterLargeD_regI = 0;
    int counterLargeD_regJ = 0;
    int counterLargeD_regK = 0;
    int counterLargeD_regL = 0;
    int counterLargeD_regM = 0;
    int counterLargeD_regN = 0;
    int counterLargeD_regO = 0;
    int counterLargeD_regP = 0;
    int counterLargeD_regQ = 0;
    int counterLargeD_regR = 0;
    int counterLargeD_regS = 0;
    
    int counterLargeA_regA = 0;
    int counterLargeA_regB = 0;
    int counterLargeA_regC = 0;
    int counterLargeA_regD = 0;
    int counterLargeA_regE = 0;
    int counterLargeA_regF = 0;
    int counterLargeA_regG = 0;
    int counterLargeA_regH = 0;
    int counterLargeA_regI = 0;
    int counterLargeA_regJ = 0;
    int counterLargeA_regK = 0;
    int counterLargeA_regL = 0;
    int counterLargeA_regM = 0;
    int counterLargeA_regN = 0;
    int counterLargeA_regO = 0;
    int counterLargeA_regP = 0;
    int counterLargeA_regQ = 0;
    int counterLargeA_regR = 0; 
    int counterLargeA_regS = 0;


    //couters for few region 
    int counterFewD_regA = 0;
    int counterFewD_regB = 0;
    int counterFewD_regC = 0;
    int counterFewD_regD = 0;
    int counterFewD_regE = 0;
    int counterFewD_regF = 0;
    int counterFewD_regG = 0;
    int counterFewD_regH = 0;
    int counterFewD_regI = 0;
    int counterFewD_regJ = 0;
    int counterFewD_regK = 0;
    int counterFewD_regL = 0;  

    int counterFewA_regA = 0;
    int counterFewA_regB = 0;
    int counterFewA_regC = 0;
    int counterFewA_regD = 0;
    int counterFewA_regE = 0;
    int counterFewA_regF = 0;
    int counterFewA_regG = 0;
    int counterFewA_regH = 0;
    int counterFewA_regI = 0;
    int counterFewA_regJ = 0;
    int counterFewA_regK = 0;
    int counterFewA_regL = 0;  


};

#endif 
