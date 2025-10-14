#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLegend.h>
#include <vector>
#include <TCanvas.h>
#include <TPDF.h>
#include "Event.h" 
#include "Monunfold.h"
//#include "RooUnfoldResponse.h"
//#include "RooUnfoldBayes.h"


#include "CutSet.h"
#include "constants.h"

    //int nubin = 100;  

Monunfold::Monunfold(CutSet a, const std::string& targetName)
    : cut1(a), targetName(targetName),

    h_Q2(new TH1F(("U_Q2_" + targetName).c_str(), "UQ2", nubin, QminX, QmaxX)),
    h_xb(new TH1F(("U_xb_" + targetName).c_str(), "Uxb", nubin, xminX, xmaxX)),
    h_y(new TH1F(("U_y_" + targetName).c_str(), "Uy", nubin, yminX, ymaxX)),
    h_nu(new TH1F(("U_nu_" + targetName).c_str(), "Unu", nubin, numinX, numaxX)),
    h_W2(new TH1F(("U_W2_" + targetName).c_str(), "UW2", nubin, WminX, WmaxX)),
    h_z(new TH1F(("U_z_" + targetName).c_str(), "Uz", nubin, zminX, zmaxX)),
    h_pt2(new TH1F(("U_pt2_" + targetName).c_str(), "Upt2", nubin, pt2minX, pt2maxX)),
    h_phih(new TH1F(("U_phih_" + targetName).c_str(), "Uphih", nubin, phihminX, phihmaxX)),
    h_vertexZ(new TH1F(("U_targetVz_" + targetName).c_str(), "Uvertex4target", 100, -10, 0)),
    h_vertexY(new TH1F(("U_targetVy_" + targetName).c_str(), "Uvertex4target", 100, -10, 10)),
    h_vertexX(new TH1F(("U_targetVx_" + targetName).c_str(), "Uvertex4target", 100, -10, 10)),
    
    h_Q2MC(new TH1F(("U_Q2_MC" + targetName).c_str(), "Q2MC", nubin, QminX, QmaxX)),
    h_xbMC(new TH1F(("U_xb_MC" + targetName).c_str(), "xbMC", nubin, xminX, xmaxX)),
    h_yMC(new TH1F(("U_y_MC" + targetName).c_str(), "yMC", nubin, yminX, ymaxX)),
    h_nuMC(new TH1F(("U_nu_MC" + targetName).c_str(), "nuMC", nubin, numinX, numaxX)),
    h_W2MC(new TH1F(("U_W2_MC" + targetName).c_str(), "W2MC", nubin, WminX, 20)),
    h_zMC(new TH1F(("U_z_MC" + targetName).c_str(), "zMC", nubin, zminX, zmaxX)),
    h_pt2MC(new TH1F(("U_pt2_MC" + targetName).c_str(), "pt2MC", nubin, pt2minX, pt2maxX)),
    h_phihMC(new TH1F(("U_phih_MC" + targetName).c_str(), "phihMC", nubin, phihminX, phihmaxX)),
    h_vertexZMC(new TH1F(("U_targetVz_MC" + targetName).c_str(), "vertex4targetMC", 100, -20, 10)),
    h_vertexYMC(new TH1F(("U_targetVy_MC" + targetName).c_str(), "vertex4targetMC", 100, -10, 10)),
    h_vertexXMC(new TH1F(("U_targetVx_MC" + targetName).c_str(), "vertex4targetMC", 100, -10, 10)),
    h_pt2zMC(new TH2F(("U_pt2z_MC" + targetName).c_str(), "pt2zMC", nubin, pt2minX, pt2maxX, nubin, zminX, zmaxX)),

    h_xBthREC(new TH2F(("U_xBthREC_" + targetName).c_str(), "xBthREC", nubin, 0,30, nubin, xminX, xmaxX)),
    h_xBthMC(new TH2F(("U_xBthMC_" + targetName).c_str(), "xBthMC", nubin, 0,30, nubin, xminX, xmaxX)),
    h_xQREC(new TH2F(("U_xQ2REC_" + targetName).c_str(), "xQ2REC", nubin, xminX, xmaxX, nubin, QminX, QmaxX)),
    h_xQMC(new TH2F(("U_xQ2MC_" + targetName).c_str(), "xQ2MC", nubin, xminX, xmaxX, nubin, QminX, QmaxX)),

    h_Q2comp(new TH2F(("U_Q2comp_" + targetName).c_str(), "Q2comp", nubin, QminX, QmaxX, nubin, QminX, QmaxX)),
    h_xbcomp(new TH2F(("U_xbcomp_" + targetName).c_str(), "xbcomp", nubin, xminX, xmaxX, nubin, xminX, xmaxX)),
    h_ycomp(new TH2F(("U_ycomp_" + targetName).c_str(), "ycomp", nubin, yminX, ymaxX, nubin, yminX, ymaxX)),
    h_nucomp(new TH2F(("U_nucomp_" + targetName).c_str(), "nucomp", nubin, numinX, numaxX, nubin, numinX, numaxX)),
    h_W2comp(new TH2F(("U_W2comp_" + targetName).c_str(), "W2comp", nubin, WminX, 20, nubin, WminX, 20)),
    h_zcomp(new TH2F(("U_zcomp_" + targetName).c_str(), "zcomp", nubin, zminX, zmaxX, nubin, zminX, zmaxX)),
    h_pt2comp(new TH2F(("U_pt2comp_" + targetName).c_str(), "pt2comp", nubin, pt2minX, pt2maxX, nubin, pt2minX, pt2maxX)),
    h_phihcomp(new TH2F(("U_phihcomp_" + targetName).c_str(), "phihcomp", nubin, phihminX, phihmaxX, nubin, phihminX, phihmaxX)),
    h_vertexZcomp(new TH2F(("U_targetVzcomp_" + targetName).c_str(), "vertex4targetcomp", 100, -10, 10.0, 100, -10, 10.0)),

    h_px_el(new TH1F(("px_ele_" + targetName).c_str(), "px_ele", nubin, 0, 3)),
    h_py_el(new TH1F(("py_ele_" + targetName).c_str(), "py_ele", nubin, 0, 3)),
    h_pz_el(new TH1F(("pz_ele_" + targetName).c_str(), "pz_ele", nubin, 0, 10)),
    h_ptot_el(new TH1F(("ptot_ele_" + targetName).c_str(), "ptot_ele", nubin, 0, 10)),
    h_Delta_px_el(new TH1F(("px_ele_" + targetName).c_str(), "px_ele", nubin, -1, 1)),
    h_Delta_py_el(new TH1F(("py_ele_" + targetName).c_str(), "py_ele", nubin, -1, 1)),
    h_Delta_pz_el(new TH1F(("pz_ele_" + targetName).c_str(), "pz_ele", nubin, -1, 1)),
    h_Delta_ptot_el(new TH1F(("ptot_ele_" + targetName).c_str(), "ptot_ele", nubin, -1, 1)),

    h_px_pi(new TH1F(("px_had_" + targetName).c_str(), "px_had", nubin, 0, 3)),
    h_py_pi(new TH1F(("py_had_" + targetName).c_str(), "py_had", nubin, 0, 3)),
    h_pz_pi(new TH1F(("pz_had_" + targetName).c_str(), "pz_had", nubin, 0, 10)),
    h_ptot_pi(new TH1F(("ptot_had_" + targetName).c_str(), "ptot_had", nubin, 0, 10)),
    h_theta_el(new TH1F(("theta_el" + targetName).c_str(), "theta", nubin, 0, 30)),
    h_phi_el(new TH1F(("phi_el" + targetName).c_str(), "phi", nubin, 0, 360)),
    h_theta_pi(new TH1F(("theta_pi" + targetName).c_str(), "theta", nubin, 0, 150)),
    h_Delta_theta_el(new TH1F(("theta_el" + targetName).c_str(), "theta", nubin, -2,2)),
    h_Delta_phi_el(new TH1F(("phi_el" + targetName).c_str(), "phi", nubin, -5,5)),

    h_phi_pi(new TH1F(("phi_pi" + targetName).c_str(), "phi", nubin, 0, 360)),
    h_px_elMC(new TH1F(("px_ele_MC" + targetName).c_str(), "px_eleMC", nubin, 0, 2)),
    h_py_elMC(new TH1F(("py_ele_MC" + targetName).c_str(), "py_eleMC", nubin, 0, 2)),
    h_pz_elMC(new TH1F(("pz_ele_MC" + targetName).c_str(), "pz_eleMC", nubin, 0, 10)),
    h_ptot_elMC(new TH1F(("ptot_ele_MC" + targetName).c_str(), "ptot_eleMC", nubin, 0, 10)),
    h_px_piMC(new TH1F(("px_had_MC" + targetName).c_str(), "px_hadMC", nubin, 0, 2)),
    h_py_piMC(new TH1F(("py_had_MC" + targetName).c_str(), "py_hadMC", nubin, 0, 2)),
    h_pz_piMC(new TH1F(("pz_had_MC" + targetName).c_str(), "pz_hadMC", nubin, 0, 10)),
    h_ptot_piMC(new TH1F(("ptot_had_MC" + targetName).c_str(), "ptot_hadMC", nubin, 0, 10)),
    h_theta_elMC(new TH1F(("theta_elMC" + targetName).c_str(), "thetaMC", nubin, 0, 30)),
    h_phi_elMC(new TH1F(("phi_elMC" + targetName).c_str(), "phiMC", nubin, 0, 360)),
    h_theta_piMC(new TH1F(("theta_piMC" + targetName).c_str(), "thetaMC", nubin, 0, 150)),
    h_phi_piMC(new TH1F(("phi_piMC" + targetName).c_str(), "phiMC", nubin, 0, 360)),
    h_E_el(new TH1F(("E_el" + targetName).c_str(), "E", nubin, 0, 10)),
    h_E_pi(new TH1F(("E_pi" + targetName).c_str(), "E", nubin, 0, 10)),
    h_E_elMC(new TH1F(("E_elMC" + targetName).c_str(), "EelMC", nubin, 0, 10)),
    h_E_piMC(new TH1F(("E_piMC" + targetName).c_str(), "EpiMC", nubin, 0, 10)),
    h_deltaphipi(new TH1F(("U_deltaphipi_" + targetName).c_str(), "deltaphipi", 100, -5, 5)),
    h_deltathetapi(new TH1F(("U_deltathetapi_" + targetName).c_str(), "deltathetapi", 100, -5, 5)),

    h_px_elcomp(new TH2F(("U_px_elcomp_" + targetName).c_str(), "px_elcomp", nubin, 0, 2, nubin, 0, 2)),
    h_py_elcomp(new TH2F(("U_py_elcomp_" + targetName).c_str(), "py_elcomp", nubin, 0, 2, nubin, 0, 2)),
    h_pz_elcomp(new TH2F(("U_pz_elcomp_" + targetName).c_str(), "pz_elcomp", nubin, 0, 10, nubin, 0, 10)),
    h_ptot_elcomp(new TH2F(("U_ptot_elcomp_" + targetName).c_str(), "ptot_elcomp", nubin, 0, 10, nubin, 0, 10)),
    h_px_picomp(new TH2F(("U_px_picomp_" + targetName).c_str(), "px_picomp", nubin, 0, 2, nubin, 0, 2)),
    h_py_picomp(new TH2F(("U_py_picomp_" + targetName).c_str(), "py_picomp", nubin, 0, 2, nubin, 0, 2)),
    h_pz_picomp(new TH2F(("U_pz_picomp_" + targetName).c_str(), "pz_picomp", nubin, 0, 10, nubin, 0, 10)),
    h_ptot_picomp(new TH2F(("U_ptot_picomp_" + targetName).c_str(), "ptot_picomp", nubin, 0, 10, nubin, 0, 10)),
    h_theta_elcomp(new TH2F(("U_theta_elcomp_" + targetName).c_str(), "theta_elcomp", nubin, 0, 30, nubin, 0, 30)),
    h_phi_elcomp(new TH2F(("U_phi_elcomp_" + targetName).c_str(), "phi_elcomp", nubin, 0, 360, nubin, 0, 360)),
    h_theta_picomp(new TH2F(("U_theta_picomp_" + targetName).c_str(), "theta_picomp", nubin, 0, 150, nubin, 0, 150)),
    h_phi_picomp(new TH2F(("U_phi_picomp_" + targetName).c_str(), "phi_picomp", nubin, 0, 360, nubin, 0, 360)),
    h_E_elcomp(new TH2F(("U_E_elcomp_" + targetName).c_str(), "E_elcomp", nubin, 0, 10, nubin, 0, 10)),
    h_E_picomp(new TH2F(("U_E_picomp_" + targetName).c_str(), "E_picomp", nubin, 0, 10, nubin, 0, 10)),
    h_vxcomp(new TH2F(("U_vxcomp_" + targetName).c_str(), "vxcomp", nubin, -10, 10, nubin, -10, 10)),
    h_vycomp(new TH2F(("U_vycomp_" + targetName).c_str(), "vycomp", nubin, -10, 10, nubin, -10, 10)), 
    h_chi2_el(new TH1F(("U_chi2_el_" + targetName).c_str(), "chi2_el", nubin, 0, 10)),
    h_chi2_pi(new TH1F(("U_chi2_pi_" + targetName).c_str(), "chi2_pi", nubin, 0, 10)),
    h_px_piDelta(new TH1F (("U_Delta_px_pi_" + targetName).c_str(), "Delta_px_pi", nubin, -1, 1)),
    h_py_piDelta(new TH1F (("U_Delta_py_pi_" + targetName).c_str(), "Delta_py_pi", nubin, -1, 1)),
    h_pz_piDelta(new TH1F (("U_Delta_pz_pi_" + targetName).c_str(), "Delta_pz_pi", nubin, -1, 1)),
    h_ptot_piDelta(new TH1F (("U_Delta_ptot_pi_" + targetName).c_str(), "Delta_ptot_pi", nubin, -1, 1)),
    h_theta_piDelta(new TH1F (("U_Delta_theta_pi_" + targetName).c_str(), "Delta_theta_pi", nubin, -10, 10)),
    h_phi_piDelta(new TH1F (("U_Delta_phi_pi_" + targetName).c_str(), "Delta_phi_pi", nubin, -100, 100)),
    h_E_piDelta(new TH1F (("U_Delta_E_pi_" + targetName).c_str(), "Delta_E_pi", nubin, -1, 1)),

    //h_xB_thetaelMC(Monunfold::HistoMC()),  // already initialized
    //h_xB_thetaelREC(Monunfold::HistoREC()),  // already initialized
    h_xB_thetaelMC(Monunfold::MakeHistoMC( Form("h_xB_thetaelMC_%s", targetName.c_str()), Form("MC (%s): x_{B} vs #theta_{e};x_{B};#theta_{e} [deg]", targetName.c_str()))),
    h_xB_thetaelREC(Monunfold::MakeHistoSIM( Form("h_xB_thetaelREC_%s", targetName.c_str()), Form("REC (%s): x_{B} vs #theta_{e};x_{B};#theta_{e} [deg]", targetName.c_str()))),

    tEv_(new TTree(("tEv_" + targetName).c_str(), "Event Tree")),

    //h_xB_thetaelMC(new TH2F(("h_xB_thetaelMC_"+ targetName).c_str(),"MC: x_{B} vs #theta_{e};x_{B};#theta_{e} [deg]", nxMC, xEdgesMC.data(), nyMC, thEdgesMC.data())),
    //h_xB_thetaelREC(new TH2F(("h_xB_thetaelREC_"+ targetName).c_str(),"REC: x_{B} vs #theta_{e};x_{B};#theta_{e} [deg]", nxREC, xEdgesREC.data(), nyREC, thEdgesREC.data())),
    //h_xB_thetaelREC(new TH2F(("h_xB_thetaelREC_"+ targetName).c_str(),"REC: x_{B} vs #theta_{e};x_{B};#theta_{e} [deg]", nxREC, xEdgesREC.data(), 4, 0 , 30)),
    h_thetaelMC_1D(new TH1F(("U_thetaelMC_1D_" + targetName).c_str(), "thetaelMC_1D", 100, 0, 30.0)),
    h_thetaelREC_1D(new TH1F(("U_thetaelREC_1D_" + targetName).c_str(), "thetaelREC_1D", 100, 0, 30.0)),
    h_thetaelcalcMC_1D(new TH1F(("U_thetaelcalcMC_1D_" + targetName).c_str(), "thetaelcalcMC_1D", 100, 0, 30.0)),
    h_thetaelcalcREC_1D(new TH1F(("U_thetaelcalcREC_1D_" + targetName).c_str(), "thetaelcalcREC_1D", 100, 0, 30.0)),
    h_thxbuniform(new TH2F(("U_thxbuniform_" + targetName).c_str(), "thxbuniform", 100, 0,1, 100, 0, 30.0)),


    h_evtnbrdiff(new TH1F(("U_evtnbrdiff_" + targetName).c_str(), "evtnbrdiff", 100, -10, 10.0)),




      counterel_R(0) {
        tEv_->Branch("evnbr", &evnum_, "evtnbr/I"); 
        tEv_->Branch("xb_MC", &Br_xbMC, "xb_MC/D");
        tEv_->Branch("th_el_MC", &Br_thMC, "th_el_MC/D");
        tEv_->Branch("xb_REC", &Br_xbREC, "xb_REC/D");
        tEv_->Branch("th_el_REC", &Br_thREC, "th_el_REC/D");
        tEv_->Branch("Q2_MC", &Br_Q2MC, "Q2_MC/D");
        tEv_->Branch("Q2_REC", &Br_Q2REC, "Q2_REC/D");
        tEv_->Branch("y_MC", &Br_yMC, "y_MC/D");
        tEv_->Branch("y_REC", &Br_yREC, "y_REC/D");
        tEv_->Branch("W2_MC", &Br_W2MC, "W2_MC/D");
        tEv_->Branch("W2_REC", &Br_W2REC, "W2_REC/D");
        tEv_->Branch("nu_MC", &Br_nuMC, "nu_MC/D");
        tEv_->Branch("nu_REC", &Br_nuREC, "nu_REC/D");


    // Add more histograms as needed

}


//void Monunfold::SetResponseIndex(const Event& event){
//    int indexsetter = 0;
//    if (cut1.PassCutsDetectors(event)) {
//        if (cut1.PassCutsElectrons(event)==true) {
//            int Q2evt = event.GetQ2();
//            int xbevt = event.Getxb();
//            
//            for (const Particle& hadron : event.GetHadrons()) {
//                if (cut1.PassCutsHadrons(hadron)) {
//                    if (hadron.GetPID() == Constants::PION_PLUS_PID) {
//                    int zevt = hadron.Getz();
//                    int pt2evt = hadron.Getpt2();
//                    int phihevt = hadron.Getphih();
//                    //if (check values of Q2, xb,phih, pt2, z in that order to check their window. depending on the window you will assign an index ) {
//                    }
//                }
//            }
//        }
//    }
//}
//


void Monunfold::CheckLargeBins(const Event& event){
    int targetType = event.GetTargetType();
    if ( cut1.PassCutsElectrons(event)==true && cut1.PassCutsDetectors(event)==true ) { 
        if (event.Getxb()>Constants::LargexblowAB && event.Getxb()<Constants::LargexbhighAB){   //for regions A and B 
            if (event.GetQ2()>Constants::LARGEQlowA && event.GetQ2()<Constants::LARGEQhighA) {      //REGION A
                counterLargeD_regA ++;
                //std::cout << "bump REGION A" << std::endl;
            }
            if (event.GetQ2()>Constants::LARGEQlowB && event.GetQ2()<Constants::LARGEQhighB) {      //REGION B
                counterLargeD_regB ++;
            }
        }
        else if (event.Getxb()>Constants::LARGEQhighB && event.Getxb()<Constants::LargexbhighCDE ){ //for regions C D and E
            if (event.GetQ2()>Constants::LARGEQlowC && event.GetQ2()<Constants::LARGEQhighC) {      //REGION C
                counterLargeD_regC ++;
            }
            if (event.GetQ2()>Constants::LARGEQlowD && event.GetQ2()<Constants::LARGEQhighD) {      //REGION D
                counterLargeD_regD ++;
            }
            if (event.GetQ2()>Constants::LARGEQlowE && event.GetQ2()<Constants::LARGEQhighE) {      //REGION E
                counterLargeD_regE ++;
            }
        }    
        else if (event.Getxb()>Constants::LargexbhighCDE && event.Getxb()<Constants::LargexbhighFGH){    //for region F G and H
            if (event.GetQ2()>Constants::LARGEQlowF && event.GetQ2()<Constants::LARGEQhighF) {      //REGION F
                counterLargeD_regF ++;
            }
            if (event.GetQ2()>Constants::LARGEQlowG && event.GetQ2()<Constants::LARGEQhighG) {      //REGION G
                counterLargeD_regG ++;
            }
            if (event.GetQ2()>Constants::LARGEQlowH && event.GetQ2()<Constants::LARGEQhighH) {      //REGION H
                counterLargeD_regH ++;
            }
        } 
        else if (event.Getxb()>Constants::LargexbhighFGH && event.Getxb()<Constants::LargexbhighIJK){    // for regions I J and K 
            if (event.GetQ2()>Constants::LARGEQlowI && event.GetQ2()<Constants::LARGEQhighI) {      //REGION I
                counterLargeD_regI ++;
            }
            if (event.GetQ2()>Constants::LARGEQlowJ && event.GetQ2()<Constants::LARGEQhighJ) {      //REGION J
                counterLargeD_regJ ++;
            }
            if (event.GetQ2()>Constants::LARGEQlowK && event.GetQ2()<Constants::LARGEQhighK) {      //REGION K
                counterLargeD_regK ++;
            }
        }   
        else if (event.Getxb()>Constants::LargexbhighIJK && event.Getxb()<Constants::LargexbhighLMN){    //for regions L M and N
            if (event.GetQ2()>Constants::LARGEQlowL && event.GetQ2()<Constants::LARGEQhighL) {      //REGION L
                counterLargeD_regL ++;
            }
            if (event.GetQ2()>Constants::LARGEQlowM && event.GetQ2()<Constants::LARGEQhighM) {      //REGION M
                counterLargeD_regM ++;
            }
            if (event.GetQ2()>Constants::LARGEQlowN && event.GetQ2()<Constants::LARGEQhighN) {      //REGION N
                counterLargeD_regN ++;
            }
        } 
        else if (event.Getxb()>Constants::LargexbhighLMN && event.Getxb()<Constants::LargexbhighOPQ){    //for regions O P and Q
            if (event.GetQ2()>Constants::LARGEQlowO && event.GetQ2()<Constants::LARGEQhighO) {      //REGION O
                counterLargeD_regO ++;
            }
            if (event.GetQ2()>Constants::LARGEQlowP && event.GetQ2()<Constants::LARGEQhighP) {      //REGION P
                counterLargeD_regP ++;
            }
            if (event.GetQ2()>Constants::LARGEQlowQ && event.GetQ2()<Constants::LARGEQhighQ) {      //REGION Q
                counterLargeD_regQ ++;
            }
        } 
        else if (event.Getxb()>Constants::LargexbhighOPQ && event.Getxb()<Constants::LargexbhighRS){    //for region R S
            if (event.GetQ2()>Constants::LARGEQlowR && event.GetQ2()<Constants::LARGEQhighR) {      //REGION R
                counterLargeD_regR ++;
            }
            if (event.GetQ2()>Constants::LARGEQlowS && event.GetQ2()<Constants::LARGEQhighS) {      //REGION S
                counterLargeD_regS ++;
            }
        }
    }     
}

void Monunfold::PrintRegionCounters( ){
    std::cout << "RegA nb = " << counterLargeD_regA << std::endl;
    std::cout << "RegB nb = " << counterLargeD_regB << std::endl;
    std::cout << "RegC nb = " << counterLargeD_regC << std::endl;
    std::cout << "RegD nb = " << counterLargeD_regD << std::endl;
    std::cout << "RegE nb = " << counterLargeD_regE << std::endl;
    std::cout << "RegF nb = " << counterLargeD_regF << std::endl;
    std::cout << "RegG nb = " << counterLargeD_regG << std::endl;
    std::cout << "RegH nb = " << counterLargeD_regH << std::endl;
    std::cout << "RegI nb = " << counterLargeD_regI << std::endl;
    std::cout << "RegJ nb = " << counterLargeD_regJ << std::endl;
    std::cout << "RegK nb = " << counterLargeD_regK << std::endl;
    std::cout << "RegL nb = " << counterLargeD_regL << std::endl;
    std::cout << "RegM nb = " << counterLargeD_regM << std::endl;
    std::cout << "RegN nb = " << counterLargeD_regN << std::endl;
    std::cout << "RegO nb = " << counterLargeD_regO << std::endl;
    std::cout << "RegP nb = " << counterLargeD_regP << std::endl;
    std::cout << "RegQ nb = " << counterLargeD_regQ << std::endl;
    std::cout << "RegR nb = " << counterLargeD_regR << std::endl;
    std::cout << "RegS nb = " << counterLargeD_regS << std::endl;
}


int Monunfold::SetResponseIndex(const Event& event) {
    int indexsetter = 0;

    if (!cut1.PassCutsDetectors(event)) return 0;
    if (!cut1.PassCutsElectrons(event)) return 0;
    //how to handle if return is 0. Do not add int matrx

    double Q2evt = event.GetQ2();
    double xbevt = event.Getxb();

    int regionID = -1;
    if (Q2evt > 1 && Q2evt <= 2.56 && xbevt > 0.06 && xbevt <= 0.2) regionID = 0; //regA  //regionID acts like a counter useful when building the idx using regs 
    else if (Q2evt > 2.56 && Q2evt <= 5 && xbevt > 0.2 && xbevt <= 0.5) regionID = 1; //regB
    else if (Q2evt > 1 && Q2evt <= 1.22 && xbevt > 0.06 && xbevt <= 0.2) regionID = 2; //regC
    else if (Q2evt > 1.22 && Q2evt <= 3.21 && xbevt > 0.2 && xbevt <= 0.5) regionID = 3; //regD
    else return 0;

    for (const Particle& hadron : event.GetHadrons()) {
        if (!cut1.PassCutsHadrons(hadron)) continue;
        if (hadron.GetPID() != Constants::PION_PLUS_PID) continue;

        double zevt = hadron.Getz();
        double pt2evt = hadron.Getpt2();
        double phihevt = hadron.Getphih();

        int zbin = int((zevt - 0.3) / 0.1);
        if (zbin < 0 || zbin >= 4) continue;

        int pt2bin = int(pt2evt / 0.3);
        if (pt2bin < 0 || pt2bin >= 4) continue;

        int phibin = int(phihevt / 90.0);
        if (phibin < 0 || phibin >= 4) continue;

        indexsetter = zbin + 4 * pt2bin + 16 * phibin + 64 * regionID + 1;
        break;  // We found the index for the first valid pion, return it
    }

    return indexsetter;
}


// ===================== Unfolding (fine A..S) =====================

int Monunfold::U_RegionIndexFine_(const Event& e) const {
  const double xb = e.Getxb();
  const double Q2 = e.GetQ2();

  // A,B (xB in AB)
  if (xb > Constants::LargexblowAB && xb < Constants::LargexbhighAB) {
    if (Q2 > Constants::LARGEQlowA && Q2 < Constants::LARGEQhighA) return 0; // A
    if (Q2 > Constants::LARGEQlowB && Q2 < Constants::LARGEQhighB) return 1; // B
    return -1;
  }
  // C,D,E
  if (xb > Constants::LargexbhighAB && xb < Constants::LargexbhighCDE) {
    if (Q2 > Constants::LARGEQlowC && Q2 < Constants::LARGEQhighC) return 2;
    if (Q2 > Constants::LARGEQlowD && Q2 < Constants::LARGEQhighD) return 3;
    if (Q2 > Constants::LARGEQlowE && Q2 < Constants::LARGEQhighE) return 4;
    return -1;
  }
  // F,G,H
  if (xb > Constants::LargexbhighCDE && xb < Constants::LargexbhighFGH) {
    if (Q2 > Constants::LARGEQlowF && Q2 < Constants::LARGEQhighF) return 5;
    if (Q2 > Constants::LARGEQlowG && Q2 < Constants::LARGEQhighG) return 6;
    if (Q2 > Constants::LARGEQlowH && Q2 < Constants::LARGEQhighH) return 7;
    return -1;
  }
  // I,J,K
  if (xb > Constants::LargexbhighFGH && xb < Constants::LargexbhighIJK) {
    if (Q2 > Constants::LARGEQlowI && Q2 < Constants::LARGEQhighI) return 8;
    if (Q2 > Constants::LARGEQlowJ && Q2 < Constants::LARGEQhighJ) return 9;
    if (Q2 > Constants::LARGEQlowK && Q2 < Constants::LARGEQhighK) return 10;
    return -1;
  }
  // L,M,N
  if (xb > Constants::LargexbhighIJK && xb < Constants::LargexbhighLMN) {
    if (Q2 > Constants::LARGEQlowL && Q2 < Constants::LARGEQhighL) return 11;
    if (Q2 > Constants::LARGEQlowM && Q2 < Constants::LARGEQhighM) return 12;
    if (Q2 > Constants::LARGEQlowN && Q2 < Constants::LARGEQhighN) return 13;
    return -1;
  }
  // O,P,Q
  if (xb > Constants::LargexbhighLMN && xb < Constants::LargexbhighOPQ) {
    if (Q2 > Constants::LARGEQlowO && Q2 < Constants::LARGEQhighO) return 14;
    if (Q2 > Constants::LARGEQlowP && Q2 < Constants::LARGEQhighP) return 15;
    if (Q2 > Constants::LARGEQlowQ && Q2 < Constants::LARGEQhighQ) return 16;
    return -1;
  }
  // R,S
  if (xb > Constants::LargexbhighOPQ && xb < Constants::LargexbhighRS) {
    if (Q2 > Constants::LARGEQlowR && Q2 < Constants::LARGEQhighR) return 17;
    if (Q2 > Constants::LARGEQlowS && Q2 < Constants::LARGEQhighS) return 18;
    return -1;
  }
  return -1;
}

void Monunfold::U_InitUnfold(const std::string& tag, const CutSet& dataCuts) {
  U_tag_ = tag;
  U_dataCuts_ = dataCuts;
  const int nFine = 19; // A..S

  U_h_true_.reset(new TH1D(("U_h_true_"+tag).c_str(), ("MC truth ["+tag+"]").c_str(), nFine, 0.5, nFine+0.5));
  U_h_meas_.reset(new TH1D(("U_h_meas_"+tag).c_str(), ("REC meas ["+tag+"]").c_str(),  nFine, 0.5, nFine+0.5));
  U_h_data_.reset(new TH1D(("U_h_data_"+tag).c_str(), ("RGD data ["+tag+"]").c_str(),  nFine, 0.5, nFine+0.5));

  U_h_true_->SetDirectory(nullptr);
  U_h_meas_->SetDirectory(nullptr);
  U_h_data_->SetDirectory(nullptr);

  static const char* labs = "ABCDEFGHIJKLMNOPQRS";
  for (int i=1;i<=nFine;i++) {
    U_h_true_->GetXaxis()->SetBinLabel(i, std::string(1,labs[i-1]).c_str());
    U_h_meas_->GetXaxis()->SetBinLabel(i, std::string(1,labs[i-1]).c_str());
    U_h_data_->GetXaxis()->SetBinLabel(i, std::string(1,labs[i-1]).c_str());
  }

  // Response built on the same binnings uncomment
  //U_response_.reset(new RooUnfoldResponse(U_h_meas_.get(), U_h_true_.get()));
}

void Monunfold::U_FillSimPair(const Event& rec, const Event& mc, double w) {
  // Apply your detector-level analysis cuts to REC (as you do for RGD)
  if (!(U_dataCuts_.PassCutsElectrons(rec) && U_dataCuts_.PassCutsDetectors(rec))) return;

  const int k_meas = U_RegionIndexFine_(rec);
  const int k_true = U_RegionIndexFine_(mc);
  if (k_true>=0) U_h_true_->Fill(k_true+1, w);
  if (k_meas>=0 && k_true>=0) {
    //U_response_->Fill(k_meas+1, k_true+1, w);
    U_h_meas_->Fill(k_meas+1, w);
  } else if (k_true>=0 && k_meas<0) {
    //U_response_->Miss(k_true+1, w);
  } else if (k_true<0 && k_meas>=0) {
    //U_response_->Fake(k_meas+1, w);
    U_h_meas_->Fill(k_meas+1, w);
  }
}

void Monunfold::U_FillSimTruthOnly(const Event& mc, double w) {
  const int k_true = U_RegionIndexFine_(mc);
  if (k_true>=0) {
    U_h_true_->Fill(k_true+1, w);
    //U_response_->Miss(k_true+1, w);
  }
}

void Monunfold::U_FillSimRecoOnly(const Event& rec, double w) {
  if (!(U_dataCuts_.PassCutsElectrons(rec) && U_dataCuts_.PassCutsDetectors(rec))) return;
  const int k_meas = U_RegionIndexFine_(rec);
  if (k_meas>=0) {
    U_h_meas_->Fill(k_meas+1, w);
    //U_response_->Fake(k_meas+1, w);
  }
}

void Monunfold::U_FillData(const Event& data, double w) {
  if (!(U_dataCuts_.PassCutsElectrons(data) && U_dataCuts_.PassCutsDetectors(data))) return;
  const int k_meas = U_RegionIndexFine_(data);
  if (k_meas>=0) U_h_data_->Fill(k_meas+1, w);
}

//TH1D* Monunfold::U_UnfoldBayes(int nIter) {
//  //RooUnfoldBayes unfold(U_response_.get(), U_h_data_.get(), nIter);
//  unfold.SetVerbose(0);
//  TH1D* h_unf = (TH1D*) unfold.Hreco(RooUnfold::kCovariance);
//  h_unf->SetName(("U_h_unfold_"+U_tag_).c_str());
//  h_unf->SetTitle(("Unfolded ["+U_tag_+"]").c_str());
//  U_cov_ = std::make_unique<TMatrixD>(unfold.ErecoV());
//  return h_unf;
//}
//
//TH1D* Monunfold::U_Refold(const TH1D* htruth_like) const {
//  TH1* href = U_response_->ApplyToTruth(htruth_like, ("U_hrefold_"+U_tag_).c_str());
//  href->SetTitle(("Refolded ["+U_tag_+"]").c_str());
//  return (TH1D*) href;
//}

void Monunfold::U_SaveAll(const std::string& fname, TH1D* h_unfold) {
  TFile f(fname.c_str(), "RECREATE");
  if (U_h_true_)  U_h_true_->Write();
  if (U_h_meas_)  U_h_meas_->Write();
  if (U_h_data_)  U_h_data_->Write();
  //if (U_response_) U_response_->Write(("U_response_"+U_tag_).c_str());
  //if (h_unfold) h_unfold->Write();
  //if (U_cov_)   U_cov_->Write(("U_cov_unfold_"+U_tag_).c_str());
  f.Close();
}


void Monunfold::FillOnlyVz(const Event& event) {
    //Temporary function to fill only vertexZ histograms
    //this used to check on vz for both (only true) True Carbon Targets
    //no need to cjeck MC banks on true data. 
    if (cut1.PassCutVzselection(event) ==true) {
        h_vertexZ->Fill(event.GetVz());
        //h_vertexZMC->Fill(event.GetVzMC());
        //h_vertexZcomp->Fill(event.GetVz(), event.GetVzMC());
    }
    //h_vertexZMC->Fill(event.GetVzMC());
    //h_vertexZcomp->Fill(event.GetVz(), event.GetVzMC());
    //then create a function to draw. Delete this function once used
}

void Monunfold::DrawOnlyVz(Monunfold& munfoldC2 , const std::string filename) {
    //Temporary function to draw only vertexZ histograms
    //this used to check on vz for both (only true) True Carbon Targets
    //no need to cjeck MC banks on true data. 
    TCanvas c1("c1", "c1", 800, 600);
    c1.Divide(1, 2);
    c1.cd(1);
    h_vertexZ->Draw();
    munfoldC2.h_vertexZ->SetLineColor(kBlack);
    munfoldC2.h_vertexZ->Draw("same");
    c1.cd(2);
    //h_vertexZMC->SetLineColor(kRed);
    h_vertexZMC->Draw("same");
    //h_vertexZcomp->Draw("colz");
    c1.Print((filename + ".pdf").c_str());
    //then create a function to draw. Delete this function once used
}

void Monunfold::FillHistogramsNoCuts(const Event& event) {
    // Fill histograms with no cuts
    h_Q2->Fill(event.GetQ2());
    h_xb->Fill(event.Getxb());
    h_y->Fill(event.Gety());
    h_nu->Fill(event.Getnu());
    h_W2->Fill(event.GetW2());
    h_vertexZ->Fill(event.GetVz());


    //h_xQ2->Fill(event.Getxb(), event.GetQ2());
    //h_calXY->Fill(event.GetCalX(), event.GetCalY());
    //h_lu->Fill(event.Getlu());
    //h_lv->Fill(event.Getlv());
    //h_lw->Fill(event.Getlw());
    //h_epcal->Fill(event.GetEpcal());
    //h_thetaelectron->Fill(event.GetThetaElectron());
    //h_rapport->Fill(event.GetAcosyada());
    ////h_eecalin->Fill(event.GetEcalin());
    ////h_epcalout->Fill(event.GetEcalout());
    ////h_calEall->Fill(event.GetEpcal(), event.GetEcalin()+event.GetEcalout());
            //h_eecalin->Fill(event.GetEcalin());
    //if (event.GetEcalin()>0.01){h_eecalin->Fill(event.GetEcalin());}
    ////h_epcalout->Fill(event.GetEcalout());
    //if (event.GetEcalout()>0.01){h_epcalout->Fill(event.GetEcalout());}
    //if (event.GetEcalout()>0.01 & event.GetEcalout()>0.01){h_calEall->Fill(event.GetEpcal(), event.GetEcalin()+event.GetEcalout());}
    ////h_calEall->Fill(event.GetEpcal(), event.GetEcalin()+event.GetEcalout());
    //h_Nphe15->Fill(event.Getnphe15());
    //h_Nphe16->Fill(event.Getnphe16());
    //h_calSector->Fill(event.GetCalSector());
    //h_helicity->Fill(event.GetHel());
    //h_helicity_raw->Fill(event.GetHelRaw());
    Particle electron = event.GetElectron();
    //h_px_el->Fill(electron.GetMomentum().X());
    //h_py_el->Fill(electron.GetMomentum().Y());
    //h_pz_el->Fill(electron.GetMomentum().Z());
    //h_ptot_el->Fill(electron.GetMomentum().P());
    //h_theta_el->Fill(electron.GetMomentum().Theta()*180/Constants::PI);
    //h_phi_el->Fill(electron.GetMomentum().Phi()*180/Constants::PI +180);
    //h_polcoord_el->Fill(electron.GetMomentum().Theta()*180/Constants::PI, electron.GetMomentum().Phi()*180/Constants::PI +180);
    //h_E_el->Fill(electron.GetMomentum().E());
    //h_E_el_theta->Fill(electron.GetMomentum().Theta()*180/Constants::PI, electron.GetMomentum().E());
    //h_E_el_phi->Fill(electron.GetMomentum().Phi()*180/Constants::PI +180, electron.GetMomentum().E());

    
    
    for (const Particle& hadron : event.GetHadrons()) {
        if (hadron.GetPID() == Constants::PION_PLUS_PID  ){   //adding this condition for pion+ and erasing the condit at evtprocessr
        h_z->Fill(hadron.Getz());
        h_pt2->Fill(hadron.Getpt2());
        h_phih->Fill(hadron.Getphih());
        //h_px_pi->Fill(hadron.GetMomentum().X());
        //h_py_pi->Fill(hadron.GetMomentum().Y());
        //h_pz_pi->Fill(hadron.GetMomentum().Z());
        //h_ptot_pi->Fill(hadron.GetMomentum().P());
        //h_theta_pi->Fill(hadron.GetMomentum().Theta()*180/Constants::PI);
        //h_phi_pi->Fill(hadron.GetMomentum().Phi()*180/Constants::PI +180);
        //h_polcoord_pi->Fill(hadron.GetMomentum().Theta()*180/Constants::PI, hadron.GetMomentum().Phi()*180/Constants::PI +180);
        //h_E_pi->Fill(hadron.GetMomentum().E());
        //h_E_pi_theta->Fill(hadron.GetMomentum().Theta()*180/Constants::PI, hadron.GetMomentum().E());
        //h_E_pi_phi->Fill(hadron.GetMomentum().Phi()*180/Constants::PI +180, hadron.GetMomentum().E());
        //h_vertexZ_pi->Fill(hadron.GetParticleVertexZ());
        ////std::cout<<"hadron vertex z = "<<hadron.GetParticleVertexZ()<<std::endl;
        ////std::cout<<"electron vertex z = "<<event.GetVz()<<std::endl;
//
        //h_DeltaVz->Fill(event.GetVz()-hadron.GetParticleVertexZ());
        //if (event.GetHel() > 0){
        //    h_phih_plus->Fill(hadron.Getphih());
        //}
        //else if (event.GetHel() < 0){
        //    h_phih_minus->Fill(hadron.Getphih());
        //}
        //
        //
        ////
        ////std::cout << " kinhad " << hadron.Getz()<<" , " << hadron.Getpt2()<<" , " << hadron.Getphih()<<std::endl;  
        ////
        }
    }
}

void Monunfold::FillDISforUnfoldMC(const Event& event) {
    //std::cout<<"entering MC function"<<std::endl;
    //if (cut1.PassCutsDetectors(event)) {
        //std::cout<<"event passed det cuts"<<std::endl;
        if (cut1.PassCutsElectronsMC(event)==true) {
            // dont forget to convert from rad to deg
            //h_xB_thetaelMC->Fill(event.GetxbMC(), event.GetThetaElectronMC()* 180.0 / Constants::PI);
            h_vertexZMC->Fill(event.GetVzMC());
            h_xB_thetaelMC->Fill(event.GetxbMC(), event.GetThetaElectronMC()* 180.0 / Constants::PI);
            h_thetaelMC_1D->Fill(event.MCelectron.GetMomentum().Theta()*180/Constants::PI);
            h_thetaelcalcMC_1D->Fill(event.GetThetaElectronMC()* 180.0 / Constants::PI);
                      //std::cout<<"xb mc = "<<event.GetxbMC()<<" theta el mc = "<<event.MCelectron.GetPID()<<std::endl;
        }
    //}
}
void Monunfold::SetrangesMC(){
    h_xB_thetaelMC->GetXaxis()->SetRangeUser(0.0, 1.0);
    h_xB_thetaelMC->GetYaxis()->SetRangeUser(0.0, 30.0);
}
void Monunfold::SetrangesREC(){
    h_xB_thetaelREC->GetXaxis()->SetRangeUser(0.0, 1.0);
    h_xB_thetaelREC->GetYaxis()->SetRangeUser(0.0, 30.0);
}

void Monunfold::debughisto(){
    //h_xB_thetaelREC->
    std::cout << "nbins x = " << h_xB_thetaelREC->GetNbinsX() << std::endl;
    std::cout << "nbins y = " << h_xB_thetaelREC->GetNbinsY() << std::endl;
    std::cout << "suuposed nbr of bins = " << NX_REC*NY_REC << std::endl;
    std::cout << "supp bins X = " << NX_REC <<   std::endl;
    std::cout << "supp bins Y = " << NY_REC <<   std::endl;
}
void Monunfold::FillDISforUnfoldREC(const Event& event) {
    if (cut1.PassCutsDetectors(event)) {
        if (cut1.PassCutsElectrons(event)==true) {
            h_vertexZ->Fill(event.GetVz());
            //h_xB_thetaelREC->Fill(event.Getxb(), event.GetThetaElectron()* 180.0 / Constants::PI);
            h_xB_thetaelREC->Fill(event.Getxb(), event.electron.GetMomentum().Theta()* 180.0 / Constants::PI);
            //std::cout<<"xb rec = "<<event.Getxb()<<" theta el rec = "<<event.GetThetaElectron()<<std::endl;
            h_thetaelREC_1D->Fill(event.electron.GetMomentum().Theta()*180/Constants::PI);
            h_thetaelcalcREC_1D->Fill(event.GetThetaElectron()* 180.0 / Constants::PI);
            h_thxbuniform->Fill(event.Getxb(), event.electron.GetMomentum().Theta()*180/Constants::PI);

        }
    }
}
void Monunfold::saveDISforUnfoldRoot(const std::string& filenameDIS) {
    TFile* rootFile = new TFile((filenameDIS + ".root").c_str(), "RECREATE");
    if (h_xB_thetaelMC)  h_xB_thetaelMC->Write();
    if (h_xB_thetaelREC)  h_xB_thetaelREC->Write();
    if (h_thetaelMC_1D)  h_thetaelMC_1D->Write();
    if (h_thetaelREC_1D)  h_thetaelREC_1D->Write();
    if (h_thetaelcalcMC_1D)  h_thetaelcalcMC_1D->Write();
    if (h_thetaelcalcREC_1D)  h_thetaelcalcREC_1D->Write();
    if (h_thxbuniform)  h_thxbuniform->Write();
    if (h_vertexZMC)  h_vertexZMC->Write();
    if (h_vertexZ)  h_vertexZ->Write();
    rootFile->Close();

    delete rootFile; 
}

void Monunfold::PrintFAKE(){
    std::cout<<"Number of 'FAKE' events = "<<counterFAKEMC<<std::endl;
    std::cout<<"Number of 'MISS' events = "<<counterMISSREC<<std::endl;
    std::cout<<"Number of 'RESPONSE' events = "<<counterMATCHREC<<std::endl;
    std::cout<<"Number of 'TOTAL' events = "<<counterTOTAL<<std::endl;
}



void Monunfold::ProperFillRECMC(const Event& event_MC, const Event& event_REC, int option ){
    //ig folwoing the option system. should consider only evts that pass REC and MC cuts I suppose... 
    //is this useful ? 
    double currentxbMC = -.1;
    double currentQ2MC = -.1;
    double currentthMC = -.1;
    double currentyMC = -.1;
    double currentW2MC = -.1;
    double currentnuMC = -.1;
    double currentxbREC = -.1;
    double currentQ2REC = -.1;
    double currentthREC = -.1;
    double currentyREC = -.1;
    double currentW2REC = -.1;
    double currentnuREC = -.1;
    //defining them here helps reset the values to 0 in each call of the function (evt) 
    if (cut1.PassCutsElectronsMC(event_MC)==true) {
            currentxbMC = event_MC.GetxbMC();
            currentQ2MC = event_MC.GetQ2MC();
            currentyMC = event_MC.GetyMC();
            currentW2MC = event_MC.GetW2MC();
            currentnuMC = event_MC.GetnuMC();
            currentthMC = event_MC.MCelectron.GetMomentum().Theta()*  180.0 / Constants::PI;
            counterPassMCproper++;
        //std::cout<<"[bump] MC"<<std::endl;
    }
    if (option == true ){   
        if ( cut1.PassCutsElectrons(event_REC)==true){
            currentQ2REC = event_REC.GetQ2();
            currentxbREC = event_REC.Getxb();
            currentyREC = event_REC.Gety();
            currentW2REC = event_REC.GetW2();
            currentnuREC = event_REC.Getnu();
            currentthREC = event_REC.electron.GetMomentum().Theta()* 180.0 / Constants::PI;

        }
                
    }

    double sumMC = currentQ2MC + currentxbMC + currentyMC + currentW2MC + currentnuMC;
    double sumREC = currentQ2REC + currentxbREC + currentyREC + currentW2REC + currentnuREC;  

    if (sumREC > 0.0 ){
        h_Q2->Fill(currentQ2REC);
        h_xb->Fill(currentxbREC);
        h_y->Fill(currentyREC);
        h_nu->Fill(currentnuREC);
        h_W2->Fill(currentnuREC);
        h_xQREC->Fill(currentxbREC, currentQ2REC);
        h_xBthREC->Fill(currentxbREC, currentthREC);
        Br_xbREC = currentxbREC;
        Br_Q2REC = currentQ2REC;
        Br_yREC = currentyREC;
        Br_nuREC = currentnuREC;
        Br_W2REC = currentW2REC;
        
    }
    if (sumMC > 0.0 ){
        h_Q2MC->Fill(currentxbMC);
        h_xbMC->Fill(currentQ2MC);
        h_yMC->Fill(currentyMC);
        h_nuMC->Fill(currentnuMC);
        h_W2MC->Fill(currentW2MC);
        h_xQMC->Fill(currentxbMC, currentQ2MC);
        h_xBthMC->Fill(currentxbMC, currentthMC);
        Br_xbMC = currentxbMC;
        Br_Q2MC = currentQ2MC;
        Br_yMC = currentyMC;
        Br_nuMC = currentnuMC;
        Br_W2MC = currentW2MC;
    }
    tEv_->Fill();
    
}

void Monunfold::ProperSaveRECMC(const std::string& filenameRECMC) {
    std::cout<<"[COUNT] passingxb MC = "<<counterPassMCproper<<std::endl;
    
    TFile* rootFile = new TFile((filenameRECMC + ".root").c_str(), "RECREATE");
    if (h_Q2)  h_Q2->Write();
    if (h_xb)  h_xb->Write();
    if (h_y)  h_y->Write();
    if (h_nu)  h_nu->Write();
    if (h_W2)  h_W2->Write();

    if (h_Q2MC)  h_Q2MC->Write();
    if (h_xbMC)  h_xbMC->Write();
    if (h_yMC)  h_yMC->Write();
    if (h_nuMC)  h_nuMC->Write();
    if (h_W2MC)  h_W2MC->Write();

    if (h_xQREC)  h_xQREC->Write();
    if (h_xQMC)  h_xQMC->Write();
    if (h_xBthREC)  h_xBthREC->Write();
    if (h_xBthMC)  h_xBthMC->Write();

    if (!tEv_) return;
    rootFile->cd();
    tEv_->Write();  // writes "tEv_<targetName>"

    rootFile->Close();

    delete rootFile; 
}


void Monunfold::FillTreeEvt(const Event& event_MC  , const Event& event_REC , int option  ){
    //this version is WORKING
    //add cuts, you can definnitelt handle different cuts for bot hevts, needs to recheck cut handling in MC
    
    //if option is true, register value, if option is false, then set REC to 0 
    //( will later be treated as miiss hen checking the TTRee) 
 
    //reset values of the variables before branching, sometimes we need REC valeues at zer to rgister MISS
    double currentxbMC = 0.0;
    double currentthMC = 0.0;
    
    if (cut1.PassCutsElectronsMC(event_MC)==true){
        //std::cout<<"[BUMP] event passed mc cuts"<<std::endl;
        counterTOTAL++; //not a real total, just total of considered MC evts 

        currentthMC = event_MC.MCelectron.GetMomentum().Theta()*180/Constants::PI;
        currentxbMC = event_MC.GetxbMC();
        double currentxbREC = 0.0;
        double currentthREC = 0.0;
        if (option == true ){   
            if (cut1.PassCutsDetectors(event_REC)==true && cut1.PassCutsElectrons(event_REC)==true){
//                std::cout<<"[BUMP] event passed rec cuts"<<std::endl;
                currentthREC = event_REC.electron.GetMomentum().Theta()* 180.0 / Constants::PI;
                currentxbREC = event_REC.Getxb();
                counterMATCHREC++;

            }
        }
        if (currentthREC == 0 && currentxbREC == 0){
            counterMISSREC++;
        }
        Br_xbREC = currentxbREC;
        Br_thREC = currentthREC;
        Br_xbMC = currentxbMC;
        Br_thMC = currentthMC;


        tEv_->Fill();

    }
    

}


void Monunfold::WriteTTree(const std::string& filenameTREE) {
    TFile* rootFile = new TFile((filenameTREE + ".root").c_str(), "RECREATE");
    if (!tEv_) return;
    rootFile->cd();
    tEv_->Write();  // writes "tEv_<targetName>"

    rootFile->Close();
    delete rootFile; 


}

void Monunfold::FillHistogramswCuts(const Event& event) {
    if (cut1.PassCutsDetectors(event)) {

        if (cut1.PassCutsElectrons(event)==true) {
        // Fill histograms after passing detector cuts
        h_vertexZ->Fill(event.GetVz()); // Fill vertex Z histogram

        // Fill other histograms related to electron variables
        h_Q2->Fill(event.GetQ2());
        h_xb->Fill(event.Getxb());
        h_y->Fill(event.Gety());
        h_nu->Fill(event.Getnu());
        h_W2->Fill(event.GetW2());
        
        // Fill electron momentum histograms
        Particle electron = event.GetElectron();
        h_px_el->Fill(electron.GetMomentum().X());
        h_py_el->Fill(electron.GetMomentum().Y());
        h_pz_el->Fill(electron.GetMomentum().Z());
        h_ptot_el->Fill(sqrt(electron.GetMomentum().Mag2()));
        h_theta_el->Fill(electron.GetMomentum().Theta() * 180.0 / Constants::PI);
        h_phi_el->Fill(electron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
        h_chi2_el->Fill(electron.Getchi2());
        //h_polcoord_el->Fill(electron.GetMomentum().Theta() * 180.0 / Constants::PI, electron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
        //h_E_el->Fill(electron.GetMomentum().E());
        //h_E_el_theta->Fill(electron.GetMomentum().Theta() * 180.0 / Constants::PI, electron.GetMomentum().E());
        //h_E_el_phi->Fill(electron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0, electron.GetMomentum().E());

        // Loop over hadrons and fill histograms for those passing hadron cuts
        for (const Particle& hadron : event.GetHadrons()) {
            if (cut1.PassCutsHadrons(hadron)) {
                if (hadron.GetPID() == Constants::PION_PLUS_PID) {
                    h_z->Fill(hadron.Getz());
                    h_pt2->Fill(hadron.Getpt2());
                    h_phih->Fill(hadron.Getphih());
                    h_px_pi->Fill(hadron.GetMomentum().X());
                    h_py_pi->Fill(hadron.GetMomentum().Y());
                    h_pz_pi->Fill(hadron.GetMomentum().Z());
                    h_ptot_pi->Fill(hadron.GetMomentum().P());
                    h_theta_pi->Fill(hadron.GetMomentum().Theta() * 180.0 / Constants::PI);
                    h_phi_pi->Fill(hadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
                    h_chi2_pi->Fill(hadron.Getchi2());
                    //h_polcoord_pi->Fill(hadron.GetMomentum().Theta() * 180.0 / Constants::PI, hadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
                    //h_E_pi->Fill(hadron.GetMomentum().E());
                    //h_E_pi_theta->Fill(hadron.GetMomentum().Theta() * 180.0 / Constants::PI, hadron.GetMomentum().E());
                    //h_E_pi_phi->Fill(hadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0, hadron.GetMomentum().E());
                    //h_vertexZ_pi->Fill(hadron.GetParticleVertexZ());
                    //h_DeltaVz->Fill(event.GetVz() - hadron.GetParticleVertexZ());
                }
            }
        }
        }
    }
}


void Monunfold::FillHistogramswCutsMC(const Event& event) {

    //std::cout<< "bump MC starts"<<std::endl;
    if (cut1.PassCutsElectronsMC(event)) {

        // Fill histograms after passing detector cuts for MC data
        h_vertexZMC->Fill(event.GetVzMC()); // Fill MC vertex Z histogram

        // Fill other histograms related to MC electron variables
        h_Q2MC->Fill(event.GetQ2MC());
        h_xbMC->Fill(event.GetxbMC());
        h_yMC->Fill(event.GetyMC());
        h_nuMC->Fill(event.GetnuMC());
        h_W2MC->Fill(event.GetW2MC());
        //std::cout<< "bump MC el cuts"<<std::endl;
        // Fill MC electron momentum histograms
        Particle electron = event.GetElectron();
        //h_px_el->Fill(electron.GetMomentum().X()); // Assuming the same histogram as in real data
        //h_py_el->Fill(electron.GetMomentum().Y());
        //h_pz_el->Fill(electron.GetMomentum().Z());
        //h_ptot_el->Fill(sqrt(electron.GetMomentum().Mag2()));
        //h_theta_el->Fill(electron.GetMomentum().Theta() * 180.0 / Constants::PI);
        //h_phi_el->Fill(electron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
        //h_polcoord_el->Fill(electron.GetMomentum().Theta() * 180.0 / Constants::PI, electron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
        //h_E_el->Fill(electron.GetMomentum().E());
        //h_E_el_theta->Fill(electron.GetMomentum().Theta() * 180.0 / Constants::PI, electron.GetMomentum().E());
        //h_E_el_phi->Fill(electron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0, electron.GetMomentum().E());

        // Loop over MC hadrons and fill histograms for those passing hadron cuts
        for (const Particle& hadron : event.GetMCHadrons()) {
            if (cut1.PassCutsHadrons(hadron)) {
                if (hadron.GetPID() == Constants::PION_PLUS_PID) {
                    std::cout<< "bump MC hadron cuts"<<std::endl;
                    h_zMC->Fill(hadron.GetzMC());
                    h_pt2MC->Fill(hadron.Getpt2MC());
                    h_phihMC->Fill(hadron.GetphihMC());
                    h_pt2zMC->Fill(hadron.Getpt2MC(), hadron.GetzMC());
                    //h_px_pi->Fill(hadron.GetMomentum().X()); // Assuming the same histogram as in real data
                    //h_py_pi->Fill(hadron.GetMomentum().Y());
                    //h_pz_pi->Fill(hadron.GetMomentum().Z());
                    //h_ptot_pi->Fill(hadron.GetMomentum().P());
                    //h_theta_pi->Fill(hadron.GetMomentum().Theta() * 180.0 / Constants::PI);
                    //h_phi_pi->Fill(hadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
                    //h_polcoord_pi->Fill(hadron.GetMomentum().Theta() * 180.0 / Constants::PI, hadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
                    //h_E_pi->Fill(hadron.GetMomentum().E());
                    //h_E_pi_theta->Fill(hadron.GetMomentum().Theta() * 180.0 / Constants::PI, hadron.GetMomentum().E());
                    //h_E_pi_phi->Fill(hadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0, hadron.GetMomentum().E());
                    //h_vertexZMC->Fill(hadron.GetParticleVertexZ()); // Fill MC hadron vertex Z histogram
                    //h_DeltaVz->Fill(event.GetVz() - hadron.GetParticleVertexZ());
                }
            }
        }
    }
}


void Monunfold::FillHistogramsNoCutsMC(const Event& event) {            //EVERYTHING WORKS DONT DELETE
    // Fill histograms with no cuts
    h_Q2MC->Fill(event.GetQ2MC());
    h_xbMC->Fill(event.GetxbMC());
    h_yMC->Fill(event.GetyMC());
    h_nuMC->Fill(event.GetnuMC());
    h_W2MC->Fill(event.GetW2MC());
    h_vertexZMC->Fill(event.GetVzMC());
    Particle electron = event.GetElectron();
    for (const Particle& MChadron : event.GetMCHadrons()) {
        if (MChadron.GetPID() == Constants::PION_PLUS_PID  ){   
        h_zMC->Fill(MChadron.GetzMC());
        h_pt2MC->Fill(MChadron.Getpt2MC());
        h_phihMC->Fill(MChadron.GetphihMC());
        h_pt2zMC->Fill(MChadron.Getpt2MC(), MChadron.GetzMC());
        }
    }
}


void Monunfold::FillHistComp(const Event& eventsim, const Event& eventmc){
    //h_vxcomp->Fill(eventsim.GetVx(), eventmc.GetVxMC());
    //h_vycomp->Fill(eventsim.GetVy(), eventmc.GetVyMC());

    if (eventsim.GetElectron().GetMomentum().P()> 2 ){

        h_px_el->Fill(eventsim.GetElectron().GetMomentum().X());
        h_py_el->Fill(eventsim.GetElectron().GetMomentum().Y());
        h_pz_el->Fill(eventsim.GetElectron().GetMomentum().Z());
        h_ptot_el->Fill(eventsim.GetElectron().GetMomentum().P());
        h_theta_el->Fill(eventsim.GetElectron().GetMomentum().Theta() * 180.0 / Constants::PI);
        h_phi_el->Fill(eventsim.GetElectron().GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
        h_E_el->Fill(eventsim.GetElectron().GetMomentum().E());
    }
    if (eventmc.GetElectronMC().GetMomentum().P()> 2 ){

        h_px_elMC->Fill(eventmc.GetElectronMC().GetMomentum().X());
        h_py_elMC->Fill(eventmc.GetElectronMC().GetMomentum().Y());
        h_pz_elMC->Fill(eventmc.GetElectronMC().GetMomentum().Z());
        h_ptot_elMC->Fill(eventmc.GetElectronMC().GetMomentum().P());
        h_theta_elMC->Fill(eventmc.GetElectronMC().GetMomentum().Theta() * 180.0 / Constants::PI);
        h_phi_elMC->Fill(eventmc.GetElectronMC().GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
        h_E_elMC->Fill(eventmc.GetElectronMC().GetMomentum().E());
    }
    if (eventsim.GetElectron().GetMomentum().P()> 2 && eventmc.GetElectronMC().GetMomentum().P()> 2 ){
        double deltaphi_el = abs(eventsim.GetElectron().GetMomentum().Phi()  - eventmc.GetElectronMC().GetMomentum().Phi() )*180;
        double deltatheta_el = abs(eventsim.GetElectron().GetMomentum().Theta()  - eventmc.GetElectronMC().GetMomentum().Theta() )*180; 
        if ( deltaphi_el+ deltatheta_el < 10){
            h_Q2comp->Fill(eventsim.GetQ2(), eventmc.GetQ2MC());
            h_xbcomp->Fill(eventsim.Getxb(), eventmc.GetxbMC());
            h_ycomp->Fill(eventsim.Gety(), eventmc.GetyMC());
            h_nucomp->Fill(eventsim.Getnu(), eventmc.GetnuMC());
            h_W2comp->Fill(eventsim.GetW2(), eventmc.GetW2MC());
            h_vertexZcomp->Fill(eventsim.GetVz(), eventmc.GetVzMC());

            h_px_elcomp->Fill(eventsim.GetElectron().GetMomentum().X(), eventmc.GetElectronMC().GetMomentum().X());
            h_py_elcomp->Fill(eventsim.GetElectron().GetMomentum().Y(), eventmc.GetElectronMC().GetMomentum().Y());
            h_pz_elcomp->Fill(eventsim.GetElectron().GetMomentum().Z(), eventmc.GetElectronMC().GetMomentum().Z());
            h_ptot_elcomp->Fill(eventsim.GetElectron().GetMomentum().P(), eventmc.GetElectronMC().GetMomentum().P());
            h_Delta_px_el->Fill(eventmc.GetElectronMC().GetMomentum().X()- eventsim.GetElectron().GetMomentum().X());
            h_Delta_py_el->Fill(eventmc.GetElectronMC().GetMomentum().Y()- eventsim.GetElectron().GetMomentum().Y());
            h_Delta_pz_el->Fill(eventmc.GetElectronMC().GetMomentum().Z()- eventsim.GetElectron().GetMomentum().Z());
            h_Delta_ptot_el->Fill(eventmc.GetElectronMC().GetMomentum().P()- eventsim.GetElectron().GetMomentum().P()); 
            h_theta_elcomp->Fill(eventsim.GetElectron().GetMomentum().Theta() * 180.0 / Constants::PI, eventmc.GetElectronMC().GetMomentum().Theta() * 180.0 / Constants::PI);
            h_phi_elcomp->Fill(eventsim.GetElectron().GetMomentum().Phi() * 180.0 / Constants::PI + 180.0, eventmc.GetElectronMC().GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
            h_Delta_phi_el->Fill( ( eventmc.GetElectronMC().GetMomentum().Phi() - eventsim.GetElectron().GetMomentum().Phi()   )*180 );
            //h_Delta_theta_el->Fill(eventmc.GetElectronMC().GetMomentum().Theta() * 180.0 / Constants::PI - eventsim.GetElectron().GetMomentum().Theta() * 180.0 / Constants::PI  );
            h_Delta_theta_el->Fill(( eventmc.GetElectronMC().GetMomentum().Theta() - eventsim.GetElectron().GetMomentum().Theta()  )*180  );
            h_E_elcomp->Fill(eventsim.GetElectron().GetMomentum().E(), eventmc.GetElectronMC().GetMomentum().E());
        }
    }


    for (const Particle& MChadron: eventmc.GetMCHadrons()){
        if (MChadron.GetPID() == Constants::PION_PLUS_PID  ){   
            h_px_piMC->Fill(MChadron.GetMomentum().X());
            h_py_piMC->Fill(MChadron.GetMomentum().Y());
            h_pz_piMC->Fill(MChadron.GetMomentum().Z());
            h_ptot_piMC->Fill(MChadron.GetMomentum().P());
            h_theta_piMC->Fill(MChadron.GetMomentum().Theta() * 180.0 / Constants::PI);
            h_phi_piMC->Fill(MChadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
            h_E_piMC->Fill(MChadron.GetMomentum().E());
            h_phihMC->Fill(MChadron.GetphihMC());
        }
    }
    for (const Particle& hadron : eventsim.GetHadrons()) {
        if (hadron.GetPID() == Constants::PION_PLUS_PID  ){   
        

        h_px_pi->Fill(hadron.GetMomentum().X());
        h_py_pi->Fill(hadron.GetMomentum().Y());
        h_pz_pi->Fill(hadron.GetMomentum().Z());
        h_ptot_pi->Fill(hadron.GetMomentum().P());
        h_theta_pi->Fill(hadron.GetMomentum().Theta() * 180.0 / Constants::PI);
        h_E_pi->Fill(hadron.GetMomentum().E());
        h_phi_pi->Fill(hadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
        h_phih->Fill(hadron.Getphih());

        }
    }
    for (const Particle& hadron : eventsim.GetHadrons()) {
        for (const Particle& MChadron: eventmc.GetMCHadrons()){
            if (MChadron.GetPID() == Constants::PION_PLUS_PID && hadron.GetPID() == Constants::PION_PLUS_PID ){   
                // Cut on hadrons can be added
                //cut on Delta phi and Delta thetacan be added
                //double deltaphi_pi = 
                //double deltatheta_pi
                h_phihcomp->Fill(hadron.Getphih(), MChadron.GetphihMC());
                h_deltaphipi->Fill(abs(hadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0) - abs (MChadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0));
                h_deltathetapi->Fill(hadron.GetMomentum().Theta() * 180.0 / Constants::PI - MChadron.GetMomentum().Theta() * 180.0 / Constants::PI);
                h_px_picomp->Fill(hadron.GetMomentum().X(), MChadron.GetMomentum().X());
                h_py_picomp->Fill(hadron.GetMomentum().Y(), MChadron.GetMomentum().Y());
                h_pz_picomp->Fill(hadron.GetMomentum().Z(), MChadron.GetMomentum().Z());
                h_ptot_picomp->Fill(hadron.GetMomentum().P(), MChadron.GetMomentum().P());
                h_theta_picomp->Fill(hadron.GetMomentum().Theta() * 180.0 / Constants::PI, MChadron.GetMomentum().Theta() * 180.0 / Constants::PI);
                h_phi_picomp->Fill(hadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0, MChadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0);
                h_E_picomp->Fill(hadron.GetMomentum().E(), MChadron.GetMomentum().E());
                h_zcomp->Fill(hadron.Getz(), MChadron.GetzMC());
                h_pt2comp->Fill(hadron.Getpt2(), MChadron.Getpt2MC());
                h_px_piDelta->Fill(hadron.GetMomentum().X()-MChadron.GetMomentum().X());
                h_py_piDelta->Fill(hadron.GetMomentum().Y()-MChadron.GetMomentum().Y());
                h_pz_piDelta->Fill(hadron.GetMomentum().Z()-MChadron.GetMomentum().Z());
                h_ptot_piDelta->Fill(hadron.GetMomentum().P()-MChadron.GetMomentum().P());
                h_theta_piDelta->Fill(hadron.GetMomentum().Theta() * 180.0 / Constants::PI - MChadron.GetMomentum().Theta() * 180.0 / Constants::PI);
                h_phi_piDelta->Fill(hadron.GetMomentum().Phi() * 180.0 / Constants::PI  - MChadron.GetMomentum().Phi() * 180.0 / Constants::PI );
                //std::cout << " " << std::endl;
                //std::cout << "phiREC = " << hadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0  << std::endl;
                //std::cout << "phiMC = " << MChadron.GetMomentum().Phi() * 180.0 / Constants::PI + 180.0 << std::endl;
                h_E_piDelta->Fill(hadron.GetMomentum().E()-MChadron.GetMomentum().E());



            }


        }
    }


    //h_zcomp->Fill(eventsim.Getz(), eventmc.GetzMC());
    //h_pt2comp->Fill(eventsim.Getpt2(), eventmc.Getpt2MC());
    //h_phihcomp->Fill(eventsim.Getphih(), eventmc.GetphihMC());
    //h_vertexZcomp->Fill(eventsim.GetVz(), eventmc.GetVzMC());
}


void Monunfold::FillHistCompwCuts(const Event& eventsim, const Event& eventmc) {
    double tempQ2sim = 0;
    double tempQ2mc = 0;
    //if (cut1.PassCutsDetectors(eventsim)){
    //if (cut1.PassCutsElectrons(eventsim)==true) {
    if (eventsim.GetQ2()>1.5 ){
        if (eventsim.Gety() > 0.25 && eventsim.Gety() < 0.85 ){
            if (eventsim.GetW2()>2.5 ){
                h_Q2comp->Fill(eventsim.GetQ2(), eventmc.GetQ2MC());
                h_xbcomp->Fill(eventsim.Getxb(), eventmc.GetxbMC());
                h_ycomp->Fill(eventsim.Gety(), eventmc.GetyMC());
                h_nucomp->Fill(eventsim.Getnu(), eventmc.GetnuMC());
                h_W2comp->Fill(eventsim.GetW2(), eventmc.GetW2MC());
                h_vertexZcomp->Fill(eventsim.GetVz(), eventmc.GetVzMC());
                h_evtnbrdiff->Fill(eventsim.GetEvtnbr()-eventmc.GetEvtnbr());
                for (const Particle& MChadron: eventmc.GetMCHadrons()){
                    for (const Particle& hadron : eventsim.GetHadrons()) {
                        if (MChadron.GetPID() == Constants::PION_PLUS_PID  ){   
                            if (cut1.PassCutsHadrons(hadron) == true){ 
                                h_zcomp->Fill(hadron.Getz(), MChadron.GetzMC());
                                h_pt2comp->Fill(hadron.Getpt2(), MChadron.Getpt2MC());
                                h_phihcomp->Fill(hadron.Getphih(), MChadron.GetphihMC());
                            }
                        }
                    }
                }
            }
        }
    }
//    }
}


//void Monunfold::createResponseMatrix(){
//    RooUnfoldResponse responseof2D(h_xB_thetaelREC,h_xB_thetaelMC);
//}



void Monunfold::Fill(const Event& event) {
    // Fill histograms with event data
    h_Q2->Fill(event.GetQ2());
    //std::cout << "Q2: " << event.GetQ2() << std::endl;
    h_xb->Fill(event.Getxb());
    h_y->Fill(event.Gety());
    h_nu->Fill(event.Getnu());
    h_W2->Fill(event.GetW2());
    //h_z->Fill(event.Getz());
    //h_pt2->Fill(event.Getpt2());
    //h_phih->Fill(event.Getphih());
    h_vertexZ->Fill(event.GetVz());
    //h_pid->Fill(event.Getpid());
    //
    //// fill MC histograms with event data
    h_Q2MC->Fill(event.GetQ2MC());
    h_xbMC->Fill(event.GetxbMC());
    h_yMC->Fill(event.GetyMC());
    h_nuMC->Fill(event.GetnuMC());
    h_W2MC->Fill(event.GetW2MC());
    h_vertexZMC->Fill(event.GetVzMC());
    //
    //// fill comparison histograms with event data
    h_Q2comp->Fill(event.GetQ2(), event.GetQ2MC());
    h_xbcomp->Fill(event.Getxb(), event.GetxbMC());
    h_ycomp->Fill(event.Gety(), event.GetyMC());
    h_nucomp->Fill(event.Getnu(), event.GetnuMC());
    h_W2comp->Fill(event.GetW2(), event.GetW2MC());
    h_vertexZcomp->Fill(event.GetVz(), event.GetVzMC());
    
    // issue: MChadron and hadron are two diff loops. Solution is to treat them separately then merge
    //collect z values from event>hadron and event>MChadron (they come from same event)
    std::vector<double> z_values;
    std::vector<double> zMC_values;
    std::vector<double> pt2_values;
    std::vector<double> pt2MC_values;
    std::vector<double> phih_values;
    std::vector<double> phihMC_values;

    //1st loop in hadrons 
    for (const Particle& hadron : event.GetHadrons()) {
        if (hadron.GetPID() == Constants::PION_PLUS_PID  ){
            z_values.push_back(hadron.Getz());
            pt2_values.push_back(hadron.Getpt2());
            phih_values.push_back(hadron.Getphih());
        }
    }
    
    //2nd loop over MChadrons

    for (const Particle& MChadron : event.GetMCHadrons()) {
        if (MChadron.GetPID() == Constants::PION_PLUS_PID  ){
            zMC_values.push_back(MChadron.GetzMC());
            pt2MC_values.push_back(MChadron.Getpt2MC());
            phihMC_values.push_back(MChadron.GetphihMC());
        }
    }

    // select the smallest vector to fill both histo with max those n values
    for (size_t i = 0; i < std::min(z_values.size(), zMC_values.size()); ++i) {
        h_zcomp->Fill(z_values[i], zMC_values[i]);
        h_z->Fill(z_values[i]);
        h_zMC->Fill(zMC_values[i]);
    }
    for (size_t i = 0; i < std::min(pt2_values.size(), pt2MC_values.size()); ++i) {
        h_pt2comp->Fill(pt2_values[i], pt2MC_values[i]);
        h_pt2->Fill(pt2_values[i]);
        h_pt2MC->Fill(pt2MC_values[i]);
    }
    for (size_t i = 0; i < std::min(phih_values.size(), phihMC_values.size()); ++i) {
        h_phihcomp->Fill(phih_values[i], phihMC_values[i]);
        h_phih->Fill(phih_values[i]);
        h_phihMC->Fill(phihMC_values[i]);
    }


}  

    

void Monunfold::DrawHistoRec(const std::string filename) {
    // Draw histograms and save them to a PDF file and in a root file 
    TFile* rootFile = new TFile((filename + ".root").c_str(), "RECREATE");
    TCanvas cREC("cREC", "cREC");
    cREC.Divide(3, 3);
    cREC.cd(1);
    h_Q2->Draw();
    h_Q2->Write();
    cREC.cd(2);
    h_xb->Draw();
    h_xb->Write();
    cREC.cd(3);
    h_y->Draw();
    h_y->Write();
    cREC.cd(4);
    h_nu->Draw();
    h_nu->Write();
    cREC.cd(5);
    h_W2->Draw();
    h_W2->Write();
    cREC.cd(6);
    h_z->Draw();
    h_z->Write();
    cREC.cd(7);
    h_pt2->Draw();
    h_pt2->Write();
    cREC.cd(8);
    h_phih->Draw();
    h_phih->Write();
    cREC.cd(9);
    h_vertexZ->Draw();
    h_vertexZ->Write();
    cREC.Print((filename + ".pdf").c_str());
    rootFile->Close();
    delete rootFile; 
}


void Monunfold::SaveHistRoot(const std::string& filenameREC) {
    TFile* rootFile = new TFile((filenameREC + ".root").c_str(), "RECREATE");
    //outputFile.Open(filenameREC.c_str(), "RECREATE");

    // Write histograms related to real data
    h_vertexZ->Write();
    h_Q2->Write();
    h_xb->Write();
    h_y->Write();
    h_nu->Write();
    h_W2->Write();
    //h_px_el->Write();
    //h_py_el->Write();
    //h_pz_el->Write();
    //h_ptot_el->Write();
    //h_theta_el->Write();
    //h_phi_el->Write();
    //h_polcoord_el->Write();
    //h_E_el->Write();
    //h_E_el_theta->Write();
    //h_E_el_phi->Write();
    h_z->Write();
    h_pt2->Write();
    h_phih->Write();
    //h_px_pi->Write();
    //h_py_pi->Write();
    //h_pz_pi->Write();
    //h_ptot_pi->Write();
    //h_theta_pi->Write();
    //h_phi_pi->Write();
    //h_polcoord_pi->Write();
    //h_E_pi->Write();
    //h_E_pi_theta->Write();
    //h_E_pi_phi->Write();
    //h_vertexZ_pi->Write();
    //h_DeltaVz->Write();

    //    outputFile.Close();
    h_chi2_el->Write();
    h_chi2_pi->Write();
    rootFile->Close();
    delete rootFile; 


}


void Monunfold::DrawHistoMC(const std::string filename) {
    // Draw histograms and save them to a PDF file
    TCanvas cMC("cMC", "cMC");
    cMC.Divide(3, 3);
    cMC.cd(1);
    h_Q2MC->SetLineColor(kRed);
    h_Q2MC->Draw();
    cMC.cd(2);
    h_xbMC->SetLineColor(kRed);
    h_xbMC->Draw();
    cMC.cd(3);
    h_yMC->SetLineColor(kRed);
    h_yMC->Draw();
    cMC.cd(4);
    h_nuMC->SetLineColor(kRed);
    h_nuMC->Draw();
    cMC.cd(5);
    h_W2MC->SetLineColor(kRed);
    h_W2MC->Draw();
    cMC.cd(6);
    h_zMC->SetLineColor(kRed);
    h_zMC->Draw();
    cMC.cd(7);
    h_pt2MC->SetLineColor(kRed);
    h_pt2MC->Draw();
    cMC.cd(8);
    h_phihMC->SetLineColor(kRed);
    h_phihMC->Draw();
    cMC.cd(9);
    h_vertexZMC->SetLineColor(kRed);
    h_vertexZMC->Draw();
    cMC.Print((filename + ".pdf").c_str());

}


void Monunfold::SaveHistMCRoot(const std::string& filenameMC) {
    TFile* rootFile = new TFile((filenameMC + ".root").c_str(), "RECREATE");

    //outputFile.Open(filenameMC.c_str(), "UPDATE"); // "UPDATE" mode to append to an existing file or create a new one

    // Write histograms related to MC simulated data
    h_vertexZMC->Write();
    h_Q2MC->Write();
    h_xbMC->Write();
    h_yMC->Write();
    h_nuMC->Write();
    h_W2MC->Write();
    //h_px_el->Write(); // Assuming these are the same as real data for simplicity
    //h_py_el->Write();
    //h_pz_el->Write();
    //h_ptot_el->Write();
    //h_theta_el->Write();
    //h_phi_el->Write();
    //h_polcoord_el->Write();
    //h_E_el->Write();
    //h_E_el_theta->Write();
    //h_E_el_phi->Write();
    h_zMC->Write();
    h_pt2MC->Write();
    h_phihMC->Write();
    //h_px_pi->Write(); // Assuming these are the same as real data for simplicity
    //h_py_pi->Write();
    //h_pz_pi->Write();
    //h_ptot_pi->Write();
    //h_theta_pi->Write();
    //h_phi_pi->Write();
    //h_polcoord_pi->Write();
    //h_E_pi->Write();
    //h_E_pi_theta->Write();
    //h_E_pi_phi->Write();
    //h_vertexZMC->Write(); // Assuming this is the MC hadron vertex Z histogram
    //h_DeltaVz->Write();

    //outputFile.Close();
        rootFile->Close();
    delete rootFile; 

}


void Monunfold::DrawHistograms(const std::string filename) {
    // Draw histograms and save them to a PDF file
    TCanvas c1("c1", "c1", 800, 600);
    c1.Divide(3, 3);
    c1.cd(1);
    //h_Q2->Draw();
    h_Q2MC->SetLineColor(kRed);
    h_Q2MC->Draw();
    c1.cd(2);
    //h_xb->Draw();
    h_xbMC->SetLineColor(kRed);
    h_xbMC->Draw();
    c1.cd(3);
    //h_y->Draw();
    h_yMC->SetLineColor(kRed);
    h_yMC->Draw();
    c1.cd(4);
    //h_nu->Draw();
    h_nuMC->SetLineColor(kRed);
    h_nuMC->Draw();
    c1.cd(5);
    //h_W2->Draw();
    h_W2MC->SetLineColor(kRed);
    h_W2MC->Draw();
    c1.cd(6);
    //h_z->Draw();
    h_zMC->SetLineColor(kRed);
    h_zMC->Draw();

    c1.cd(7);
    //h_pt2->Draw();
    h_pt2MC->SetLineColor(kRed);
    h_pt2MC->Draw();
    c1.cd(8);
    //h_phih->Draw();
    h_phihMC->SetLineColor(kRed);
    h_phihMC->Draw();
    c1.cd(9);
    //h_vertexZ->Draw();
    h_vertexZMC->SetLineColor(kRed);
    h_vertexZMC->Draw();
    c1.Print((filename + ".pdf").c_str());

}
//void Monunfold::DrawCompRECMC(const std::string& filename) {
//    TCanvas c2("c2", "c2", 800, 600);
//    c2.Divide(3, 3);
//    
//    // Define X and Y axis labels for each histogram with units
//    std::vector<std::pair<std::string, std::string>> axisLabels = {
//        {"Q^{2} REC (GeV^{2})", "Q^{2} MC "},        // Q2
//        {"x_{B} REC", "x_{B} MC"},          // xb
//        {"y REC", "y MC"},                 // y
//        {"\\nu REC (GeV)", "\\nu MC"},     // nu
//        {"W^{2} REC (GeV^{2})", "W^{2} MC"}, // W2
//        {"z REC", "z MC"},                 // z
//        {"p_{T}^{2} REC (GeV^{2})", "p_{T}^{2} MC "}, // pt2
//        {"\\phi_{h} REC (deg)", "\\phi_{h} MC"},   // phih
//        {"Vertex Z REC (cm)", "Vertex Z MC "} // vertexZ
//    };
//
//    std::vector<TH2*> histograms = {
//        h_Q2comp,
//        h_xbcomp,
//        h_ycomp,
//        h_nucomp,
//        h_W2comp,
//        h_zcomp,
//        h_pt2comp,
//        h_phihcomp,
//        h_vertexZcomp
//    };
//    
//    for (int i = 0; i < 4; ++i) {
//        c2.cd(i + 1);
//        histograms[i]->SetXTitle(axisLabels[i].first.c_str());
//        histograms[i]->SetYTitle(axisLabels[i].second.c_str());
//        histograms[i]->Draw("colz");
//    }
//
//    c2.Print((filename + ".pdf").c_str());
//    
//    TFile rootFile((filename + ".root").c_str(), "RECREATE");
//    for (auto hist : histograms) {
//        hist->Write();
//    }
//    rootFile.Close();
//}
//


void Monunfold::DrawCompRECMC(const std::string& filename) {
// Draw histograms and save them to a PDF file
    TCanvas c1("c3", "c3", 800, 600);
    c1.Divide(3, 3);
    c1.cd(1);
    //h_Q2->Draw();
    h_Q2comp->Draw("colz");
    c1.cd(2);
    //h_xb->Draw();
    h_xbcomp->Draw("colz");
    c1.cd(3);
    //h_y->Draw();
    h_ycomp->Draw("colz");
    c1.cd(4);
    //h_evtnbrdiff->Draw();
    //h_nu->Draw();
    h_nucomp->Draw("colz");
    c1.cd(5);
    //h_W2->Draw();
    h_W2comp->Draw("colz");
    c1.cd(6);
    //h_z->Draw();  
    h_zcomp->Draw("colz");
    c1.cd(7);
    //h_pt2->Draw();
    h_pt2comp->Draw("colz");
    c1.cd(8);
    //h_phih->Draw();
    h_phihcomp->Draw("colz");
    c1.cd(9);
    //h_vertexZ->Draw();
    h_vertexZcomp->Draw("colz");

    //c1.Print((filename + ".pdf").c_str());
    c1.Print((filename + ".pdf(").c_str());
TCanvas c2("c2", "c2", 800, 600);
    c2.Divide(3, 3);
    c2.cd(1);
    h_px_elcomp->Draw("colz");
    c2.cd(2);
    h_py_elcomp->Draw("colz");
    c2.cd(3);
    h_pz_elcomp->Draw("colz");
    c2.cd(4);
    h_ptot_elcomp->Draw("colz");
    c2.cd(5);
    h_theta_elcomp->Draw("colz");
    c2.cd(6);
    h_phi_elcomp->Draw("colz");
    c2.cd(7);
    h_px_picomp->Draw("colz");
    c2.cd(8);
    h_py_picomp->Draw("colz");
    c2.cd(9);
    h_pz_picomp->Draw("colz");
    c2.Print((filename + ".pdf").c_str());
TCanvas c3("c3", "c3", 800, 600);
    c3.Divide(3, 3);
    c3.cd(1);
    h_ptot_picomp->Draw("colz");
    c3.cd(2);
    h_theta_picomp->Draw("colz");
    c3.cd(3);
    h_phi_picomp->Draw("colz");
    c3.cd(4);
    h_E_elcomp->Draw("colz");
    c3.cd(5);
    h_E_picomp->Draw("colz");
    c3.cd(6);
    h_phihcomp->Draw("colz");
    c3.cd(7);
    h_pt2zMC->Draw("colz");


    c3.Print((filename + ".pdf").c_str());
TCanvas c4("c4", "c4", 800, 600);
    c4.Divide(3, 3);
    c4.cd(1);
    h_px_piDelta->Draw();
    c4.cd(2);
    h_py_piDelta->Draw();
    c4.cd(3);
    h_pz_piDelta->Draw();
    c4.cd(4);
    h_ptot_piDelta->Draw();
    c4.cd(5);
    h_theta_piDelta->Draw();
    c4.cd(6);
    h_phi_piDelta->Draw();
    c4.cd(7);
    h_E_piDelta->Draw();
    c4.Print((filename + ".pdf)").c_str());

TCanvas c5("c5", "c5", 800, 600);
    c5.Divide(3, 3);
    c5.cd(1);
    h_Delta_px_el->Draw();
    c5.cd(2);
    h_Delta_py_el->Draw();
    c5.cd(3);
    h_Delta_pz_el->Draw();
    c5.cd(4);
    h_Delta_ptot_el->Draw();
    c5.cd(5);
    h_Delta_theta_el->Draw();
    c5.cd(6);
    h_Delta_phi_el->Draw();
    c5.cd(7);
    h_E_elcomp->Draw("colz");   
    c5.Print((filename + ".pdf").c_str());




    


    //c2.Print((filename + ".pdf)").c_str());



}

void Monunfold::DrawMomentainSim(const std::string& filename){
    TCanvas c1("c4", "c4", 800, 600);
    c1.Divide(4, 4);
    c1.cd(1);
    h_px_elMC->SetLineColor(kRed);
    h_px_elMC->Draw();
    h_px_el->Draw("same");

    c1.cd(2);
    h_py_elMC->SetLineColor(kRed);
    h_py_elMC->Draw();
    h_py_el->Draw("same");
    c1.cd(3);
    h_pz_elMC->SetLineColor(kRed);
    h_pz_elMC->Draw();
    h_pz_el->Draw("same");
    c1.cd(4);
    h_ptot_elMC->SetLineColor(kRed);
    h_ptot_elMC->Draw();
    h_ptot_el->Draw("same");
    c1.cd(5);
    h_theta_elMC->SetLineColor(kRed);
    h_theta_elMC->Draw();
    h_theta_el->Draw("same");
    c1.cd(6);
    h_phi_elMC->SetLineColor(kRed);
    h_phi_elMC->Draw();
    h_phi_el->Draw("same");
    c1.cd(7);
    h_px_piMC->SetLineColor(kRed);
    h_px_piMC->Draw();
    h_px_pi->Draw("same");
    c1.cd(8);
    h_py_piMC->SetLineColor(kRed);
    h_py_piMC->Draw();
    h_py_pi->Draw("same");
    c1.cd(9);
    h_pz_piMC->SetLineColor(kRed);
    h_pz_piMC->Draw();
    h_pz_pi->Draw("same");
    c1.cd(10);
    h_ptot_piMC->SetLineColor(kRed);
    h_ptot_piMC->Draw();
    h_ptot_pi->Draw("same");
    c1.cd(11);
    h_theta_piMC->SetLineColor(kRed);
    h_theta_piMC->Draw();
    h_theta_pi->Draw("same");
    c1.cd(12);
//    h_phi_pi->Draw();
//    h_phi_piMC->SetLineColor(kRed);
//    h_phi_piMC->Draw("same");

    h_phihMC->SetLineColor(kRed);
    h_phihMC->Draw();
    h_phih->Draw("same");
    c1.cd(13);
    h_E_elMC->SetLineColor(kRed);
    h_E_elMC->Draw();
    h_E_el->Draw("same");
    c1.cd(14);
    h_E_piMC->SetLineColor(kRed);
    h_E_piMC->Draw();
    h_E_pi->Draw("same");
    c1.cd(15);
    h_deltathetapi->Draw();
    c1.cd(16);
    h_deltaphipi->Draw();





    c1.Print((filename + ".pdf").c_str());
}





Monunfold::~Monunfold() {
    // Delete dynamically allocated histograms
    delete h_Q2;
    delete h_xb;
    delete h_y;
    delete h_nu;
    delete h_W2;
    //delete h_z;
    //delete h_pt2;
    //delete h_phih;
    //delete h_vertexZ;
    //delete h_pid;
    
}
