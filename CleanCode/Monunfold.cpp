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



    h_evtnbrdiff(new TH1F(("U_evtnbrdiff_" + targetName).c_str(), "evtnbrdiff", 100, -10, 10.0)),




      counterel_R(0) {
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
    if (targetType == 0 && cut1.PassCutsElectrons(event)==true && cut1.PassCutsDetectors(event)==true ) { 
        if (event.Getxb()>Constants::LargexblowAB && event.Getxb()<Constants::LargexbhighAB){   //for regions A and B 
            if (event.GetQ2()>Constants::LARGEQlowA && event.GetQ2()<Constants::LARGEQhighA) {      //REGION A
                counterLargeD_regA ++;
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
    if (cut1.PassCutsDetectors(event)) {
        // Fill histograms after passing detector cuts for MC data
        h_vertexZMC->Fill(event.GetVz()); // Fill MC vertex Z histogram

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
