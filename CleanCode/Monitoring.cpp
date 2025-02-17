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
#include "Monitoring.h"
#include "CutSet.h"
#include "constants.h"
#include <TStyle.h>


    //int nubin = 100;  

Monitoring::Monitoring(CutSet a, const std::string& targetName)
    : cut1(a), targetName(targetName),
    R_nu_el(new TH1F (("R_nuel_" + targetName).c_str(), "R_nuel", Rbin , numinR, numaxR)),
    R_nu_had(new TH1F (("R_nuhad_" + targetName).c_str(), "R_nuhad", Rbin , numinR, numaxR)),
    R_z(new TH1F (("R_z_" + targetName).c_str(), "R_z", Rbin , zminR, zmaxR)),
    R_pt2(new TH1F (("R_pt2_" + targetName).c_str(), "R_pt2", Rbin , pt2minR, pt2maxR)),

    h_Q2(new TH1F(("Q2_" + targetName).c_str(), "Q2", nubin, QminX, QmaxX)),
    h_xb(new TH1F(("xb_" + targetName).c_str(), "xb", nubin, xminX, xmaxX)),
    h_y(new TH1F(("y_" + targetName).c_str(), "y", nubin, yminX, ymaxX)),
    h_nu(new TH1F(("nu_" + targetName).c_str(), "nu", nubin, numinX, numaxX)),
    h_W2(new TH1F(("W2_" + targetName).c_str(), "W2", nubin, WminX, WmaxX)),
    h_z(new TH1F(("z_" + targetName).c_str(), "z", nubin, zminX, zmaxX)),
    h_pt2(new TH1F(("pt2_" + targetName).c_str(), "pt2", nubin, pt2minX, pt2maxX)),
    h_phih(new TH1F(("phih_" + targetName).c_str(), "phih", nubin, phihminX, phihmaxX)),
    h_vertexZ(new TH1F(("targetVz_" + targetName).c_str(), "vertexz4target", 100, -10, 0.0)),
    h_vertexY(new TH1F(("targetVy_" + targetName).c_str(), "vertexy4target", 100, -10, 0.0)),
    h_vertexX(new TH1F(("targetVx_" + targetName).c_str(), "vertexx4target", 100, -10, 0.0)),
    h_vertexZ_pi(new TH1F(("targetVz_pi_" + targetName).c_str(), "vertex4target", 100, -10, 0.0)),
    h_DeltaVz(new TH1F(("DeltaVz_" + targetName).c_str(), "DeltaVz", 100, -5, 5)),
    
    h_Q2MC(new TH1F(("Q2_MC" + targetName).c_str(), "Q2MC", nubin, QminX, QmaxX)),
    h_xbMC(new TH1F(("xb_MC" + targetName).c_str(), "xbMC", nubin, xminX, xmaxX)),
    h_yMC(new TH1F(("y_MC" + targetName).c_str(), "yMC", nubin, yminX, ymaxX)),
    h_nuMC(new TH1F(("nu_MC" + targetName).c_str(), "nuMC", nubin, numinX, numaxX)),
    h_W2MC(new TH1F(("W2_MC" + targetName).c_str(), "W2MC", nubin, WminX, WmaxX)),
    h_zMC(new TH1F(("z_MC" + targetName).c_str(), "zMC", nubin, zminX, zmaxX)),
    h_pt2MC(new TH1F(("pt2_MC" + targetName).c_str(), "pt2MC", nubin, pt2minX, pt2maxX)),
    h_phihMC(new TH1F(("phih_MC" + targetName).c_str(), "phihMC", nubin, phihminX, phihmaxX)),
    h_vertexZMC(new TH1F(("targetVz_MC" + targetName).c_str(), "vertex4ztargetMC", 100, -10, 0.0)),
    h_vertexYMC(new TH1F(("targetVy_MC" + targetName).c_str(), "vertexy4targetMC", 100, -10, 0.0)),
    h_vertexXMC(new TH1F(("targetVx_MC" + targetName).c_str(), "vertexx4targetMC", 100, -10, 0.0)),
    h_thetaelectron (new TH1F(( "thetaelectron_" + targetName).c_str(), "thetaelectron", 500, 0, 1)),
    h_rapport (new TH1F(( " rapport_beam_el" + targetName).c_str(), "rapport_beam_el", 500, 0, 1)),
    h_thetaelectronMC (new TH1F(( "thetaelectron_MC" + targetName).c_str(), "thetaelectronMC", 500, 0, 1)),
    h_rapportMC (new TH1F(( " rapport_beam_elMC" + targetName).c_str(), "rapport_beam_elMC", 500, 0, 1)),



    h_pid(new TH1F(("pid_" + targetName).c_str(), "pid", 100, -250, 250)),
    h_xQ2(new TH2F(("xQ2_" + targetName).c_str(), "xQ2", nubin, xminX, xmaxX, nubin, QminX, QmaxX)),
    h_xQ2pos(new TH2F(("xQ2pos_" + targetName).c_str(), "xQ2pos", nubin, xminX, xmaxX, nubin, QminX, QmaxX)),
    h_pt2z(new TH2F(("pt2z_" + targetName).c_str(), "pt2z", nubin, pt2minX, pt2maxX, nubin, zminX, zmaxX)),
    h_px_el(new TH1F(("px_ele_" + targetName).c_str(), "px_ele", nubin, 0, 10)),
    h_py_el(new TH1F(("py_ele_" + targetName).c_str(), "py_ele", nubin, 0, 10)),
    h_pz_el(new TH1F(("pz_ele_" + targetName).c_str(), "pz_ele", nubin, 0, 10)),
    h_ptot_el(new TH1F(("ptot_ele_" + targetName).c_str(), "ptot_ele", nubin, 0, 10)),
    h_px_pi(new TH1F(("px_pi_" + targetName).c_str(), "px_pi", nubin, 0, 10)),
    h_py_pi(new TH1F(("py_pi_" + targetName).c_str(), "py_pi", nubin, 0, 10)),
    h_pz_pi(new TH1F(("pz_pi_" + targetName).c_str(), "pz_pi", nubin, 0, 10)),
    h_ptot_pi(new TH1F(("ptot_pro_" + targetName).c_str(), "ptot_pro", nubin, 0, 10)),
    h_theta_el(new TH1F(("theta_el" + targetName).c_str(), "theta", nubin, 0, 30)),
    h_phi_el(new TH1F(("phi_el" + targetName).c_str(), "phi", nubin, 0, 360)),
    h_polcoord_el(new TH2F(("pol_el"+ targetName).c_str(), "polcoord", nubin, 0,40, nubin, 0, 360)),
    h_theta_pi(new TH1F(("theta_pi" + targetName).c_str(), "theta", nubin, 0, 150)),
    h_phi_pi(new TH1F(("phi_pi" + targetName).c_str(), "phi", nubin, 0, 360)),
    h_polcoord_pi(new TH2F(("pol_pi"+ targetName).c_str(), "polcoordphi", nubin, 0,150, nubin, 0, 360)),
    h_E_el(new TH1F(("E_el" + targetName).c_str(), "E", nubin, 0, 10)),
    h_E_el_theta(new TH2F(("E_el_theta" + targetName).c_str(), "E_theta", nubin, 0, 30, nubin, 0, 10)),
    h_E_el_phi(new TH2F(("E_el_phi" + targetName).c_str(), "E_phi", nubin, 0, 360, nubin, 0, 10)),
    h_E_pi(new TH1F(("E_pi" + targetName).c_str(), "E", nubin, 0, 10)),
    h_E_pi_theta(new TH2F(("E_pi_theta" + targetName).c_str(), "E_theta", nubin, 0, 150, nubin, 0, 10)),
    h_E_pi_phi(new TH2F(("E_pi_phi" + targetName).c_str(), "E_phi", nubin, 0, 360, nubin, 0, 10)),
    h_lu(new TH1F(("lu_el" + targetName).c_str(), "lu", nubin, 0, 450)),
    h_lv(new TH1F(("lv_el" + targetName).c_str(), "lv", nubin, 0, 350)),
    h_lw(new TH1F(("lw_el" + targetName).c_str(), "lw", nubin, 0, 350)),
    h_epcal(new TH1F(("epcal_el" + targetName).c_str(), "epcal", nubin, 0, 1)),
    h_eecalin(new TH1F(("eecalin_el" + targetName).c_str(), "eecalin", nubin, 0, 1)),
    h_epcalout(new TH1F(("epcalout_el" + targetName).c_str(), "epcalout", nubin, 0, 1)),
    h_calXY(new TH2F(("calxy_el" + targetName).c_str(), "calxy", nubin, -400, 400, nubin, -400, 400)),
    h_calEall(new TH2F(("calEall_el" + targetName).c_str(), "calEall", nubin, 0, 1, nubin, 0, 1)),
    h_Nphe15(new TH1F(("Nphe15_" + targetName).c_str(), "Nphe15", nubin, 0, 60)),
    h_Nphe16(new TH1F(("Nphe16_" + targetName).c_str(), "Nphe16", nubin, 0, 60)),
    h_calSector(new TH1F(("calSector_" + targetName).c_str(), "calSector", nubin, -1, 10)),
    h_helicity(new TH1F(("helicity_" + targetName).c_str(), "helicity", nubin, -2, 2)),
    h_helicity_raw(new TH1F(("helicity_raw_" + targetName).c_str(), "helicity_raw", nubin, -2, 2)),
    h_phih_plus(new TH1F(("phih_plus_" + targetName).c_str(), "phih_plus", 10, 0, 360)),
    h_phih_minus(new TH1F(("phih_minus_" + targetName).c_str(), "phih_minus", 10, 0, 360)),
    h_BSA(new TH1F(("BSA_" + targetName).c_str(), "BSA", 10, 0, 360)),
    h_Q2comp(new TH2F(("Q2comp_" + targetName).c_str(), "Q2comp", nubin, QminX, QmaxX,nubin, QminX, QmaxX)),
    h_chi2_el(new TH1F(("chi2_el_" + targetName).c_str(), "chi2_el", 100, -10, 10)),
    h_chi2_pi(new TH1F(("chi2_pi_" + targetName).c_str(), "chi2_pi", 100, -10, 10)),
    h_chi2_pid_pi(new TH2F(("chi2_pid_pi_" + targetName).c_str(), "chi2_pid_pi", 100, -15, 15, 100, 0, 1)),
    h_luthetael(new TH2F(("u_vs_thetael_" + targetName).c_str(), "u_vs_thetael", 100, 0, 400, 100, 0, 40)),
    h_sampl_el(new TH2F(("sampl_el_" + targetName).c_str(), "sampl_el", nubin, 0, 10, nubin, 0, 1)),
      counterel_R(0) {
    // Add more histograms as needed
}


//void Monitoring::FillQ2pre(const Event& event ){
//    h_Q2pre->Fill(event.GetQ2());
//}
//void Monitoring::Fillypre(const Event& event ){
//    h_ypre->Fill(event.Gety());
//}
//void Monitoring::Fillnupre(const Event& event ){
//    h_xbpre->Fill(event.Getxb());
//}
//void Monitoring::Fillnupre(const Event& event ){
//    h_nupre->Fill(event.Getnu());
//}
//void Monitoring::FillW2pre(const Event& event ){
//    h_W2pre->Fill(event.GetW2());
//}



void Monitoring::FillHistogramswCuts(const Event& event) {              /// good CUTS in hadron !!!!!!!
    //if (cut1.PassCutsDetectors(event)==true){
        // Fill Detector histograms after cuts



        //h_eecalin->Fill(event.GetEcalin());
        if (event.electron.GetEcalin()>0.01){h_eecalin->Fill(event.electron.GetEcalin());}
        //h_epcalout->Fill(event.GetEcalout());
        if (event.electron.GetEcalout()>0.01){h_epcalout->Fill(event.electron.GetEcalout());}
        //h_calEall->Fill(event.GetEpcal(), event.GetEcalin()+event.GetEcalout());
        if (event.electron.GetEcalout()>0.01 & event.electron.GetEcalout()>0.01){h_calEall->Fill(event.electron.GetEpcal(), event.electron.GetEcalin()+event.electron.GetEcalout());}
        h_Nphe15->Fill(event.electron.Getnphe15());
        h_Nphe16->Fill(event.electron.Getnphe16());
        h_calSector->Fill(event.electron.GetCalSector());
        h_helicity->Fill(event.GetHel());
        h_helicity_raw->Fill(event.GetHelRaw());
        if (cut1.PassCutsElectrons(event)==true) {
            // Fill Electron variable histograms after cuts on electron
            h_vertexZ->Fill(event.GetVz());     //Vz only exists when an electron is detected !!!!
                                                //add Vz for the hadron too
            h_Q2->Fill(event.electron.GetQ2());

            h_xb->Fill(event.Getxb());
            h_xQ2->Fill(event.Getxb(), event.electron.GetQ2());
            h_y->Fill(event.Gety());    
            h_nu->Fill(event.Getnu());
            h_W2->Fill(event.GetW2());
            h_chi2_el->Fill(event.electron.Getchi2());

            h_px_el->Fill(event.electron.GetMomentum().X());
            h_py_el->Fill(event.electron.GetMomentum().Y());
            h_pz_el->Fill(event.electron.GetMomentum().Z());
            h_ptot_el->Fill(event.electron.GetMomentum().P());

            h_sampl_el->Fill(event.electron.GetMomentum().P(), (event.electron.GetEpcal()/event.electron.GetMomentum().P()) );
            h_theta_el->Fill(event.electron.GetMomentum().Theta()*180/Constants::PI);
            h_phi_el->Fill(event.electron.GetMomentum().Phi()*180/Constants::PI +180);
            h_polcoord_el->Fill(event.electron.GetMomentum().Theta()*180/Constants::PI, event.electron.GetMomentum().Phi()*180/Constants::PI +180);
            h_E_el_theta->Fill(event.electron.GetMomentum().Theta()*180/Constants::PI, event.electron.GetMomentum().E());
            h_E_el_phi->Fill(event.electron.GetMomentum().Phi()*180/Constants::PI +180, event.electron.GetMomentum().E());
        h_calXY->Fill(event.electron.GetCalX(), event.electron.GetCalY());
        h_lu->Fill(event.electron.Getlu());
        h_lv->Fill(event.electron.Getlv());
        h_lw->Fill(event.electron.Getlw());
            h_luthetael->Fill(event.electron.Getlu(), event.electron.GetMomentum().Theta()*180/Constants::PI);

            for (const Particle& hadron : event.GetHadrons()) {
                if (cut1.PassCutsHadrons(hadron)==true){

                    if (hadron.GetPID() == Constants::PION_PLUS_PID  ){   //adding this condition for pion+ and erasing the condit at evtprocessr

            h_E_el->Fill(event.electron.GetMomentum().E());

                        h_epcal->Fill(event.electron.GetEpcal());
                        h_chi2_pi->Fill(hadron.Getchi2());
                        h_z->Fill(hadron.Getz());
                        h_pt2->Fill(hadron.Getpt2());
                        h_pt2z->Fill(hadron.Getpt2(), hadron.Getz());
                        h_phih->Fill(hadron.Getphih());
                        h_px_pi->Fill(hadron.GetMomentum().X());
                        h_py_pi->Fill(hadron.GetMomentum().Y());
                        h_pz_pi->Fill(hadron.GetMomentum().Z());
                        h_ptot_pi->Fill(hadron.GetMomentum().P());
                        h_theta_pi->Fill(hadron.GetMomentum().Theta()*180/Constants::PI);
                        h_phi_pi->Fill(hadron.GetMomentum().Phi()*180/Constants::PI +180);
                        h_polcoord_pi->Fill(hadron.GetMomentum().Theta()*180/Constants::PI, hadron.GetMomentum().Phi()*180/Constants::PI +180);
                        h_E_pi->Fill(hadron.GetMomentum().E());
                        h_E_pi_theta->Fill(hadron.GetMomentum().Theta()*180/Constants::PI, hadron.GetMomentum().E());
                        h_E_pi_phi->Fill(hadron.GetMomentum().Phi()*180/Constants::PI +180, hadron.GetMomentum().E());
                        h_vertexZ_pi->Fill(hadron.GetParticleVertexZ());
                        h_DeltaVz->Fill(event.GetVz()-hadron.GetParticleVertexZ());
                        h_helicity->Fill(event.GetHel());
                        h_helicity_raw->Fill(event.GetHelRaw());
                        h_chi2_pid_pi->Fill(hadron.Getchi2(),  (event.electron.GetEpcal()/event.electron.GetMomentum().P()) );

                    }
                }


            }
        }
    //}
}

/*
void Monitoring::FillHistograms(const Event& event) {
    //if (cut1.PassCutLD2Target(event)==false) return;
    h_xQ2->Fill(event.Getxb(), event.GetQ2());
    h_calXY->Fill(event.GetCalX(), event.GetCalY());
    h_lu->Fill(event.Getlu());
    h_lv->Fill(event.Getlv());
    h_lw->Fill(event.Getlw());
    h_epcal->Fill(event.GetEpcal());
        //h_eecalin->Fill(event.GetEcalin());
    if (event.GetEcalin()>0.5){h_eecalin->Fill(event.GetEcalin());}
    //h_epcalout->Fill(event.GetEcalout());
    if (event.GetEcalout()>0.02){h_epcalout->Fill(event.GetEcalout());}
    if (event.GetEcalout()>0.02 & event.GetEcalout()>0.02){h_calEall->Fill(event.GetEpcal(), event.GetEcalin()+event.GetEcalout());}
    //h_calEall->Fill(event.GetEpcal(), event.GetEcalin()+event.GetEcalout());
    h_Nphe15->Fill(event.Getnphe15());
    h_Nphe16->Fill(event.Getnphe16());
    h_calSector->Fill(event.GetCalSector());

    if (cut1.PassCutsElectrons(event)==false) return;
    //std::cout << " -> electron cuts passed successfully  " << std::endl;
    counterel_R ++;
    h_vertexZ->Fill(event.GetVz());     //Vz only exists when an electron is detected !!!!
    h_Q2->Fill(event.GetQ2());
    h_xb->Fill(event.Getxb());
    h_y->Fill(event.Gety());    
    h_nu->Fill(event.Getnu());
    h_W2->Fill(event.GetW2());

    

    //h_Q2posel->Fill(event.GetQ2());
    //h_nuposel->Fill(event.Getnu());     //throw hists in mratio.cpp(filling) TBD
    //std::cout << " kinel " << event.GetQ2()<<" , " << event.Getxb()<<" , " << event.Gety()<<" , "<< event.Getnu()<<" , "<< event.GetW2()<<std::endl;  
    for (const Particle& hadron : event.GetHadrons()) {
            //         h_z->Fill(hadron.Getz());
            // h_pt2->Fill(hadron.Getpt2());

        if (cut1.PassCutsHadrons(hadron)==true){
            
            if (hadron.GetPID() == Constants::PION_PLUS_PID  ){   //adding this condition for pion+ and erasing the condit at evtprocessr
                                                                // PID condition added for CLAS COLL        //DELETE THIS 
                                                                //should be replaced in CUTSET

            //std::cout << " -> hadron cuts passed successfully  " << std::endl;
            //if passcuts(given hadron )
            //fill histogram with hadron variables 
            h_pid->Fill(hadron.GetPID());
            //h_Q2_had->Fill(event.GetQ2());
            //h_nu_had->Fill(event.Getnu());
            //h_z_had->Fill(hadron.Getz());
            //h_pt2_had->Fill(hadron.Getpt2());
            h_z->Fill(hadron.Getz());
            h_pt2->Fill(hadron.Getpt2());
            h_phih->Fill(hadron.Getphih());
            h_xQ2pos->Fill(event.Getxb(), event.GetQ2());

            //
            //std::cout << " kinhad " << hadron.Getz()<<" , " << hadron.Getpt2()<<" , " << hadron.Getphih()<<std::endl;  
            //
            }
        }
    }
 //add histos 

}


void Monitoring::FillHistogramsNoCutsMC(const Event& event) {
    h_Q2MC->Fill(event.GetQ2MC());
    h_Q2comp->Fill(event.GetQ2(), event.GetQ2MC());
    h_Q2->Fill(event.GetQ2());
    h_xbMC->Fill(event.GetxbMC());
    h_yMC->Fill(event.GetyMC());
    h_nuMC->Fill(event.GetnuMC());
    h_W2MC->Fill(event.GetW2MC());
    h_vertexZMC->Fill(event.GetVzMC());
    h_thetaelectronMC->Fill(event.GetThetaElectronMC());
    h_rapportMC->Fill(event.GetAcosyadaMC());
    for (const Particle& MChadron : event.GetMCHadrons()) {
        if (MChadron.GetPID() == Constants::PION_PLUS_PID  ){   //adding this condition for pion+ and erasing the condit at evtprocessr
        h_zMC->Fill(MChadron.GetzMC());
        h_pt2MC->Fill(MChadron.Getpt2MC());
        h_phihMC->Fill(MChadron.GetphihMC());
        
        
        //
        //std::cout << " kinhad " << hadron.Getz()<<" , " << hadron.Getpt2()<<" , " << hadron.Getphih()<<std::endl;  
        //
        }
    }



}
*/
void Monitoring::FillHistogramsNoCuts(const Event& event) {
    // Fill histograms with no cuts
    h_calXY->Fill(event.electron.GetCalX(), event.electron.GetCalY());
    h_lu->Fill(event.electron.Getlu());
    h_lv->Fill(event.electron.Getlv());
    h_lw->Fill(event.electron.Getlw());
    h_epcal->Fill(event.electron.GetEpcal());
    //h_eecalin->Fill(event.GetEcalin());
    //h_epcalout->Fill(event.GetEcalout());
    //h_calEall->Fill(event.GetEpcal(), event.GetEcalin()+event.GetEcalout());
            //h_eecalin->Fill(event.GetEcalin());
    if (event.electron.GetEcalin()>0.01){h_eecalin->Fill(event.electron.GetEcalin());}
    //h_epcalout->Fill(event.GetEcalout());
    if (event.electron.GetEcalout()>0.01){h_epcalout->Fill(event.electron.GetEcalout());}
    if (event.electron.GetEcalout()>0.01 & event.electron.GetEcalout()>0.01){h_calEall->Fill(event.electron.GetEpcal(), event.electron.GetEcalin()+event.electron.GetEcalout());}
    //h_calEall->Fill(event.GetEpcal(), event.GetEcalin()+event.GetEcalout());
    h_Nphe15->Fill(event.electron.Getnphe15());
    h_Nphe16->Fill(event.electron.Getnphe16());
    h_calSector->Fill(event.electron.GetCalSector());
    h_helicity->Fill(event.GetHel());
    h_helicity_raw->Fill(event.GetHelRaw());


    h_Q2->Fill(event.GetQ2());
    h_xb->Fill(event.Getxb());
    h_y->Fill(event.Gety());
    h_nu->Fill(event.Getnu());
    h_W2->Fill(event.GetW2());
    h_vertexZ->Fill(event.GetVz());
    h_xQ2->Fill(event.Getxb(), event.GetQ2());
    h_thetaelectron->Fill(event.GetThetaElectron());
    h_rapport->Fill(event.GetAcosyada());

    //Particle electron = event.GetElectron();
    h_px_el->Fill(event.electron.GetMomentum().X());
    h_py_el->Fill(event.electron.GetMomentum().Y());
    h_pz_el->Fill(event.electron.GetMomentum().Z());
    h_ptot_el->Fill(event.electron.GetMomentum().P());
    h_theta_el->Fill(event.electron.GetMomentum().Theta()*180/Constants::PI);
    h_phi_el->Fill(event.electron.GetMomentum().Phi()*180/Constants::PI +180);
    h_polcoord_el->Fill(event.electron.GetMomentum().Theta()*180/Constants::PI, event.electron.GetMomentum().Phi()*180/Constants::PI +180);
    h_E_el->Fill(event.electron.GetMomentum().E());
    h_E_el_theta->Fill(event.electron.GetMomentum().Theta()*180/Constants::PI, event.electron.GetMomentum().E());
    h_E_el_phi->Fill(event.electron.GetMomentum().Phi()*180/Constants::PI +180, event.electron.GetMomentum().E());

    
    
    for (const Particle& hadron : event.GetHadrons()) {
        if (hadron.GetPID() == Constants::PION_PLUS_PID  ){   //adding this condition for pion+ and erasing the condit at evtprocessr
        h_z->Fill(hadron.Getz());
        h_pt2->Fill(hadron.Getpt2());
        h_phih->Fill(hadron.Getphih());
        h_px_pi->Fill(hadron.GetMomentum().X());
        h_py_pi->Fill(hadron.GetMomentum().Y());
        h_pz_pi->Fill(hadron.GetMomentum().Z());
        h_ptot_pi->Fill(hadron.GetMomentum().P());
        h_theta_pi->Fill(hadron.GetMomentum().Theta()*180/Constants::PI);
        h_phi_pi->Fill(hadron.GetMomentum().Phi()*180/Constants::PI +180);
        h_polcoord_pi->Fill(hadron.GetMomentum().Theta()*180/Constants::PI, hadron.GetMomentum().Phi()*180/Constants::PI +180);
        h_E_pi->Fill(hadron.GetMomentum().E());
        h_E_pi_theta->Fill(hadron.GetMomentum().Theta()*180/Constants::PI, hadron.GetMomentum().E());
        h_E_pi_phi->Fill(hadron.GetMomentum().Phi()*180/Constants::PI +180, hadron.GetMomentum().E());
        h_vertexZ_pi->Fill(hadron.GetParticleVertexZ());
        //std::cout<<"hadron vertex z = "<<hadron.GetParticleVertexZ()<<std::endl;
        //std::cout<<"electron vertex z = "<<event.GetVz()<<std::endl;

        h_DeltaVz->Fill(event.GetVz()-hadron.GetParticleVertexZ());
        if (event.GetHel() > 0){
            h_phih_plus->Fill(hadron.Getphih());
        }
        else if (event.GetHel() < 0){
            h_phih_minus->Fill(hadron.Getphih());
        }
        
        
        //
        //std::cout << " kinhad " << hadron.Getz()<<" , " << hadron.Getpt2()<<" , " << hadron.Getphih()<<std::endl;  
        //
        }
    }
}


void Monitoring::Fill_sR_Histograms(const Event& event  ) {
    //int monhel is an option if we want to have the helicity in histograms.
    int targetType = event.GetTargetType();
    int helicity = event.GetHel();  //only use this to monitor sinus

    if ( cut1.PassCutsElectrons(event)==true) {
        h_nu->Fill(event.Getnu(),helicity);

        for (const Particle& hadron : event.GetHadrons()) {
            if (cut1.PassCutsHadrons(hadron)==true){
                double phiD = hadron.Getphih();
                h_Q2->Fill(event.GetQ2(),sin(phiD)*helicity);
                h_xb->Fill(event.Getxb(),sin(phiD)*helicity);
                h_z->Fill(hadron.Getz(), sin(phiD)*helicity);
                //consider filling 3D histos for visuals
                //h_wD_Sratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), sin(phiD)*helicity) ;    //3 arguments and the WEIGHT
                //h_wD_sqSratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), sin(phiD)*sin(phiD)*helicity);    //3 arguments and the WEIGHT (pt2 squared) 4 variance
                //h_D_Sratio3D->Fill(event.Getxb(), event.GetQ2(), hadron.Getz()*helicity);    //3 arguments only counts not weight (cphi)

            }
        }
    }
}


    //else if (target == "Sn" && cuta.PassCutsElectrons(event)==true) {
    //    R_nu_el->Fill(event.Getnu());
    //    for (const Particle& hadron : event.GetHadrons()) {
    //        if (cuta.PassCutsHadrons(hadron)==true){
    //            R_nu_had->Fill(event.Getnu());
    //            R_z->Fill(hadron.Getz());
    //            R_pt2->Fill(hadron.Getpt2());
//
    //        }
    //    }
    //}

/*       DELETE THIS EFTER previous function is done 
void sratio::FillHistograms(const Event& event){
        //Retrieving helicity to consider "loss" or "cost" when value is -1. gain when value = +1
        //for 'loss' and 'gain' we'll add additionnal weight in count histograms. And adding a helicity factor to the histos already weighted. 
            counter_elLD2 ++;
            //set a counter that increases when electroncuts = passed; in order for it to be called when R is  computed in had variables (?) TBD
                }
            }
        }
        else if (targetType == 1 && cuta.PassCutsElectrons(event)==true) {
            counter_elSn++; //counter for only electrons for z and pt
            //here change the else if to just else in order to have a generic target 
            hSratio_nuA->Fill(event.Getnu());
            for (const Particle& hadron : event.GetHadrons()) {
                if (cuta.PassCutsHadrons(hadron)==true){
                    double phiA = hadron.Getphih();
                    h_wA_Sratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), sin(phiA)*helicity);    //3 arguments and the WEIGHT
                    h_wA_sqSratio->Fill(event.Getxb(), event.GetQ2(), hadron.Getz(), sin(phiA)*sin(phiA)*helicity);    //3 arguments and the WEIGHT (pt2 squared)
                    h_A_Sratio3D->Fill(event.Getxb(), event.GetQ2(), hadron.Getz()*helicity);    //3 arguments only counts not weight

                }
            }
        }
    }
*/

void Monitoring::WriteHistogramsToFile(const std::string filename) {
    //this function recreates a new rootfile everytime is called 
    //useful to have different rootfiles if different cuts were implemented
    //useful since monitoring class is called with a cutset as an argument
    //argument: title of the wanted output file 
    // save histograms to the specified ROOT file
    TFile file = TFile(filename.c_str(), "RECREATE");
    h_Q2->Write();
    h_xb->Write();
    h_y->Write();
    h_nu->Write();
    h_W2->Write();
    h_z->Write();
    h_pt2->Write();
    h_phih->Write();
    h_vertexZ->Write();
    //h_nuposel->Write();
    //h_Q2posel->Write();
    //h_Q2_had->Write();
    //h_nu_had->Write();
    //h_z_had->Write();
    //h_pt2_had->Write();
    h_pid->Write();
    TTree tree("treecounter","treecounter");
    tree.Branch("counterel_R", &counterel_R, "counterel_R/I");
    tree.Fill();
    file.Write();
    // Write other histograms here
    file.Close();
}

void Monitoring::DrawR_Histograms(const std::string filename) {
    TCanvas MonR("MonitoringR canvas", "Monitoring R Histograms");
    MonR.Divide(2, 2);
    MonR.cd(1);
    R_nu_el->SetMinimum(0);
    R_nu_el->Draw("hist");
    R_nu_el->SetTitle("nu(el) Distribution for R");
    MonR.cd(2);
    R_nu_had->SetMinimum(0);
    R_nu_had->Draw("hist");
    R_nu_had->SetTitle("nu(had) Distribution for R");
    MonR.cd(3);
    R_z->SetMinimum(0);
    R_z->Draw("hist");
    R_z->SetTitle("z Distribution for R");
    MonR.cd(4);
    R_pt2->SetMinimum(0);
    R_pt2->Draw("hist");
    R_pt2->SetTitle("pt2 Distribution for R");
    
    MonR.Print((filename + ".pdf").c_str());
    std::cout<<R_pt2->GetXaxis()->GetBinCenter(1)<<std::endl;
    std::cout<<R_pt2->GetXaxis()->GetBinCenter(2)<<std::endl;
    std::cout<<R_pt2->GetXaxis()->GetBinCenter(3)<<std::endl;
    std::cout<<R_pt2->GetXaxis()->GetBinCenter(4)<<std::endl;
    std::cout<<R_pt2->GetXaxis()->GetBinCenter(5)<<std::endl;
    

}
void Monitoring::SaveKeyHistograms() {
    // Disable the statistics box for these histograms only
    gStyle->SetOptStat(0);

    // Save the xB vs Q^2 histogram
    TCanvas canvas_xQ2("canvas_xQ2", "", 800, 600);
    h_xQ2->SetTitle("");  // Remove title
    h_xQ2->GetXaxis()->SetTitle("x_{B}");
    h_xQ2->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
    h_xQ2->GetXaxis()->SetTitleSize(0.05);  // Set X-axis title size
    h_xQ2->GetYaxis()->SetTitleSize(0.05);  // Set Y-axis title size
    h_xQ2->Draw("COLZ");
    canvas_xQ2.SaveAs("h_xQ2.png");  // Save PNG output

    // Save the z vs pT^2 histogram
    TCanvas canvas_pt2z("canvas_pt2z", "", 800, 600);
    h_pt2z->SetTitle("");  // Remove title
    h_pt2z->GetXaxis()->SetTitle("p_{t}^{2} (GeV^{2})");
    h_pt2z->GetYaxis()->SetTitle("z");
    h_pt2z->GetXaxis()->SetTitleSize(0.05);  // Set X-axis title size
    h_pt2z->GetYaxis()->SetTitleSize(0.05);  // Set Y-axis title size
    h_pt2z->Draw("COLZ");
    canvas_pt2z.SaveAs("h_pt2z.png");  // Save PNG output
}


void Monitoring::DrawHistograms(const std::string filename) {
    
    // Initialize the PDF file
    TCanvas canvas("canvas", "Monitoring Histograms", 800, 600);
    canvas.Divide(3, 3);

    // First page: All histograms with TLines
    canvas.cd(1);
    h_Q2->SetTitle("Q^{2} Distribution");
    h_Q2->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    h_Q2->GetXaxis()->SetRangeUser(0, 5);
    h_Q2->SetMinimum(0);
    h_Q2->Draw("hist");
    TLine *line_Q2 = new TLine(Constants::RcutminQ, h_Q2->GetMinimum(), Constants::RcutminQ, h_Q2->GetMaximum());
    line_Q2->SetLineStyle(2);
    line_Q2->Draw();

    canvas.cd(2);
    h_xb->SetTitle("xb Distribution");
    h_xb->GetXaxis()->SetTitle("x_{B}");
    h_xb->GetXaxis()->SetRangeUser(0, 0.6);
    h_xb->SetMinimum(0);
    h_xb->Draw("hist");

    canvas.cd(3);
    h_y->SetTitle("y Distribution");
    h_y->GetXaxis()->SetTitle("y");
    h_y->SetMinimum(0);
    h_y->Draw("hist");
    TLine *line_y = new TLine(Constants::RcutminY, h_y->GetMinimum(), Constants::RcutminY, h_y->GetMaximum());
    line_y->SetLineStyle(2);
    line_y->Draw();

    canvas.cd(4);
    h_nu->SetTitle("nu Distribution");
    h_nu->GetXaxis()->SetTitle("nu (GeV)");
    h_nu->GetXaxis()->SetRangeUser(0, 10);
    h_nu->SetMinimum(0);
    h_nu->Draw("hist");

    canvas.cd(5);
    h_W2->SetTitle("W^{2} Distribution");
    h_W2->GetXaxis()->SetTitle("W^{2} (GeV^{2})");
    h_W2->SetMinimum(0);
    h_W2->Draw("hist");
    TLine *line_W2 = new TLine(Constants::RcutminW, h_W2->GetMinimum(), Constants::RcutminW, h_W2->GetMaximum());
    line_W2->SetLineStyle(2);
    line_W2->Draw();

    canvas.cd(6);
    h_z->SetTitle("z Distribution");
    h_z->GetXaxis()->SetTitle("z");
    h_z->SetMinimum(0);
    h_z->Draw("hist");
    TLine *line_z = new TLine(Constants::RcutminZ, h_z->GetMinimum(), Constants::RcutminZ, h_z->GetMaximum());
    TLine *line_zmax = new TLine(Constants::RcutmaxZ, h_z->GetMinimum(), Constants::RcutmaxZ, h_z->GetMaximum());
    line_z->SetLineStyle(2);
    line_zmax->SetLineStyle(2);
    line_z->Draw();
    line_zmax->Draw();

    canvas.cd(7);
    h_pt2->SetTitle("pt2 Distribution");
    h_pt2->GetXaxis()->SetTitle("p_{t}^{2} (GeV^{2})");
    h_pt2->SetMinimum(0);
    h_pt2->Draw("hist");
    TLine *line_pt2 = new TLine(Constants::RcutminPt2, h_pt2->GetMinimum(), Constants::RcutminPt2, h_pt2->GetMaximum());
    TLine *line_pt2max = new TLine(Constants::RcutmaxPt2, h_pt2->GetMinimum(), Constants::RcutmaxPt2, h_pt2->GetMaximum());
    line_pt2->SetLineStyle(2);
    line_pt2max->SetLineStyle(2);
    line_pt2->Draw();
    line_pt2max->Draw();

    canvas.cd(8);
    h_phih->SetTitle("phih Distribution");
    h_phih->GetXaxis()->SetTitle("phi_{h}");
    h_phih->SetMinimum(0);
    h_phih->Draw("hist");

    canvas.cd(9);
    h_vertexZ->SetTitle("Vertex Z Distribution");
    h_vertexZ->GetXaxis()->SetTitle("Z_{vertex}");
    h_vertexZ->SetMinimum(0);
    h_vertexZ->Draw("hist");

    // Save the first page
    canvas.Print((filename + ".pdf(").c_str()); // Open PDF and save the first page

    // Second page: Selected histograms and 2D distributions
    canvas.Clear();
    canvas.Divide(3, 3);

    canvas.cd(1);
    h_Q2->SetTitle("Q^{2} Distribution");
    h_Q2->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    h_Q2->GetXaxis()->SetRangeUser(0, 5);
    h_Q2->SetMinimum(0);
    h_Q2->Draw("hist");
    TLine *line_Q2bis = new TLine(Constants::RcutminQ, h_Q2->GetMinimum(), Constants::RcutminQ, h_Q2->GetMaximum());
    line_Q2bis->SetLineStyle(2);
    line_Q2bis->Draw();

    canvas.cd(2);
    h_xb->SetTitle("x_{B} Distribution");
    h_xb->GetXaxis()->SetTitle("x_{B}");
    h_xb->GetXaxis()->SetRangeUser(0, 0.6);
    h_xb->SetMinimum(0);
    h_xb->Draw("hist");

    canvas.cd(3);
    h_z->SetTitle("z Distribution");
    h_z->GetXaxis()->SetTitle("z");
    h_z->SetMinimum(0);
    h_z->Draw("hist");
    TLine *line_zbis = new TLine(Constants::RcutminZ, h_z->GetMinimum(), Constants::RcutminZ, h_z->GetMaximum());
    TLine *line_zmaxbis = new TLine(Constants::RcutmaxZ, h_z->GetMinimum(), Constants::RcutmaxZ, h_z->GetMaximum());
    line_zbis->SetLineStyle(2);
    line_zmaxbis->SetLineStyle(2);
    line_zbis->Draw();
    line_zmaxbis->Draw();

    canvas.cd(4);
    h_pt2->SetTitle("p_{t}^{2} Distribution");
    h_pt2->GetXaxis()->SetTitle("p_{t}^{2} (GeV^{2})");
    h_pt2->SetMinimum(0);
    h_pt2->Draw("hist");
    TLine *line_pt2bis = new TLine(Constants::RcutminPt2, h_pt2->GetMinimum(), Constants::RcutminPt2, h_pt2->GetMaximum());
    TLine *line_pt2maxbis = new TLine(Constants::RcutmaxPt2, h_pt2->GetMinimum(), Constants::RcutmaxPt2, h_pt2->GetMaximum());
    line_pt2bis->SetLineStyle(2);
    line_pt2maxbis->SetLineStyle(2);
    line_pt2bis->Draw();
    line_pt2maxbis->Draw();

    canvas.cd(5);
    h_xQ2->SetTitle("x_{B} vs Q^{2}");
    h_xQ2->GetXaxis()->SetTitle("x_{B}");
    h_xQ2->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
    h_xQ2->Draw("COLZ");

    canvas.cd(6);
    h_pt2z->SetTitle("z vs p_{t}^{2}");
    h_pt2z->GetXaxis()->SetTitle("p_{t}^{2} (GeV^{2})");
    h_pt2z->GetYaxis()->SetTitle("z");
    h_pt2z->Draw("COLZ");

//    canvas.cd(7);
//    h_thetaelectron->SetTitle("Theta Electron");
//    h_thetaelectron->Draw("hist");
//
//    canvas.cd(8);
//    h_thetaelectronMC->SetTitle("Theta Electron MC");
//    h_thetaelectronMC->Draw("hist");
//
//    canvas.cd(9);
//    h_chi2_el->SetTitle("Chi2 Electron");
//    h_chi2_el->Draw("hist");

    // Save the second page and close the PDF
    canvas.Print((filename + ".pdf)").c_str()); // Close the PDF
}





void Monitoring::DrawQmonitoring(Monitoring& monTrue, const std::string filename){
    TCanvas MonQ("Monitoring Q canvas", "Monitoring Q Histograms");
    MonQ.Divide(2, 2);
    MonQ.cd(1);
    //h_rapportMC->Draw("hist");
    h_rapportMC->SetLineColor(kGreen); 
    h_rapport->Draw("hist ");
    h_rapport->SetLineColor(kRed);
    monTrue.h_rapport->Draw("hist same");
    h_rapport->SetTitle("rapport beam and el Distribution");

    MonQ.cd(2);
    //h_thetaelectronMC->Draw("hist");
    h_thetaelectronMC->SetLineColor(kGreen);
    h_thetaelectron->Draw("hist ");
    h_thetaelectron->SetLineColor(kRed);
    monTrue.h_thetaelectron->Draw("hist same");
    h_thetaelectron->SetTitle("acos (rapport)");

    MonQ.cd(3);
    h_Q2->Draw("hist");
    h_Q2->SetLineColor(kRed); // Set different color for monRec histogram
    h_Q2->SetTitle("Q2 Distribution");
    monTrue.h_Q2->Draw("hist same");
    monTrue.h_Q2->SetLineColor(kBlue);
    MonQ.Print((filename + ".pdf").c_str());
}



void Monitoring::DrawHistTrueandSIM(Monitoring& monTrue  , const std::string filename) {
    TCanvas MonC("Monitoring canvas", "Monitoring Histograms");
    MonC.Divide(3, 3);


     MonC.cd(1);
    h_Q2->SetLineColor(kRed); // Set different color for monRec histogram
    h_Q2MC->Draw("hist ");

    monTrue.h_Q2->Draw("hist same");

    h_Q2->Draw("hist same");
    h_Q2MC->SetLineColor(kGreen);
    
    monTrue.h_Q2->SetTitle("Q2 Distribution");
    monTrue.h_Q2->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    monTrue.h_Q2->GetXaxis()->SetRangeUser(1, 5);
    TLine *line_Q2_true = new TLine(Constants::RcutminQ, monTrue.h_Q2->GetMinimum(), Constants::RcutminQ, monTrue.h_Q2->GetMaximum());
    line_Q2_true->SetLineStyle(2); // Dashed line style
    line_Q2_true->Draw();

    //monRec.h_Q2->SetTitle("Q2 Distribution - Rec");
    //monRec.h_Q2->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    //monRec.h_Q2->GetXaxis()->SetRangeUser(1, 5);
    TLine *line_Q2_rec = new TLine(Constants::RcutminQ, h_Q2->GetMinimum(), Constants::RcutminQ, h_Q2->GetMaximum());
    line_Q2_rec->SetLineStyle(2); // Dashed line style
    line_Q2_rec->SetLineColor(kRed); // Set same color as histogram for consistency
    
    line_Q2_rec->Draw();

    // Create legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Position of the legend
    legend->AddEntry(monTrue.h_Q2, "True Data", "l");
    legend->AddEntry(h_Q2, "Sim Rec", "l");
    legend->AddEntry (h_Q2MC, "Sim MC", "l");
    legend->SetBorderSize(0); // Remove legend border
    legend->Draw();

    //now we should propagate the DrawHistograms function adapting it to this function 
    //we should add the same for the other histograms

    MonC.cd(2);
    h_xbMC->Draw("hist ");
    h_xb->SetLineColor(kRed); // Set different color for monRec histogram
    h_xbMC->SetLineColor(kGreen);
    h_xb->Draw("hist same");
    monTrue.h_xb->Draw("hist same");
    monTrue.h_xb->SetTitle("xb Distribution");
    monTrue.h_xb->GetXaxis()->SetTitle("x_{B}");
    monTrue.h_xb->GetXaxis()->SetRangeUser(0, 0.6);
    
    TLegend *legendx = new TLegend(0.7, 0.7, 0.9, 0.9); // Position of the legend
    legendx->AddEntry(monTrue.h_xb, "True Data", "l");
    legendx->AddEntry(h_xb, "Sim Rec", "l");
    legendx->AddEntry (h_xbMC, "Sim MC", "l");

    MonC.cd(3);
    h_yMC->Draw("hist ");
    h_yMC->SetLineColor(kGreen);
    h_y->SetLineColor(kRed); // Set different color for monRec histogram
    h_y->Draw("hist same ");
    monTrue.h_y->Draw("hist same");
    monTrue.h_y->SetTitle("y Distribution");
    monTrue.h_y->GetXaxis()->SetTitle("y");
    TLine *line_y_true = new TLine(Constants::RcutminY, monTrue.h_y->GetMinimum(), Constants::RcutminY, monTrue.h_y->GetMaximum());
    line_y_true->SetLineStyle(2); // Dashed line style
    line_y_true->Draw();
    TLine *line_y_rec = new TLine(Constants::RcutminY, h_y->GetMinimum(), Constants::RcutminY, h_y->GetMaximum());
    line_y_rec->SetLineStyle(2); // Dashed line style
    line_y_rec->SetLineColor(kRed); // Set same color as histogram for consistency
    line_y_rec->Draw();
    TLegend *legendy = new TLegend(0.7, 0.7, 0.9, 0.9); // Position of the legend
    legendy->AddEntry(monTrue.h_y, "True Data", "l");
    legendy->AddEntry(h_y, "Sim Rec", "l");
    legendy->AddEntry (h_yMC, "Sim MC", "l");

    MonC.cd(4);
    h_nuMC->Draw("hist ");
    h_nuMC->SetLineColor(kGreen);
    h_nu->SetLineColor(kRed); // Set different color for monRec histogram
    h_nu->Draw("hist same ");
    monTrue.h_nu->Draw("hist same");
    monTrue.h_nu->SetTitle("nu Distribution");
    monTrue.h_nu->GetXaxis()->SetTitle("nu (GeV)"); 
    monTrue.h_nu->GetXaxis()->SetRangeUser(3, 10);
    TLine *line_nu_true = new TLine(Constants::Rcutminnu, monTrue.h_nu->GetMinimum(), Constants::Rcutminnu, monTrue.h_nu->GetMaximum());
    line_nu_true->SetLineStyle(2); // Dashed line style
    line_nu_true->Draw();
    TLine *line_nu_rec = new TLine(Constants::Rcutminnu, h_nu->GetMinimum(), Constants::Rcutminnu, h_nu->GetMaximum());
    line_nu_rec->SetLineStyle(2); // Dashed line style
    line_nu_rec->SetLineColor(kRed); // Set same color as histogram for consistency
    line_nu_rec->Draw();
    TLegend *legendnu = new TLegend(0.7, 0.7, 0.9, 0.9); // Position of the legend
    legendnu->AddEntry(monTrue.h_nu, "True Data", "l");
    legendnu->AddEntry(h_nu, "Sim Rec", "l");
    legendnu->AddEntry (h_nuMC, "Sim MC", "l");

    MonC.cd(5);
    h_W2MC->Draw("hist ");
    h_W2MC->SetLineColor(kGreen);
    h_W2->SetLineColor(kRed); // Set different color for monRec histogram
    h_W2->Draw("hist same");
    monTrue.h_W2->Draw("hist same");
    monTrue.h_W2->SetTitle("W^{2} Distribution");
    monTrue.h_W2->GetXaxis()->SetTitle("W^{2} (GeV^{2})");
    TLine *line_W2_true = new TLine(Constants::RcutminW, monTrue.h_W2->GetMinimum(), Constants::RcutminW, monTrue.h_W2->GetMaximum());
    line_W2_true->SetLineStyle(2); // Dashed line style
    TLine *line_W2_rec = new TLine(Constants::RcutminW, h_W2->GetMinimum(), Constants::RcutminW, h_W2->GetMaximum());
    line_W2_rec->SetLineStyle(2); // Dashed line style
    line_W2_rec->SetLineColor(kRed); // Set same color as histogram for consistency
    line_W2_rec->Draw();
    TLegend *legendW2 = new TLegend(0.7, 0.7, 0.9, 0.9); // Position of the legend
    legendW2->AddEntry(monTrue.h_W2, "True Data", "l");
    legendW2->AddEntry(h_W2, "Sim Rec", "l");
    legendW2->AddEntry (h_W2MC, "Sim MC", "l");

    MonC.cd(6);
    h_zMC->Draw("hist ");
    h_zMC->SetLineColor(kGreen);
    h_z->SetLineColor(kRed); // Set different color for monRec histogram
    h_z->Draw("hist same ");
    monTrue.h_z->Draw("hist same");
    monTrue.h_z->SetTitle("z Distribution");
    monTrue.h_z->GetXaxis()->SetTitle("z");
    TLine *line_z_true = new TLine(Constants::RcutminZ, monTrue.h_z->GetMinimum(), Constants::RcutminZ, monTrue.h_z->GetMaximum());
    line_z_true->SetLineStyle(2); // Dashed line style
    line_z_true->Draw();
    TLine *line_zmax_true = new TLine(Constants::RcutmaxZ, monTrue.h_z->GetMinimum(), Constants::RcutmaxZ, monTrue.h_z->GetMaximum());
    line_zmax_true->SetLineStyle(2); // Dashed line style
    line_zmax_true->Draw();
    TLine *line_z_rec = new TLine(Constants::RcutminZ, h_z->GetMinimum(), Constants::RcutminZ, h_z->GetMaximum());
    line_z_rec->SetLineStyle(2); // Dashed line style
    line_z_rec->SetLineColor(kRed); // Set same color as histogram for consistency
    line_z_rec->Draw();
    TLine *line_zmax_rec = new TLine(Constants::RcutmaxZ, h_z->GetMinimum(), Constants::RcutmaxZ, h_z->GetMaximum());
    line_zmax_rec->SetLineStyle(2); // Dashed line style
    line_zmax_rec->SetLineColor(kRed); // Set same color as histogram for consistency
    line_zmax_rec->Draw();
    TLegend *legendz = new TLegend(0.7, 0.7, 0.9, 0.9); // Position of the legend
    legendz->AddEntry(monTrue.h_z, "True Data", "l");
    legendz->AddEntry(h_z, "Sim Rec", "l");
    legendz->AddEntry (h_zMC, "Sim MC", "l");


    MonC.cd(7);
    h_pt2MC->Draw("hist ");
    h_pt2MC->SetLineColor(kGreen);
    h_pt2->SetLineColor(kRed); // Set different color for monRec histogram
    h_pt2->Draw("hist same");
    monTrue.h_pt2->Draw("hist same");
    monTrue.h_pt2->SetTitle("pt2 Distribution");
    monTrue.h_pt2->GetXaxis()->SetTitle("p_{t}^{2} (GeV^{2})");
    TLine *line_pt2_true = new TLine(Constants::RcutminPt2, monTrue.h_pt2->GetMinimum(), Constants::RcutminPt2, monTrue.h_pt2->GetMaximum());
    line_pt2_true->SetLineStyle(2); // Dashed line style
    line_pt2_true->Draw();
    TLine *line_pt2max_true = new TLine(Constants::RcutmaxPt2, monTrue.h_pt2->GetMinimum(), Constants::RcutmaxPt2, monTrue.h_pt2->GetMaximum());
    line_pt2max_true->SetLineStyle(2); // Dashed line style
    line_pt2max_true->Draw();
    TLine *line_pt2_rec = new TLine(Constants::RcutminPt2, h_pt2->GetMinimum(), Constants::RcutminPt2, h_pt2->GetMaximum());
    line_pt2_rec->SetLineStyle(2); // Dashed line style
    line_pt2_rec->SetLineColor(kRed); // Set same color as histogram for consistency
    line_pt2_rec->Draw();
    TLine *line_pt2max_rec = new TLine(Constants::RcutmaxPt2, h_pt2->GetMinimum(), Constants::RcutmaxPt2, h_pt2->GetMaximum());
    line_pt2max_rec->SetLineStyle(2); // Dashed line style
    line_pt2max_rec->SetLineColor(kRed); // Set same color as histogram for consistency
    line_pt2max_rec->Draw();
    TLegend *legendpt2 = new TLegend(0.7, 0.7, 0.9, 0.9); // Position of the legend
    legendpt2->AddEntry(monTrue.h_pt2, "True Data", "l");
    legendpt2->AddEntry(h_pt2, "Sim Rec", "l");
    legendpt2->AddEntry (h_pt2MC, "Sim MC", "l");

    MonC.cd(8);
    h_phihMC->Draw("hist ");
    h_phihMC->SetLineColor(kGreen);
    h_phih->SetLineColor(kRed); // Set different color for monRec histogram
    h_phih->Draw("hist same ");
    monTrue.h_phih->Draw("hist same ");
    monTrue.h_phih->SetTitle("phih Distribution");
    monTrue.h_phih->GetXaxis()->SetTitle("phi_{h} (degrees)");
    TLegend *legendphi = new TLegend(0.7, 0.7, 0.9, 0.9); // Position of the legend
    legendphi->AddEntry(monTrue.h_phih, "True Data", "l");
    legendphi->AddEntry(h_phih, "Sim Rec", "l");
    legendphi->AddEntry (h_phihMC, "Sim MC", "l");

MonC.cd(9);
    h_vertexZMC->Draw("hist ");
    h_vertexZMC->SetLineColor(kGreen);
    h_vertexZ->SetLineColor(kRed); // Set different color for monRec histogram
    h_vertexZ->Draw("hist same");
    monTrue.h_vertexZ->Draw("hist same");
    h_vertexZ->SetTitle("Vertex Z Distribution");
    TLegend *legendvertexZ = new TLegend(0.7, 0.7, 0.9, 0.9); // Position of the legend
    legendvertexZ->AddEntry(monTrue.h_vertexZ, "True Data", "l");
    legendvertexZ->AddEntry(h_vertexZ, "Sim Rec", "l");
    legendvertexZ->AddEntry (h_vertexZMC, "Sim MC", "l");


    MonC.Print((filename + ".pdf").c_str());

   
}




void Monitoring::DrawHistogramsPos(const std::string targetName, const std::string filename){
    TCanvas MonC("Monitoring canvas", "Monitoring Histograms");
    MonC.Divide(2, 2);
    MonC.cd(1);
    h_xQ2->Draw("COLZ");
    h_xQ2->SetTitle(("x vs Q2 Distribution (" + targetName + ")").c_str());
    h_xQ2->GetXaxis()->SetTitle("x_{B}");
    h_xQ2->GetYaxis()->SetTitle("Q^{2}");
    h_xQ2->GetXaxis()->SetRangeUser(0, 1);
    //h_xQ2->SetBins(nubin, 0, 1); // 
    MonC.cd(2);
    h_xQ2pos->Draw("COLZ");
    h_xQ2pos->SetTitle(("x vs Q2 Distribution (" + targetName + ")").c_str());
    h_xQ2pos->GetXaxis()->SetTitle("x_{B}");
    h_xQ2pos->GetYaxis()->SetTitle("Q^{2}");
    


//    MonC.cd(3);
//    MonC.cd(4);


    MonC.Print((filename + ".pdf").c_str());
}




void Monitoring::DrawCaloHistograms(const std::string filename) {
    TCanvas MonCal("Monitoring calocanvas", "Monitoring caloHistograms");
    MonCal.Divide(3, 3);
    MonCal.cd(1);
    h_lu->Draw("hist");
    h_lu->SetTitle("lu el_LD2");
    h_lu->GetXaxis()->SetTitle("u (cm)");
    MonCal.cd(2);
    h_lv->Draw("hist");
    h_lv->SetTitle("lv el_LD2");
    h_lv->GetXaxis()->SetTitle("v (cm)");
    TLine *line_lv = new TLine(Constants::cutlv_min, h_lv->GetMinimum(), Constants::cutlv_min, h_lv->GetMaximum());
    line_lv->SetLineStyle(2); // Dashed line style
    line_lv->Draw();
    MonCal.cd(3);
    h_lw->Draw("hist");
    h_lw->SetTitle("lw el_LD2");
    h_lw->GetXaxis()->SetTitle("w (cm)");
    TLine *line_lw = new TLine(Constants::cutlw_min, h_lw->GetMinimum(), Constants::cutlw_min, h_lw->GetMaximum());
    line_lw->SetLineStyle(2); // Dashed line style
    line_lw->Draw();
    MonCal.cd(4);
    h_epcal->Draw("hist");
    h_epcal->SetTitle("epcal el_LD2");
    h_epcal->GetXaxis()->SetTitle("E (GeV)");
    TLine *line_epcal = new TLine(Constants::cutEpcal_min, h_epcal->GetMinimum(), Constants::cutEpcal_min, h_epcal->GetMaximum());
    line_epcal->SetLineStyle(2); // Dashed line style
    line_epcal->Draw();
    MonCal.cd(5);
    h_eecalin->Draw("hist");
    h_eecalin->SetTitle("Ecal_inner el_LD2");
    h_eecalin->GetXaxis()->SetTitle("E (GeV)");
    MonCal.cd(6);
    h_epcalout->Draw("hist");
    h_epcalout->SetTitle("Ecal_outer el_LD2");
    h_epcalout->GetXaxis()->SetTitle("E (GeV)");
    MonCal.cd(7);
    h_calEall->Draw("COLZ");
    h_calEall->SetTitle("Ecal vs Epcal el_LD2");
    h_calEall->GetXaxis()->SetTitle("Epcal (GeV)");
    h_calEall->GetYaxis()->SetTitle("Ecal In+Out (GeV)");
    TLine *line_calEall = new TLine(Constants::cutEpcal_min, h_calEall->GetYaxis()->GetXmin(), Constants::cutEpcal_min, h_calEall->GetYaxis()->GetXmax());
    line_calEall->SetLineStyle(2); 
    line_calEall->SetLineColor(kRed); 
    line_calEall->SetLineWidth(2); 
    line_calEall->Draw();
    //TLine *line_Ecal = new TLine(Constants::cutEpcal_min, h_epcal->GetMinimum(), Constants::cutEpcal_min, h_epcal->GetMaximum());
    MonCal.cd(8);
    h_calXY->Draw("COLZ");
    h_calXY->SetTitle("PCalXY");
    h_calXY->GetXaxis()->SetTitle("x (cm)");
    h_calXY->GetYaxis()->SetTitle("y (cm)");
    MonCal.cd(9);
    h_calSector->Draw("hist");
    h_calSector->SetTitle("calSector el_LD2");
    h_calSector->GetXaxis()->SetTitle("Sectors (1 - 6)");

    MonCal.Print((filename + ".pdf").c_str());
}


void Monitoring::DrawCherenkovHistograms(const std::string filename) {
    TCanvas MonCher("Monitoring cherencanvas", "Monitoring cherenkovHistograms");
    MonCher.Divide(2, 1);
    MonCher.cd(1);
    h_Nphe15->Draw("hist");
    h_Nphe15->SetTitle("Nphe HTCC el_LD2");
    h_Nphe15->GetXaxis()->SetTitle("Nphe 15");
    TLine *line_Nphe15 = new TLine(Constants::Nphe15min, h_Nphe15->GetMinimum(), Constants::Nphe15min, h_Nphe15->GetMaximum());
    line_Nphe15->SetLineStyle(2); // Dashed line style
    line_Nphe15->Draw();
    MonCher.cd(2);
    h_Nphe16->Draw("hist");
    h_Nphe16->GetXaxis()->SetTitle("Nphe 16");
    h_Nphe16->SetTitle("Nphe LTCC el_LD2");
    //TLine *line_Nphe16 = new TLine(Constants::Nphe16min, h_Nphe16->GetMinimum(), Constants::Nphe16min, h_Nphe16->GetMaximum());
    //line_Nphe16->SetLineStyle(2); // Dashed line style
    //line_Nphe16->Draw();
    


    MonCher.Print((filename + ".pdf").c_str());
}


void Monitoring::DrawMomentumHistograms(const std::string filename){
    TCanvas MonMom("Monitoring momcanvas", "Monitoring momHistograms");
    MonMom.Divide(4, 3);
    MonMom.cd(1);
    h_px_el->Draw("hist");
    h_px_el->SetTitle("px el_LD2");
    h_px_el->GetXaxis()->SetTitle("px (GeV)");
    MonMom.cd(2);
    h_py_el->Draw("hist");
    h_py_el->SetTitle("py el_LD2");
    h_py_el->GetXaxis()->SetTitle("py (GeV)");
    MonMom.cd(3);
    h_pz_el->Draw("hist");
    h_pz_el->SetTitle("pz el_LD2");
    h_pz_el->GetXaxis()->SetTitle("pz (GeV)");
    MonMom.cd(4);
    h_px_pi->Draw("hist");
    h_px_pi->SetTitle("px pi_LD2");
    h_px_pi->GetXaxis()->SetTitle("px (GeV)");
    MonMom.cd(5);
    h_py_pi->Draw("hist");
    h_py_pi->SetTitle("py pi_LD2");
    h_py_pi->GetXaxis()->SetTitle("py (GeV)");
    MonMom.cd(6);
    h_pz_pi->Draw("hist");
    h_pz_pi->SetTitle("pz pi_LD2");
    h_pz_pi->GetXaxis()->SetTitle("pz (GeV)");
    MonMom.cd(7);
    h_theta_el->Draw("hist");
    h_theta_el->SetTitle("theta el_LD2");
    h_theta_el->GetXaxis()->SetTitle("theta (degrees)");
    MonMom.cd(8);
    h_phi_el->Draw("hist");
    h_phi_el->SetTitle("phi el_LD2");
    h_phi_el->GetXaxis()->SetTitle("phi (degrees)");
    MonMom.cd(9);
    h_theta_pi->Draw("hist");
    h_theta_pi->SetTitle("theta pi_LD2");
    h_theta_pi->GetXaxis()->SetTitle("theta (degrees)");
    MonMom.cd(10);
    h_phi_pi->Draw("hist");
    h_phi_pi->SetTitle("phi pi_LD2");
    h_phi_pi->GetXaxis()->SetTitle("phi (degrees)");
    MonMom.cd(11);
    h_polcoord_el->Draw("COLZ");
    h_polcoord_el->SetTitle("polar el_LD2");
    h_polcoord_el->GetXaxis()->SetTitle("theta (degrees)");
    h_polcoord_el->GetYaxis()->SetTitle("phi (degrees)");
    MonMom.cd(12);
    h_polcoord_pi->Draw("COLZ");
    h_polcoord_pi->SetTitle("polar pi_LD2");
    h_polcoord_pi->GetXaxis()->SetTitle("theta (degrees)");
    h_polcoord_pi->GetYaxis()->SetTitle("phi (degrees)");

    

    MonMom.Print((filename + ".pdf").c_str());

}

void Monitoring::DrawMomentumElectronHistograms(const std::string filename){
    TCanvas MonMomEl("Monitoring momelcanvas", "Monitoring momelHistograms");
    MonMomEl.Divide(4, 2);
    MonMomEl.cd(1);
    h_px_el->Draw("hist");
    h_px_el->SetTitle("px el_LD2");
    h_px_el->GetXaxis()->SetTitle("px (GeV)");
    MonMomEl.cd(2);
    h_py_el->Draw("hist");
    h_py_el->SetTitle("py el_LD2");
    h_py_el->GetXaxis()->SetTitle("py (GeV)");
    MonMomEl.cd(3);
    h_pz_el->Draw("hist");
    h_pz_el->SetTitle("pz el_LD2");
    h_pz_el->GetXaxis()->SetTitle("pz (GeV)");
    MonMomEl.cd(4);
    //Ptotal histo 
    h_ptot_el->Draw("hist");
    h_ptot_el->SetTitle("ptot el_LD2");
    h_ptot_el->GetXaxis()->SetTitle("ptot (GeV)");
    MonMomEl.cd(5);
    h_theta_el->Draw("hist");
    h_theta_el->SetTitle("theta el_LD2");
    h_theta_el->GetXaxis()->SetTitle("theta (degrees)");
    MonMomEl.cd(6);
    h_phi_el->Draw("hist");
    h_phi_el->SetTitle("phi el_LD2");
    h_phi_el->GetXaxis()->SetTitle("phi (degrees)");
    MonMomEl.cd(7);
    h_polcoord_el->Draw("COLZ");
    h_polcoord_el->SetTitle("polar el_LD2");
    h_polcoord_el->GetXaxis()->SetTitle("theta (degrees)");
    h_polcoord_el->GetYaxis()->SetTitle("phi (degrees)");
    MonMomEl.cd(8);
    //Vertex histo plots
    h_vertexZ->Draw("hist");
    h_vertexZ->SetTitle("Vertex Z Distribution");
    h_vertexZ->GetXaxis()->SetTitle("V_{z} (cm)");
    h_vertexZ->GetXaxis()->SetRangeUser(-20,10);
    //TLine *line_vzmin = new TLine(Constants::cutvz_min, h_vertexZ->GetMinimum(), Constants::cutvz_min, h_vertexZ->GetMaximum());
    //line_vzmin->SetLineStyle(2); // Dashed line style
    //line_vzmin->Draw();
    //TLine *line_vzmax = new TLine(Constants::cutvz_max, h_vertexZ->GetMinimum(), Constants::cutvz_max, h_vertexZ->GetMaximum());
    //line_vzmax->SetLineStyle(2); // Dashed line style
    //line_vzmax->Draw();

    MonMomEl.Print((filename + ".pdf").c_str());

}

void Monitoring::DrawMomentumHadronHistograms(const std::string filename){
    TCanvas MonMomHad("Monitoring momhadcanvas", "Monitoring momhadHistograms");
    MonMomHad.Divide(4, 2);
    MonMomHad.cd(1);
    h_px_pi->Draw("hist");
    h_px_pi->SetTitle("px pi_LD2");
    h_px_pi->GetXaxis()->SetTitle("px (GeV)");
    MonMomHad.cd(2);
    h_py_pi->Draw("hist");
    h_py_pi->SetTitle("py pi_LD2");
    h_py_pi->GetXaxis()->SetTitle("py (GeV)");
    MonMomHad.cd(3);
    h_pz_pi->Draw("hist");
    h_pz_pi->SetTitle("pz pi_LD2");
    h_pz_pi->GetXaxis()->SetTitle("pz (GeV)");
    MonMomHad.cd(4);
    //Ptotal histo 
    h_ptot_pi->Draw("hist");
    h_ptot_pi->SetTitle("ptot pi_LD2");
    h_ptot_pi->GetXaxis()->SetTitle("ptot (GeV)"); 
    MonMomHad.cd(5);
    h_theta_pi->Draw("hist");
    h_theta_pi->SetTitle("theta pi_LD2");
    h_theta_pi->GetXaxis()->SetTitle("theta (degrees)");
    MonMomHad.cd(6);
    h_phi_pi->Draw("hist");
    h_phi_pi->SetTitle("phi pi_LD2");
    h_phi_pi->GetXaxis()->SetTitle("phi (degrees)");
    MonMomHad.cd(7);
    h_polcoord_pi->Draw("COLZ");
    h_polcoord_pi->SetTitle("polar pi_LD2");
    h_polcoord_pi->GetXaxis()->SetTitle("theta (degrees)");
    h_polcoord_pi->GetYaxis()->SetTitle("phi (degrees)");
    MonMomHad.cd(8);
    //Vertex histo plots
    h_vertexZ->Draw("hist");
    h_vertexZ->SetTitle("Vertex Z Distribution");
    h_vertexZ->GetXaxis()->SetTitle("V_{z} (cm)");
    h_vertexZ->GetXaxis()->SetRangeUser(-20,10);

    MonMomHad.Print((filename + ".pdf").c_str());

}


void Monitoring::DrawEnergyHistograms(const std::string filename){
    TCanvas MonEnrg("Monitoring Enrgcanvas", "Monitoring EnrgHistograms");
    MonEnrg.Divide(3, 2);
    MonEnrg.cd(1);
    h_E_el->Draw("hist");
    h_E_el->SetTitle("E el_LD2");
    h_E_el->GetXaxis()->SetTitle("E electron (GeV)");
    MonEnrg.cd(2);
    h_E_el_theta->Draw("COLZ");
    h_E_el_theta->SetTitle("E vs theta el_LD2");
    h_E_el_theta->GetXaxis()->SetTitle("theta electron (degrees)");
    MonEnrg.cd(3);
    h_E_el_phi->Draw("COLZ");
    h_E_el_phi->SetTitle("E vs phi el_LD2");
    h_E_el_phi->GetXaxis()->SetTitle("phi electron (degrees)");
    MonEnrg.cd(4);
    h_E_pi->Draw("hist");
    h_E_pi->SetTitle("E pi_LD2");
    h_E_pi->GetXaxis()->SetTitle("E pion (GeV)");
    MonEnrg.cd(5);
    h_E_pi_theta->Draw("COLZ");
    h_E_pi_theta->SetTitle("E vs theta pi_LD2");
    h_E_pi_theta->GetXaxis()->SetTitle("theta pion (degrees)");
    MonEnrg.cd(6);
    h_E_pi_phi->Draw("COLZ");
    h_E_pi_phi->SetTitle("E vs phi pi_LD2");
    h_E_pi_phi->GetXaxis()->SetTitle("phi pion (degrees)");

    MonEnrg.Print((filename + ".pdf").c_str());

}

void Monitoring::DrawVertexHistograms(const std::string filename){
    //if you want a plot of vZ for Sn and LD2, you need to call this function twice with different filenames and monLD2/monSn
    //there is no way to have both Vz of both nuclei in the same canvas
    TCanvas MonVtz("Monitoring Vtxcanvas", "Monitoring VtxHistograms");
    MonVtz.Divide(2, 2);
    MonVtz.cd(1);
    h_vertexZ->Draw("hist");
    h_vertexZ->SetTitle("Vz_el Distribution");
    h_vertexZ->GetXaxis()->SetTitle("V_{z} (cm)");
    h_vertexZ->GetXaxis()->SetRangeUser(-20,10);
    MonVtz.cd(2);
    h_vertexZ_pi->Draw("hist");
    h_vertexZ_pi->SetTitle("Vz_pi Distribution");
    h_vertexZ_pi->GetXaxis()->SetTitle("V_{z} (cm)");
    MonVtz.cd(3);
    h_DeltaVz->Draw("hist");
    h_DeltaVz->SetTitle("DeltaVz Distribution");
    h_DeltaVz->GetXaxis()->SetTitle("DeltaV_{z} (cm)");


    MonVtz.Print((filename + ".pdf").c_str());

}

void Monitoring::DrawHelicityHistograms(const std::string filename){
    TCanvas MonHel("Monitoring Heli", "Monitoring Heli");
    MonHel.Divide(2, 2);
    MonHel.cd(1);
    //h_helicity->Draw("hist");
    //h_helicity->SetTitle("Helicity Distribution");
    //h_helicity->GetXaxis()->SetTitle("Helicity");
    h_phih_plus->Draw("hist");
    h_phih_plus->SetTitle("phih Hel + Distribution");
    h_phih_plus->GetXaxis()->SetTitle("phih (degrees)");

    MonHel.cd(2);
    //h_helicity_raw->Draw("hist");
    //h_helicity_raw->SetTitle("Raw Helicity Distribution");
    //h_helicity_raw->GetXaxis()->SetTitle("Raw Helicity");
    h_phih_minus->Draw("hist");
    h_phih_minus->SetTitle("phih Hel - Distribution");
    h_phih_minus->GetXaxis()->SetTitle("phih (degrees)");

    MonHel.cd(3);
    MonHel.cd(4);

    h_BSA = h_phih_plus->GetAsymmetry(h_phih_minus); 
    h_BSA->Draw();

    MonHel.Print((filename + ".pdf").c_str());


}


void Monitoring::DrawCompRECMC(const std::string filename){
    TCanvas MonComp("MonComp", "MonComp");
    MonComp.Divide(2, 2);
    MonComp.cd(1);
    h_Q2comp->Draw("COLZ");
    MonComp.cd(3);
    h_Q2->Draw("hist");
    MonComp.cd(4);
    h_Q2MC->Draw("hist");
    
    MonComp.Print((filename + ".pdf").c_str());


}

void Monitoring::DrawMonSratio(const std::string filename){
    TCanvas MonS("MonS", "MonS");
    MonS.Divide(2, 2);
    MonS.cd(1);
    h_nu->Draw("hist");
    MonS.cd(2);
    h_Q2->Draw("hist");
    MonS.cd(3);
    h_xb->Draw("hist");
    MonS.cd(4);
    h_z->Draw("hist");
    MonS.Print((filename + ".pdf").c_str());

}


void Monitoring::SaveHistRoot(const std::string& filenameREC) {
    TFile* rootFile = new TFile((filenameREC + ".root").c_str(), "RECREATE");
    h_vertexZ->Write();
    h_Q2->Write();
    h_xb->Write();
    h_y->Write();
    h_nu->Write();
    h_W2->Write();
    h_z->Write();
    h_pt2->Write();
    h_phih->Write();
    h_pt2z->Write();

    h_calXY->Write();
    h_lu->Write();
    h_lv->Write();
    h_lw->Write();
    h_epcal->Write();
    h_Nphe15->Write();
    h_Nphe16->Write();
    h_calSector->Write();
    h_helicity->Write();
    h_helicity_raw->Write();
    h_px_el->Write();
    h_py_el->Write();
    h_pz_el->Write();
    h_ptot_el->Write();     //this
    h_theta_el->Write();
    h_phi_el->Write();
    h_polcoord_el->Write();
    h_E_el->Write();
    h_E_el_theta->Write();
    h_E_el_phi->Write();
    h_px_pi->Write();
    h_py_pi->Write();
    h_pz_pi->Write();
    h_ptot_pi->Write();     //this
    h_theta_pi->Write();    //this
    h_phi_pi->Write();      //this
    h_polcoord_pi->Write();
    h_E_pi->Write();
    h_E_pi_theta->Write();
    h_E_pi_phi->Write();
    h_vertexZ_pi->Write();
    h_DeltaVz->Write();

    h_chi2_el->Write();
    h_chi2_pi->Write();
    h_chi2_pid_pi->Write();
    h_luthetael->Write();
    h_sampl_el->Write();


    rootFile->Close();
    delete rootFile;


    //counters are here to monitor cut evolution 
    //std::cout << " counterprecuts = " <<  counterprecuts << std::endl;
    //std::cout << " counterpassVz = " <<  counterpassVz << std::endl;
    //std::cout << " counterpassQ2 = " <<  counterpassQ2 << std::endl;
    //std::cout << " counterpassy = " <<  counterpassy << std::endl;
    //std::cout << " counterpassnu = " <<  counterpassnu << std::endl;
    //std::cout << " counterpassw = " <<  counterpassw << std::endl;
    //std::cout << " counterpassz = " <<  counterpassz << std::endl;
    //std::cout << " counterpasspt2 = " <<  counterpasspt2 << std::endl;
    //std::cout << " counterpasshadid = " <<  counterpasshadid << std::endl;


}
//void Monitoring::FillHistogramsNoCuts( const Event& event){
//    int targetType = event.GetTargetType();//using a flag for targets 
/*
    if (targetType == 0 && cutd.PassCutsElectrons(event)==true) {
        
        //set a counter that increases when electroncuts = passed; in order for it to be called when R is  computed in had variables (?) TBD
        h_nuD->Fill(event.Getnu());
        for (const Particle& hadron : event.GetHadrons()) {
            if (cutd.PassCutsHadrons(hadron)==true){
                ////not using the if (==false) return statement 

                h_nu_z_pt2D->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
            }
        }
    }
    else if (targetType == 1 && cuta.PassCutsElectrons(event)==true) {
        counter_elSn++; //counter for only electrons for z and pt
        //here change the else if to just else in order to have a generic target 
        h_nuA->Fill(event.Getnu());
        for (const Particle& hadron : event.GetHadrons()) {
            if (cuta.PassCutsHadrons(hadron)==true){
                h_nu_z_pt2A->Fill(event.Getnu(), hadron.Getz(), hadron.Getpt2());
            }
        }
        //Add Here Cu && change cut for Sn here. 1st cut is passelectrons 
    }

    //Add CxC
*/
//}


Monitoring::~Monitoring() {
    // Delete dynamically allocated histograms
    delete h_Q2;
    delete h_xb;
    delete h_y;
    delete h_nu;
    delete h_W2;
    delete h_z;
    delete h_pt2;
    delete h_phih;
    delete h_vertexZ;
    delete h_pid;
    //delete h_Q2posel;
    //delete h_nuposel;
    //delete h_Q2_had;
    //delete h_nu_had;
    //delete h_z_had;
    //delete h_pt2_had;
}
