#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include <TPDF.h>
#include "Event.h" 
#include "Monitoring.h"
#include "CutSet.h"
#include "constants.h"

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
      h_vertexZ(new TH1F(("targetVz_" + targetName).c_str(), "vertex4target", 100, -40, 40)),
      h_pid(new TH1F(("pid_" + targetName).c_str(), "pid", 100, -250, 250)),
      h_xQ2(new TH2F(("xQ2_" + targetName).c_str(), "xQ2", nubin, xminX, xmaxX, nubin, QminX, QmaxX)),
      h_xQ2pos(new TH2F(("xQ2pos_" + targetName).c_str(), "xQ2pos", nubin, xminX, xmaxX, nubin, QminX, QmaxX)),
      h_px_el(new TH1F(("px_ele_" + targetName).c_str(), "px_ele", nubin, 0, 10)),
      h_py_el(new TH1F(("py_ele_" + targetName).c_str(), "py_ele", nubin, 0, 10)),
      h_pz_el(new TH1F(("pz_ele_" + targetName).c_str(), "pz_ele", nubin, 0, 10)),
      h_px_pi(new TH1F(("px_pro_" + targetName).c_str(), "px_pro", nubin, 0, 10)),
      h_py_pi(new TH1F(("py_pro_" + targetName).c_str(), "py_pro", nubin, 0, 10)),
      h_pz_pi(new TH1F(("pz_pro_" + targetName).c_str(), "pz_pro", nubin, 0, 10)),
      h_lu(new TH1F(("lu_el" + targetName).c_str(), "lu", nubin, 0, 400)),
      h_lv(new TH1F(("lv_el" + targetName).c_str(), "lv", nubin, 0, 400)),
      h_lw(new TH1F(("lw_el" + targetName).c_str(), "lw", nubin, 0, 400)),
      h_epcal(new TH1F(("epcal_el" + targetName).c_str(), "epcal", nubin, 0, 1)),
      h_eecalin(new TH1F(("eecalin_el" + targetName).c_str(), "eecalin", nubin, 0, 1)),
      h_epcalout(new TH1F(("epcalout_el" + targetName).c_str(), "epcalout", nubin, 0, 1)),
      h_calXY(new TH2F(("calxy_el" + targetName).c_str(), "calxy", nubin, -400, 400, nubin, -400, 400)),
      h_calEall(new TH2F(("calEall_el" + targetName).c_str(), "calEall", nubin, 0, 1, nubin, 0, 1)),
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

void Monitoring::FillHistograms(const Event& event) {
    //if (cut1.PassCutLD2Target(event)==false) return;
    //h_xQ2->Fill(event.Getxb(), event.GetQ2());

    if (cut1.PassCutsElectrons(event)==false) return;
    //std::cout << " -> electron cuts passed successfully  " << std::endl;
    counterel_R ++;
    h_vertexZ->Fill(event.GetVz());     //Vz only exists when an electron is detected !!!!
    h_Q2->Fill(event.GetQ2());
    h_xb->Fill(event.Getxb());
    h_y->Fill(event.Gety());    
    h_nu->Fill(event.Getnu());
    h_W2->Fill(event.GetW2());
    h_calXY->Fill(event.GetCalX(), event.GetCalY());
    h_lu->Fill(event.Getlu());
    h_lv->Fill(event.Getlv());
    h_lw->Fill(event.Getlw());
    h_epcal->Fill(event.GetEpcal());
    h_eecalin->Fill(event.GetEcalin());
    h_epcalout->Fill(event.GetEcalout());
    h_calEall->Fill(event.GetEpcal(), event.GetEcalin()+event.GetEcalout());

    
    //h_xQ2pos->Fill(event.Getxb(), event.GetQ2());

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
            //
            //std::cout << " kinhad " << hadron.Getz()<<" , " << hadron.Getpt2()<<" , " << hadron.Getphih()<<std::endl;  
            //
            }
        }
    }
 //add histos 

}


void Monitoring::Fill_R_Histograms(const Event& event, const std::string target) {
    if (target == "D" && cut1.PassCutsElectrons(event)==true) {
        R_nu_el->Fill(event.Getnu());
        for (const Particle& hadron : event.GetHadrons()) {
            if (cut1.PassCutsHadrons(hadron)==true){
                R_nu_had->Fill(event.Getnu());
                R_z->Fill(hadron.Getz());
                R_pt2->Fill(hadron.Getpt2());
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
}

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


void Monitoring::DrawHistograms(const std::string filename) {
    TCanvas MonC("Monitoring canvas", "Monitoring Histograms");
    MonC.Divide(3, 3);
    MonC.cd(1);
    h_Q2->Draw("hist");
    h_Q2->SetTitle("Q2 Distribution");
    MonC.cd(2);
    h_xb->Draw("hist");
    h_xb->SetTitle("xb Distribution");
    MonC.cd(3);
    h_y->Draw("hist");
    h_y->SetTitle("y Distribution");
    MonC.cd(4);
    h_nu->Draw("hist");
    h_nu->SetTitle("nu Distribution");
    MonC.cd(5);
    h_W2->Draw("hist");
    h_W2->SetTitle("W2 Distribution");
    MonC.cd(6);
    h_z->Draw("hist");
    h_z->SetTitle("z Distribution");
    MonC.cd(6)->SetLogy();  // Set logarithmic scale on y-axis for h_z
    //add log scale on y axis to z distribution
    MonC.cd(7);
    h_pt2->Draw("hist");
    h_pt2->SetTitle("pt2 Distribution");
    MonC.cd(7)->SetLogy();  // Set logarithmic scale on y-axis for h_z
    //add log scale on y axis to pt distribution



    MonC.cd(8);
    h_phih->Draw("hist");
    h_phih->SetTitle("phih Distribution");
    MonC.cd(9);
    h_vertexZ->Draw("hist");
    h_vertexZ->SetTitle("Vertex Z Distribution");


    MonC.Print((filename + ".pdf").c_str());
}

void Monitoring::DrawHistogramsPos(const std::string filename) {
    TCanvas MonC("Monitoring canvas", "Monitoring Histograms");
    MonC.Divide(2, 2);
    MonC.cd(1);
    h_xQ2->Draw("COLZ");
    h_xQ2->SetTitle("x vs Q2 Distribution (LD2)");
    h_xQ2->GetXaxis()->SetTitle("x_{B}");
    h_xQ2->GetYaxis()->SetTitle("Q^{2}");
    MonC.cd(2);
    h_xQ2pos->Draw("COLZ");
    h_xQ2pos->SetTitle("x vs Q2 Distribution (LD2)");
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
    MonCal.cd(2);
    h_lv->Draw("hist");
    h_lv->SetTitle("lv el_LD2");
    MonCal.cd(3);
    h_lw->Draw("hist");
    h_lw->SetTitle("lw el_LD2");
    MonCal.cd(4);
    h_epcal->Draw("hist");
    h_epcal->SetTitle("epcal el_LD2");
    MonCal.cd(5);
    h_eecalin->Draw("hist");
    h_eecalin->SetTitle("Ecal_inner el_LD2");
    MonCal.cd(6);
    h_epcalout->Draw("hist");
    h_epcalout->SetTitle("Ecal_outer el_LD2");
    MonCal.cd(7);
    h_calEall->Draw("COLZ");
    h_calEall->SetTitle("Ecal vs Epcal el_LD2");
    MonCal.cd(8);
    h_calXY->Draw("COLZ");
    h_calXY->SetTitle("PCalXY");


    MonCal.Print((filename + ".pdf").c_str());
}



void Monitoring::FillHistogramsNoCuts( const Event& event){
    int targetType = event.GetTargetType();//using a flag for targets 
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
}


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