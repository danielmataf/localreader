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
    if (cut1.PassCutsElectrons(event)==false) return;
    //std::cout << " -> electron cuts passed successfully  " << std::endl;
    counterel_R ++;
    h_vertexZ->Fill(event.GetVz());     //Vz only exists when an electron is detecte !!!!
    h_Q2->Fill(event.GetQ2());
    h_xb->Fill(event.Getxb());
    h_y->Fill(event.Gety());
    h_nu->Fill(event.Getnu());
    h_W2->Fill(event.GetW2());
    //h_Q2posel->Fill(event.GetQ2());
    //h_nuposel->Fill(event.Getnu());     //throw hists in mratio.cpp(filling) TBD
    //std::cout << " kinel " << event.GetQ2()<<" , " << event.Getxb()<<" , " << event.Gety()<<" , "<< event.Getnu()<<" , "<< event.GetW2()<<std::endl;  
    for (const Particle& hadron : event.GetHadrons()) {
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
    MonC.cd(7);
    h_pt2->Draw("hist");
    h_pt2->SetTitle("pt2 Distribution");
    MonC.cd(8);
    h_phih->Draw("hist");
    h_phih->SetTitle("phih Distribution");
    MonC.cd(9);
    h_vertexZ->Draw("hist");
    h_vertexZ->SetTitle("Vertex Z Distribution");


    MonC.Print((filename + ".pdf").c_str());
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