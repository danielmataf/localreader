#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include <TPDF.h>
#include "Event.h" 
#include "Monitoring.h"
#include "CutSet.h"

    //int nubin = 100;  
  
Monitoring::Monitoring(CutSet a)//:   h_Q2("Q2", "Q2", nubin, QminX, QmaxX),
                                    //h_xb("xb", "xb", nubin, xminX, xmaxX),
                                    //h_y ("y" , "y" , nubin, yminX, ymaxX),
                                    //h_nu("nu", "nu", nubin,numinX,numaxX),
                                    //h_W2("W2", "W2", nubin, WminX, WmaxX),
                                    //h_z    ("z", "z", nubin, zminX, zmaxX),
                                    //h_pt2  ("pt2", "pt2", nubin, pt2minX, pt2maxX),
                                    //h_phih ("phih", "phih", nubin, phihminX, phihmaxX)
                                    //add histos here if needed
                                    //outputFile("monitoring_output.root", "RECREATE")


  {
    cut1 = a;
    
    //create the root file 
    //outputFile = new TFile("histQ.root", "RECREATE");

}

//h_xbpre
//h_ypre
//h_nupre
//h_W2pre
//h_zpre
//h_pt2pre
//h_phihpre

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
    //
    //std::cout << " -> electron cuts passed successfully  " << std::endl;
    //
    counterel_R ++;
    h_vertexZ->Fill(event.GetVz()); //Vz only exists when an electron is detecte !!!!
    h_Q2->Fill(event.GetQ2());
    h_xb->Fill(event.Getxb());
    h_y->Fill(event.Gety());
    h_nu->Fill(event.Getnu());
    h_W2->Fill(event.GetW2());
    //
    //std::cout << " kinel " << event.GetQ2()<<" , " << event.Getxb()<<" , " << event.Gety()<<" , "<< event.Getnu()<<" , "<< event.GetW2()<<std::endl;  
    //
    for (const Particle& hadron : event.GetHadrons()) {
        if (cut1.PassCutsHadrons(hadron)==true){
            //
            //std::cout << " -> hadron cuts passed successfully  " << std::endl;
            //
            //if passcuts(given hadron )
            //fill histogram with hadron variables  
            h_z->Fill(hadron.Getz());
            h_pt2->Fill(hadron.Getpt2());
            h_phih->Fill(hadron.Getphih());
            //
            //std::cout << " kinhad " << hadron.Getz()<<" , " << hadron.Getpt2()<<" , " << hadron.Getphih()<<std::endl;  
            //
        }
    }
 //add histos 

}

void Monitoring::WriteHistogramsToFile(const std::string filename) {
    //this function recreates a new rootfile everytime is called 
    //useful to have different rootfiles if different cuts were implemented
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
