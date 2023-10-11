#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include "Event.h" 
#include "Monitoring.h"
#include "CutSet.h"

    //int nubin = 100;  
    int QminX = 0;
    int QmaxX = 6;
    int nubin = 100;

        int phibin= 10;
    int Rbin = 10;
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
  
Monitoring::Monitoring(CutSet a):   //h_Q2("Q2", "Q2", nubin, QminX, QmaxX),
                                    h_xb("xb", "xb", nubin, xminX, xmaxX),
                                    h_y ("y" , "y" , nubin, yminX, ymaxX),
                                    h_nu("nu", "nu", nubin,numinX,numaxX),
                                    h_W2("W2", "W2", nubin, WminX, WmaxX),
                                    h_z    ("z", "z", nubin, zminX, zmaxX),
                                    h_pt2  ("pt2", "pt2", nubin, pt2minX, pt2maxX),
                                    h_phih ("phih", "phih", nubin, phihminX, phihmaxX)
                                    //add histos here if needed
                                    //outputFile("monitoring_output.root", "RECREATE")


  {
    cut1 = a;
    
    //create the root file 
    //outputFile = new TFile("histQ.root", "RECREATE");

}



void Monitoring::FillHistograms(const Event& event) {
    if (cut1.PassCutsElectrons(event)==false) return;
    std::cout << " -> electron cuts passed successfully  " << std::endl;
    h_Q2->Fill(event.GetQ2());
    h_xb.Fill(event.Getxb());
    h_y.Fill(event.Gety());
    h_nu.Fill(event.Getnu());
    h_W2.Fill(event.GetW2());
    std::cout << " kinel " << event.GetQ2()<<" , " << event.Getxb()<<" , " << event.Gety()<<" , "<< event.Getnu()<<" , "<< event.GetW2()<<std::endl;  
    for (const Particle& hadron : event.GetHadrons()) {
        if (cut1.PassCutsHadrons(hadron)==true){
            std::cout << " -> hadron cuts passed successfully  " << std::endl;
            //if passcuts(given hadron )
            //fill histogram with hadron variables  
            h_z.Fill(hadron.Getz());
            h_pt2.Fill(hadron.Getpt2());
            h_phih.Fill(hadron.Getphih());
            std::cout << " kinhad " << hadron.Getz()<<" , " << hadron.Getpt2()<<" , " << hadron.Getphih()<<std::endl;  

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
    h_xb.Write();
    h_y.Write();
    h_nu.Write();
    h_W2.Write();
    h_z.Write();
    h_pt2.Write();
    h_phih.Write();
    //file.Write();
    // Write other histograms here
    file.Close();
}

