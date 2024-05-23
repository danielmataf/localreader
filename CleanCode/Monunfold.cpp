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

    h_Q2(new TH1F(("U_Q2_" + targetName).c_str(), "Q2", nubin, QminX, QmaxX)),
    h_xb(new TH1F(("U_xb_" + targetName).c_str(), "xb", nubin, xminX, xmaxX)),
    h_y(new TH1F(("U_y_" + targetName).c_str(), "y", nubin, yminX, ymaxX)),
    h_nu(new TH1F(("U_nu_" + targetName).c_str(), "nu", nubin, numinX, numaxX)),
    h_W2(new TH1F(("U_W2_" + targetName).c_str(), "W2", nubin, WminX, WmaxX)),
    h_z(new TH1F(("U_z_" + targetName).c_str(), "z", nubin, zminX, zmaxX)),
    h_pt2(new TH1F(("U_pt2_" + targetName).c_str(), "pt2", nubin, pt2minX, pt2maxX)),
    h_phih(new TH1F(("U_phih_" + targetName).c_str(), "phih", nubin, phihminX, phihmaxX)),
    h_vertexZ(new TH1F(("U_targetVz_" + targetName).c_str(), "vertex4target", 100, -10, 0)),
    
    h_Q2MC(new TH1F(("U_Q2_MC" + targetName).c_str(), "Q2MC", nubin, QminX, QmaxX)),
    h_xbMC(new TH1F(("U_xb_MC" + targetName).c_str(), "xbMC", nubin, xminX, xmaxX)),
    h_yMC(new TH1F(("U_y_MC" + targetName).c_str(), "yMC", nubin, yminX, ymaxX)),
    h_nuMC(new TH1F(("U_nu_MC" + targetName).c_str(), "nuMC", nubin, numinX, numaxX)),
    h_W2MC(new TH1F(("U_W2_MC" + targetName).c_str(), "W2MC", nubin, WminX, WmaxX)),
    h_zMC(new TH1F(("U_z_MC" + targetName).c_str(), "zMC", nubin, zminX, zmaxX)),
    h_pt2MC(new TH1F(("U_pt2_MC" + targetName).c_str(), "pt2MC", nubin, pt2minX, pt2maxX)),
    h_phihMC(new TH1F(("U_phih_MC" + targetName).c_str(), "phihMC", nubin, phihminX, phihmaxX)),
    h_vertexZMC(new TH1F(("U_targetVz_MC" + targetName).c_str(), "vertex4targetMC", 100, -20, 10)),

    h_Q2comp(new TH2F(("U_Q2comp_" + targetName).c_str(), "Q2comp", nubin, QminX, QmaxX, nubin, QminX, QmaxX)),
    h_xbcomp(new TH2F(("U_xbcomp_" + targetName).c_str(), "xbcomp", nubin, xminX, xmaxX, nubin, xminX, xmaxX)),
    h_ycomp(new TH2F(("U_ycomp_" + targetName).c_str(), "ycomp", nubin, yminX, ymaxX, nubin, yminX, ymaxX)),
    h_nucomp(new TH2F(("U_nucomp_" + targetName).c_str(), "nucomp", nubin, numinX, numaxX, nubin, numinX, numaxX)),
    h_W2comp(new TH2F(("U_W2comp_" + targetName).c_str(), "W2comp", nubin, WminX, 15, nubin, WminX, 15)),
    h_zcomp(new TH2F(("U_zcomp_" + targetName).c_str(), "zcomp", nubin, zminX, zmaxX, nubin, zminX, zmaxX)),
    h_pt2comp(new TH2F(("U_pt2comp_" + targetName).c_str(), "pt2comp", nubin, pt2minX, pt2maxX, nubin, pt2minX, pt2maxX)),
    h_phihcomp(new TH2F(("U_phihcomp_" + targetName).c_str(), "phihcomp", nubin, phihminX, phihmaxX, nubin, phihminX, phihmaxX)),
    h_vertexZcomp(new TH2F(("U_targetVzcomp_" + targetName).c_str(), "vertex4targetcomp", 100, -10, 10.0, 100, -10, 10.0)),




      counterel_R(0) {
    // Add more histograms as needed
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
    c1.Divide(1, 1);
    c1.cd(1);
    h_vertexZ->Draw();
    munfoldC2.h_vertexZ->SetLineColor(kBlack);
    munfoldC2.h_vertexZ->Draw("same");
    //h_vertexZMC->SetLineColor(kRed);
    //h_vertexZMC->Draw("same");
    //h_vertexZcomp->Draw("colz");
    c1.Print((filename + ".pdf").c_str());
    //then create a function to draw. Delete this function once used
}



void Monunfold::Fill(const Event& event) {
    // Fill histograms with event data
    h_Q2->Fill(event.GetQ2());
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
        z_values.push_back(hadron.Getz());
        pt2_values.push_back(hadron.Getpt2());
        phih_values.push_back(hadron.Getphih());

    }
    
    //2nd loop over MChadrons
    for (const Particle& MChadron : event.GetMCHadrons()) {
        zMC_values.push_back(MChadron.GetzMC());
        pt2MC_values.push_back(MChadron.Getpt2MC());
        phihMC_values.push_back(MChadron.GetphihMC());
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

void Monunfold::DrawHistograms(const std::string filename) {
    // Draw histograms and save them to a PDF file
    TCanvas c1("c1", "c1", 800, 600);
    c1.Divide(3, 3);
    c1.cd(1);
    h_Q2->Draw();
    h_Q2MC->SetLineColor(kRed);
    h_Q2MC->Draw("same");
    c1.cd(2);
    h_xb->Draw();
    h_xbMC->SetLineColor(kRed);
    h_xbMC->Draw("same");
    c1.cd(3);
    h_y->Draw();
    h_yMC->SetLineColor(kRed);
    h_yMC->Draw("same");
    c1.cd(4);
    h_nu->Draw();
    h_nuMC->SetLineColor(kRed);
    h_nuMC->Draw("same");
    c1.cd(5);
    h_W2->Draw();
    h_W2MC->SetLineColor(kRed);
    h_W2MC->Draw("same");
    c1.cd(6);
    h_z->Draw();
    h_zMC->SetLineColor(kRed);
    h_zMC->Draw("same");

    c1.cd(7);
    h_pt2->Draw();
    h_pt2MC->SetLineColor(kRed);
    h_pt2MC->Draw("same");
    c1.cd(8);
    h_phih->Draw();
    h_phihMC->SetLineColor(kRed);
    h_phihMC->Draw("same");
    c1.cd(9);
    h_vertexZ->Draw();
    h_vertexZMC->SetLineColor(kRed);
    h_vertexZMC->Draw();
    c1.Print((filename + ".pdf").c_str());

}

void Monunfold::DrawCompRECMC(const std::string filename){
    TCanvas c2("c2", "c2", 800, 600);
    c2.Divide(3, 3);
    c2.cd(1);
    h_Q2comp->Draw("colz");
    c2.cd(2);
    h_xbcomp->Draw("colz");
    c2.cd(3);
    h_ycomp->Draw("colz");
    c2.cd(4);
    h_nucomp->Draw("colz");
    c2.cd(5);
    h_W2comp->Draw("colz");
    c2.cd(6);
    h_zcomp->Draw("colz");
    c2.cd(7);
    h_pt2comp->Draw("colz");
    c2.cd(8);
    h_phihcomp->Draw("colz");
    c2.cd(9);
    h_vertexZcomp->Draw("colz");
    c2.Print((filename + ".pdf").c_str());
}



Monunfold::~Monunfold() {
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
    
}
