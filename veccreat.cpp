#include <TArrow.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TLorentzVector.h>

#define _USE_MATH_DEFINES
#include <math.h>



//TLorentzVector *v_MC_scapit  = new TLorentzVector;

//function that receives ...........? as input, and has to read the hipo data bank to getfloats to fill a TLorentzVector
//returns a TLorentzvector with Px,y,z and mass 



TLorentzVector fillVector(TLorentzVector* vec1,double mass_particle, double evt) {
    // Create a TLorentzVector object with the given components
    //TLorentzVector myVector(px, py, pz, E);
    
    vec1->SetPx(REC_Particle_bank.getFloat("px", evt));
    vec1->SetPy(REC_Particle_bank.getFloat("py", evt));
    vec1->SetPz(REC_Particle_bank.getFloat("pz", evt));
    TVector3 p_vec1 = vec1->Vect();
    vec1->SetVectM(p_vec1,mass_particle);
    return vec1;
}


    

TLorentzVector fillVector2(TLorentzVector* vec1,double mass_particle, double evt) {
    // Create a TLorentzVector object with the given components
    //TLorentzVector myVector(px, py, pz, E);
    double Px = REC_Particle_bank.getFloat("px", evt) 
    double Py = REC_Particle_bank.getFloat("py", evt)
    double Pz = REC_Particle_bank.getFloat("pz", evt)
    double E = std::sqrt(Px*Px + Py*Py + Pz*Pz + mass_particle*mass_particle);
    vec1->SetPxPyPzE(Px, Py, Pz, E);
    return vec1;
}