#include <TLorentzVector.h>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include "reader.h"
#include "Event.h"
#include "Particle.h"

    
Event::Event() : electron(TLorentzVector(0.0, 0.0, 0.0, 0.0), 0) {
    IncidentLepton.SetPxPyPzE(0.0, 0.0, 11.0, 11.0);
    m_D = m_n + m_p;
}

void Event::AddElectron(const TLorentzVector& electronMomentum) {
        electron = Particle(electronMomentum, 11);  // Update the electron member
        electron.SetMomentum(electronMomentum);
        //CalcKinematics( );
}

//set electron variables in event
//and pion variables in event 


void Event::AddHadron(const TLorentzVector& hadronMomentum, int pid) {
    hadrons.push_back(Particle(hadronMomentum,  pid));
    hadrons.back().CalcHadronKin(electron.GetMomentum(),hadronMomentum);
        std::cout << "-------------------testing electron momentum "<< electron.GetMomentum().P() <<std::endl;

        // ERROR IS IN THE  EventReader::ProcessEvent FUNCTION
        // ELECTRON MOMENTA ARE SET AFTER THE HADRON MOMENTA
        // THE FUNCTION CALLS THE Event::AddHadron BEFORE CALLING Event::AddElectron
        // THEREFORE ELECTRON MOMENTUM IS EMTPY FOR ANY HADRON ALWAYS!!!! 
        //              Fix this 2morrow 
        //change the Addelectron function and with it make a vector of particles as before.
        // then change the method of electron selection. Cant be anymore in the ProcessEvent  



    //call the funbction calchadronkin 
    //need to store the electronmomentum in the particle so it can be accesed trhoug electron.something
    //probably electron.GeteMomentum (?) that retunrs a vector and a fct SeteMomentum has to be called in AddElectron I suppose 
}

int Event::GetEventIndex() const {
    return eventIndex;
}

//const std::vector<Particle>& Event::GetElectrons() const {
//    return electrons;
//}

const std::vector<Particle>& Event::GetHadrons() const {
    return hadrons;
}

void Event::Print() {   //add int v=0 as argument 4 different types of verbose TBD
    //std::cout << "Event Index: " << event.GetEventIndex() << std::endl;

    std::cout << "Electrons:" << std::endl;
    //for (const Particle& electron : event.GetElectrons()) {
        std::cout << "  Particle ID: " << electron.GetPID() << std::endl;
        //std::cout << "  Momentum: (" << electron.GetMomentum().Px() << ", "
        //          << electron.GetMomentum().Py() << ", " << electron.GetMomentum().Pz() << ")" << std::endl;
        std::cout << " Total Momentum: " << electron.GetMomentum().P() << std::endl;
        std::cout << " Q2 value : " << electron.GetQ2()<< std::endl;
        std::cout << " nu value : " << electron.Getnu()<< std::endl;
        std::cout << " y  value : " << electron.Gety()<< std::endl;
        std::cout << " W2 value : " << electron.GetW2()<< std::endl;
        std::cout << " xb value : " << electron.Getxb()<< std::endl;

    
    //}

    // Print hadrons' momentum
    std::cout << "Hadrons:" << std::endl;
    for (const Particle& hadron : GetHadrons()) {
        //if (hadron.GetPID()== 211){
            std::cout << "  ____________ " <<  std::endl;
            std::cout << "  Particle ID: " << hadron.GetPID() << std::endl;
            //std::cout << "  Momentum: (" << hadron.GetMomentum().Px() << ", "
            //          << hadron.GetMomentum().Py() << ", " << hadron.GetMomentum().Pz() << ")" << std::endl;
            std::cout << " Total Momentum: " << hadron.GetMomentum().P() << std::endl;
            std::cout << " z value : " << hadron.Getz()<< std::endl;
            std::cout << " pt2 value : " << hadron.Getpt2()<< std::endl;
            std::cout << " phih value : " << hadron.Getphih()<< std::endl;

        //}
    }
}

int Event::CalcKinematics(){
    //only argument is the scattered electron vector
    //!!!!!Argument was changed to no argument importing vector of scattered lepton 

    double theta_e = acos(IncidentLepton.Vect().Dot(electron.GetMomentum().Vect()) / (IncidentLepton.Vect().Mag() * electron.GetMomentum().Vect().Mag()));
    Q2 = 4 * IncidentLepton.E() * electron.GetMomentum().E() * pow(sin(theta_e / 2), 2);
    nu  = IncidentLepton.E() - electron.GetMomentum().E();
    y = nu / IncidentLepton.E(); 
    w2 =  sqrt(m_D* m_D + 2 * m_D * nu - Q2);
    xb = Q2 / (2 * m_D * nu);
    electron.SetKinVariables(Q2, nu, y, w2, xb);
    //add option for the target mass 
    return 0;

}



//Create fct calcs kinematix 
//assign elect to variables (same w/ hadron + le variables)
//calc here, stock  them in class particle 


/*
double Event::cordtheta(TLorentzVector* vec){
    //argument is scattered lepton 
    double theta_el = vec->Theta();
    return theta_el;
}

double Event::cordphi(TLorentzVector* vec3){
    //function can work for any particle represented with a vector
    double phi_el = vec3->Phi();
    return phi_el;
}
*/