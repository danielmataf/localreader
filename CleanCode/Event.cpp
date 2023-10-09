#include <TLorentzVector.h>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include "reader.h"
#include "Event.h"
#include "Particle.h"
#include "constants.h"
    
Event::Event() : electron(TLorentzVector(0.0, 0.0, 0.0, 0.0), 0) {
    m_D = Constants::MASS_NEUTRON + Constants::MASS_PROTON ;
}

void Event::AddElectron(const TLorentzVector& electronMomentum) {
        electron = Particle(electronMomentum, 11);  
        electron.SetMomentum(electronMomentum);
}



void Event::AddHadron(const TLorentzVector& hadronMomentum, int pid) {
    hadrons.push_back(Particle(hadronMomentum,  pid));
}

int Event::GetEventIndex() const {
    return eventIndex;
}

const std::vector<Particle>& Event::GetHadrons() const {
    return hadrons;
}

void Event::Print() {   //add int v=0 as argument 4 different types of verbose TBD

    std::cout << "Electrons:" << std::endl;
    //for (const Particle& electron : event.GetElectrons()) {
        std::cout << "  Particle ID: " << electron.GetPID() << std::endl;
        std::cout << "  Momentum: (" << electron.GetMomentum().Px() << ", "
                  << electron.GetMomentum().Py() << ", " << electron.GetMomentum().Pz() << ")" << std::endl;
        std::cout << " Total Momentum: " << electron.GetMomentum().P() << std::endl;
        std::cout << " Q2 value : " << Q2<< std::endl;
        std::cout << " nu value : " << nu<< std::endl;
        std::cout << " y  value : " << y<< std::endl;
        std::cout << " W2 value : " << w2<< std::endl;
        std::cout << " xb value : " << xb<< std::endl;

    
    //}
    std::cout << "Hadrons:" << std::endl;
    for (const Particle& hadron : GetHadrons()) {
        //if (hadron.GetPID()== 211){
            std::cout << "  ____________ " <<  std::endl;
            std::cout << "  Particle ID: " << hadron.GetPID() << std::endl;
            std::cout << " Total Momentum: " << hadron.GetMomentum().P() << std::endl;
            std::cout << " z value : " << hadron.Getz()<< std::endl;
            std::cout << " pt2 value : " << hadron.Getpt2()<< std::endl;
            std::cout << " phih value : " << hadron.Getphih()<< std::endl;

        //}
    }
}

int Event::CalcKinematics(){
    double theta_e = acos(Constants::elBeam.Vect().Dot(electron.GetMomentum().Vect()) / (Constants::elBeam.Vect().Mag() * electron.GetMomentum().Vect().Mag()));
    Q2 = 4 * Constants::elBeam.E() * electron.GetMomentum().E() * pow(sin(theta_e / 2), 2);
    nu  = Constants::elBeam.E() - electron.GetMomentum().E();
    y = nu / Constants::elBeam.E(); 
    w2 =  sqrt(m_D* m_D + 2 * m_D * nu - Q2);
    xb = Q2 / (2 * m_D * nu);
    electron.SetKinVariables(Q2, nu, y, w2, xb);    //useless 
    return 0;

}

//void SetKinVariables(double , double  , double , double, double);


void Event::calcAll(){
    CalcKinematics( );
    for (Particle& hadron : hadrons) {
        hadron.CalcHadronKin(electron.GetMomentum(), hadron.GetMomentum());
        
    }
}

double Event::GetQ2() const {
        return Q2;
}
double Event::Getnu() const {
    return nu;
}
double Event::Gety() const {
    return y;
} 
double Event::GetW2() const{
    return w2;
}
double Event::Getxb() const{
    return xb;
}
//double Event::Getz()     const {
//    return hadron.GetZ(); 
//}
//double Event::Getpt2()   const {
//    return hadron.Getpt2() ; 
//}
//double Event::Getphih()  const {
//    return hadron.Getphih(); 
//}