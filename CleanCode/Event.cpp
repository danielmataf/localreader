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
    
Event::Event() : electron(TLorentzVector(0.0, 0.0, 0.0, 0.0), 0, 0, 0.0, 0.0),
                 MCelectron(TLorentzVector(0.0, 0.0, 0.0, 0.0), 0, 0, 0.0,0.0){
    m_D = (Constants::MASS_NEUTRON + Constants::MASS_PROTON)/2 ;

}

void Event::AddElectron(const TLorentzVector& electronMomentum, int row, double vertexZ, double chi2row ) {
        electron = Particle(electronMomentum, 11 , row, vertexZ, chi2row);   
        electron.SetMomentum(electronMomentum);
        
}

void Event::AddMCElectron(const TLorentzVector& electronMomentum, int row, double vertexZ) {
        MCelectron = Particle(electronMomentum, 11 , row, vertexZ, 0.0);    
        //last argument is 0, because the chi2 is not available in the MC data  
        MCelectron.SetMomentum(electronMomentum);
        
}


void Event::AddHadron(const TLorentzVector& hadronMomentum, int pid, int row, double vertexZ,  double chi2row) {
    hadrons.push_back(Particle(hadronMomentum,  pid, row, vertexZ, chi2row));
}

void Event::AddMCHadron(const TLorentzVector& MChadronMomentum, int MCpid, int MCrow, double MCvertexZ) {
    MChadrons.push_back(Particle(MChadronMomentum,  MCpid, MCrow, MCvertexZ, 0.0));
}

int Event::GetEventIndex() const {
    return eventIndex;
}

const std::vector<Particle>& Event::GetHadrons() const {
    return hadrons;
}
const std::vector<Particle>& Event::GetMCHadrons() const {
    return MChadrons;
}
const Particle& Event::GetElectron() const {
    return electron;
}

const Particle& Event::GetElectronMC() const {
    return MCelectron;
}

void Event::SetVertexZ(double vertexz){
    vz = vertexz;
}
void Event::SetVertexZMC(double vertexzMC){
    MCvz = vertexzMC;
}

void Event::SetVertexX(double vertexx){
    vx = vertexx;
}
void Event::SetVertexY(double vertexy){
    vy = vertexy;
}

double Event::GetVz()const {
    return vz;
}
double Event::GetVzMC()const {
    return MCvz;
}

double Event::GetVx()const {
    return vx;
}double Event::GetVy()const {
    return vy;
}

void Event::SetHel(int input_hel){
    hel = input_hel;
}
void Event::SetHelRaw(int input_helraw){
    helraw = input_helraw;
}

void Event::Setevtnbr(int input_evtnbr){
    evtnbr = input_evtnbr;
}


void Event::Print() {   //add int v=0 as argument 4 different types of verbose TBD
    std::cout << "   " << std::endl;
    std::cout << "Electrons:" << std::endl;
    //for (const Particle& electron : event.GetElectrons()) {
        std::cout << "  Particle ID: " << electron.GetPID() << std::endl;
        std::cout << "  Momentum: (" << electron.GetMomentum().Px() << ", "
                  << electron.GetMomentum().Py() << ", " << electron.GetMomentum().Pz() << ")" << std::endl;
        std::cout << "  Polar Angles: (" << electron.GetMomentum().Theta() << ", "
                  << electron.GetMomentum().Phi() << ")" << std::endl;
        std::cout << " Total Momentum: " << electron.GetMomentum().P() << std::endl;
        std::cout << " Q2 value : " << Q2<< std::endl;
        std::cout << " nu value : " << nu<< std::endl;
        std::cout << " y  value : " << y<< std::endl;
        std::cout << " W2 value : " << w2<< std::endl;
        std::cout << " xb value : " << xb<< std::endl;

    
    //}
    //std::cout << "Hadrons:" << std::endl;
    //for (const Particle& hadron : GetHadrons()) {
    //    //if (hadron.GetPID()== 211){
    //        std::cout << "  ____________ " <<  std::endl;
    //        std::cout << "  Particle ID: " << hadron.GetPID() << std::endl;
    //        std::cout << " Total Momentum: " << hadron.GetMomentum().P() << std::endl;
    //        std::cout << " Polar Angles: (" << hadron.GetMomentum().Theta() << ", "
    //                  << hadron.GetMomentum().Phi() << ")" << std::endl;
    //        std::cout << " z value : " << hadron.Getz()<< std::endl;
    //        std::cout << " pt2 value : " << hadron.Getpt2()<< std::endl;
    //        std::cout << " phih value : " << hadron.Getphih()<< std::endl;
//
    //    //}
    //}
}
void Event::PrintMC() {   //add int v=0 as argument 4 different types of verbose TBD

    std::cout << "MC Electrons: -------------------" << std::endl;
    //for (const Particle& electron : event.GetElectrons()) {
        std::cout << "MC  Particle ID: " << MCelectron.GetPID() << std::endl;
        std::cout << "MC  Momentum: (" << MCelectron.GetMomentum().Px() << ", "
                  << MCelectron.GetMomentum().Py() << ", " << MCelectron.GetMomentum().Pz() << ")" << std::endl;
        std::cout << "MC  Polar Angles: (" << MCelectron.GetMomentum().Theta() << ", "
                  << MCelectron.GetMomentum().Phi() << ")" << std::endl;
        std::cout << "MC Total Momentum: " << MCelectron.GetMomentum().P() << std::endl;
        std::cout << "MC Q2 value : " << MCQ2<< std::endl;
        std::cout << "MC nu value : " << MCnu<< std::endl;
        std::cout << "MC y  value : " << MCy<< std::endl;
        std::cout << "MC W2 value : " << MCw2<< std::endl;
        std::cout << "MC xb value : " << MCxb<< std::endl;

    
    //}
    //std::cout << "MC Hadrons:" << std::endl;
    //for (const Particle& MChadron : GetHadrons()) {
    //    //if (hadron.GetPID()== 211){
    //        std::cout << "  ____________ " <<  std::endl;
    //        std::cout << " MC  Particle ID: " <<  MChadron.GetPID() << std::endl;
    //        std::cout << " MC Total Momentum: " <<MChadron.GetMomentum().P() << std::endl;
    //        std::cout << " MC Polar Angles: (" << MChadron.GetMomentum().Theta() << ", "
    //                  << MChadron.GetMomentum().Phi() << ")" << std::endl;
    //        std::cout << " MCz value : " <<   MChadron.GetzMC()<< std::endl;
    //        std::cout << " MCpt2 value : " << MChadron.Getpt2MC()<< std::endl;
    //        std::cout << " MCphih value : " <<MChadron.GetphihMC()<< std::endl;
//
    //    //}
    //}
}

void Event::CalcPol(Particle& particle) {
    //making a generalized polar angle function so it can be appplied to any particle.
    //call this function in the calcAll function
    double theta = particle.GetMomentum().Theta();      //from root utilities 
    double phi = particle.GetMomentum().Phi();
    
    //store the angles in the Particle class
    particle.SetTheta(theta);
    particle.SetPhi(phi);
}

int Event::CalcKinematics(){
    double theta_e = acos(Constants::elBeam.Vect().Dot(electron.GetMomentum().Vect()) / (Constants::elBeam.Vect().Mag() * electron.GetMomentum().Vect().Mag()));
    thetaelectron = theta_e;
    acosyada=Constants::elBeam.Vect().Dot(electron.GetMomentum().Vect()) / (Constants::elBeam.Vect().Mag() * electron.GetMomentum().Vect().Mag());
    Q2 = 4 * Constants::elBeam.E() * electron.GetMomentum().E() * pow(sin(theta_e / 2), 2);
    nu  = Constants::elBeam.E() - electron.GetMomentum().E();
    y = nu / Constants::elBeam.E(); 
    w2 =  m_D* m_D + 2 * m_D * nu - Q2;
    xb = Q2 / (2 * m_D * nu);
    electron.SetKinVariables(Q2, nu, y, w2, xb);    //useless 
    return 0;

}
int Event::CalcMCKinematics(){
            //std::cout << "Beam Energy: " << Constants::elBeam.E() << std::endl;
        //std::cout << "Electron Energy: " << MCelectron.GetMomentum().E() << std::endl;
        //std::cout << "beam momentum: " << Constants::elBeam.Vect().Mag() << std::endl;
        //std::cout << "electron momentum: " << MCelectron.GetMomentum().Vect().Mag() << std::endl;
        //std::cout <<"electron Px: " << MCelectron.GetMomentum().Px() << std::endl;
        //std::cout <<"electron Py: " << MCelectron.GetMomentum().Py() << std::endl;
        //std::cout <<"electron Pz: " << MCelectron.GetMomentum().Pz() << std::endl;
    double theta_eMC = acos(Constants::elBeam.Vect().Dot(MCelectron.GetMomentum().Vect()) / (Constants::elBeam.Vect().Mag() * MCelectron.GetMomentum().Vect().Mag()));
    thetaelectronMC = theta_eMC;
    acosyadaMC = Constants::elBeam.Vect().Dot(MCelectron.GetMomentum().Vect()) / (Constants::elBeam.Vect().Mag() * MCelectron.GetMomentum().Vect().Mag());
    MCQ2 =  4 * Constants::elBeam.E() * MCelectron.GetMomentum().E() * pow(sin(theta_eMC / 2), 2);
/* 
   if (MCQ2 < 0.8){
        std::cout << " " << std::endl;
        std::cout << "MCQ2 below threshold!" << std::endl;
        std::cout << "Beam Energy: " << Constants::elBeam.E() << std::endl;
        std::cout << "Electron Energy: " << MCelectron.GetMomentum().E() << std::endl;
        std::cout << "beam momentum: " << Constants::elBeam.Vect().Mag() << std::endl;
        std::cout << "electron momentum: " << MCelectron.GetMomentum().Vect().Mag() << std::endl;
        std::cout <<"electron Px: " << MCelectron.GetMomentum().Px() << std::endl;
        std::cout <<"electron Py: " << MCelectron.GetMomentum().Py() << std::endl;
        std::cout <<"electron Pz: " << MCelectron.GetMomentum().Pz() << std::endl;
        std::cout << "theta_eMC: " << theta_eMC << std::endl;
        std::cout << "acosyadaMC: " << acosyadaMC << std::endl;
        std::cout << "MCQ2: " << MCQ2 << std::endl;
        
    }
*/
    MCnu  = Constants::elBeam.E() - MCelectron.GetMomentum().E();
    MCy = MCnu / Constants::elBeam.E(); 
    MCw2 =  m_D* m_D + 2 * m_D * MCnu - MCQ2;
    MCxb = MCQ2 / (2 * m_D * MCnu);
    MCelectron.SetKinVariables(MCQ2, MCnu, MCy, MCw2, MCxb);    //useless
    return 0;

}



void Event::calcAll(){
    CalcKinematics( );
    CalcPol(electron);
    for (Particle& hadron : hadrons) {
        // Add Cut PID
        hadron.CalcHadronKin(electron.GetMomentum(), hadron.GetMomentum());
        CalcPol(hadron);
        
    }
}

void Event::calcMCAll(){
    CalcMCKinematics( );
    for (Particle& MChadron : MChadrons) {
        MChadron.CalcMCHadronKin(MCelectron.GetMomentum(), MChadron.GetMomentum());
    }
}

double Event::GetQ2() const {
        return Q2;
}
double Event::GetQ2MC() const {
        return MCQ2;
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
double Event::GetnuMC() const {
    return MCnu   ;
}
double Event::GetyMC() const {
    return MCy;
} 
double Event::GetW2MC() const{
    return MCw2;
}
double Event::GetxbMC() const{
    return MCxb;
}
double Event::GetThetaElectron() const{
    return thetaelectron;
}

double Event::GetAcosyada() const{
    return acosyada;
}
double Event::GetThetaElectronMC() const{
    return thetaelectronMC;
}

double Event::GetAcosyadaMC() const{
    return acosyadaMC;
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


int Event::GetTargetType() const {
    return target_type; 
}
void Event::SetTargetType(int type) {
    target_type = type;
    //0=D, 1=Sn, 2=Cu   
}

int Event::GetHel() const {
    return hel; 
}
int Event::GetHelRaw() const {
    return helraw; 
}

int Event::GetEvtnbr() const {
    return evtnbr; 
}