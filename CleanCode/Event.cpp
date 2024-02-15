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
    
Event::Event() : electron(TLorentzVector(0.0, 0.0, 0.0, 0.0), 0, 0) {
    m_D = (Constants::MASS_NEUTRON + Constants::MASS_PROTON)/2 ;

}

void Event::AddElectron(const TLorentzVector& electronMomentum, int row) {
        electron = Particle(electronMomentum, 11 , row);  
        electron.SetMomentum(electronMomentum);
        
}



void Event::AddHadron(const TLorentzVector& hadronMomentum, int pid, int row) {
    hadrons.push_back(Particle(hadronMomentum,  pid, row));
}

int Event::GetEventIndex() const {
    return eventIndex;
}

const std::vector<Particle>& Event::GetHadrons() const {
    return hadrons;
}

const Particle& Event::GetElectron() const {
    return electron;
}

void Event::SetVertexZ(double vertexz){
    vz = vertexz;
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

double Event::GetVx()const {
    return vx;
}double Event::GetVy()const {
    return vy;
}

void Event::SetCalSector(int sector){
    cal_sector = sector;
}
void Event::Setlu(double u){
    lu = u;
}
void Event::Setlv(double v){
    lv = v;
}
void Event::Setlw(double w){
    lw = w;
}
void Event::SetCalX(double calx){
    calox = calx;
}
void Event::SetCalY(double caly){
    caloy = caly;
}
void Event::SetCalZ(double calz){
    caloz = calz;
}
void Event::SetEpcal(double epcal){
    Epcal = epcal;
}
void Event::SetEcalin(double ecalin){
    Ecalin = ecalin;
}
void Event::SetEcalout(double ecalout){
    Ecalout = ecalout;
}

void Event::Setnphe15(double nphe_15){
    nphe15 = nphe_15;
}
void Event::Setnphe16(double nphe_16){
    nphe16 = nphe_16;
}


void Event::Print() {   //add int v=0 as argument 4 different types of verbose TBD

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
    std::cout << "Hadrons:" << std::endl;
    for (const Particle& hadron : GetHadrons()) {
        //if (hadron.GetPID()== 211){
            std::cout << "  ____________ " <<  std::endl;
            std::cout << "  Particle ID: " << hadron.GetPID() << std::endl;
            std::cout << " Total Momentum: " << hadron.GetMomentum().P() << std::endl;
            std::cout << " Polar Angles: (" << hadron.GetMomentum().Theta() << ", "
                      << hadron.GetMomentum().Phi() << ")" << std::endl;
            std::cout << " z value : " << hadron.Getz()<< std::endl;
            std::cout << " pt2 value : " << hadron.Getpt2()<< std::endl;
            std::cout << " phih value : " << hadron.Getphih()<< std::endl;

        //}
    }
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
    Q2 = 4 * Constants::elBeam.E() * electron.GetMomentum().E() * pow(sin(theta_e / 2), 2);
    nu  = Constants::elBeam.E() - electron.GetMomentum().E();
    y = nu / Constants::elBeam.E(); 
    w2 =  m_D* m_D + 2 * m_D * nu - Q2;
    //W = sqrt(Mn*Mn + 2*Mn*gamnu - Q2);
    xb = Q2 / (2 * m_D * nu);
    electron.SetKinVariables(Q2, nu, y, w2, xb);    //useless 
    return 0;

}

//void SetKinVariables(double , double  , double , double, double);


void Event::calcAll(){
    CalcKinematics( );
    CalcPol(electron);
    for (Particle& hadron : hadrons) {
        hadron.CalcHadronKin(electron.GetMomentum(), hadron.GetMomentum());
        CalcPol(hadron);
        
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
double Event::Getlu( )const{
    return lu;
}
double Event::Getlv( )const{
    return lv;
}
double Event::Getlw( )const{
    return lw;
}
double Event::GetEpcal( )const{
    return Epcal;
}
double Event::GetEcalin( )const{
    return Ecalin;
}
double Event::GetEcalout( )const{
    return Ecalout;
}
double Event::GetCalX() const{
    return calox;
}
double Event::GetCalY() const{
    return caloy;
}
double Event::GetCalZ() const{
    return caloz;
}
double Event::GetCalSector( )const{
    return cal_sector;
}
double Event::Getnphe15( ) const{
    return nphe15;
}
double Event::Getnphe16( ) const{
    return nphe16;
}


int Event::GetTargetType() const {
    return target_type; 
}
void Event::SetTargetType(int type) {
    target_type = type;
    //0=D, 1=Sn, 2=Cu   
}
