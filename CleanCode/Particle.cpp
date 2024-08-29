 #include <TLorentzVector.h>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "Particle.h"
#include "constants.h"



    Particle::Particle(const TLorentzVector& momentum, int pid, int row, double particleVertexZ)
        : momentum(momentum), pid(pid), particleRow(row), particleVertexZ(particleVertexZ) {

        }
    const TLorentzVector& Particle::GetMomentum() const {
        return momentum;
    }
    //    double Particle::GetVz() const {
    //    return particleVertexZ;
    //}
        double Particle::GetParticleVertexZ() const {
        return particleVertexZ;
    }


    int Particle::GetEventIndex() const {
        return eventIndex;
    }

    int Particle::GetPID() const {
        return pid;
    }
    void Particle::SetQ2(double q2) {
        Q2 = q2;
    }

    void Particle::Setnu(double v) {
        nu = v;
    }
    void Particle::SetKinVariables(double a, double v, double c, double d, double e){
        Q2 = a;
        nu = v;
        y  = c;
        w2 = d;
        xb = e; 

    }
    void Particle::SetVx(double vertX) {
        vx = vertX;
    }
    void Particle::SetVy(double vertY) {
        vy = vertY;
    }
    void Particle::SetVz(double particleVertexZ) {
        this-> particleVertexZ = particleVertexZ;
        
    }

    void Particle::SetParticleRow(int row) {
        particleRow = row;
    }

    int Particle::GetParticleRow() const {
        return particleRow  ;
    }


   
    double Particle::GetQ2() const {
        return Q2;
    }
    double Particle::Getnu() const {
        return nu;
    }
    double Particle::Gety() const {
        return y;
    } 
    double Particle::GetW2() const{
        return w2;
    }
    double Particle::Getxb() const{
        return xb;
    }
    double Particle::Getz()     const {
        return z; 
    }
    double Particle::Getpt2()   const {
        return pt2; 
    }
    double Particle::Getphih()  const {
        return phih; 
    }

    double Particle::GetzMC()     const {
        return MCz; 
    }
    double Particle::Getpt2MC()   const {
        return MCpt2; 
    }
    double Particle::GetphihMC()  const {
        return MCphih; 
    }


    void Particle::SetTheta(double t) {
        theta = t;
    }

    void Particle::SetPhi(double p) {
        phi = p;
    }

    void Particle::SetMomentum( TLorentzVector momentum){

        this-> momentum = momentum;
    }     
    //void Particle::SetPolarAngle( double theta, double phi){
    //    this->theta = theta;
    //    this->phi = phi;
    //}
    //void Particle::SetVertexZ ( double vertz) {
    //    this->
    //}
    double Particle::GetVx() const {
        return vx;
    }
    double Particle::GetVy() const {
        return vy;
    }



    double Particle::GetPx(TLorentzVector momentum) const {
        return momentum.Px();
    }

    double Particle::GetPy(TLorentzVector momentum) const {
        return momentum.Py();
    }

    double Particle::GetPz(TLorentzVector momentum) const {
        return momentum.Pz();
    }

    //double Particle::GetPhi(TLorentzVecotr momentum) const {
    //    return momentum.Phi();
    //}
    //double Particle::GetTheta(TLorentzVector momentum) const {
    //    return momentum.Theta();
    //}




                                //vec1=scatterdelec (scal)    vec2=scatteredhadron (scapip)
    int Particle::CalcHadronKin(TLorentzVector vec1 ,TLorentzVector vec2){
        TLorentzVector vipho = Constants::elBeam - vec1; 
        TVector3 n1 = vipho.Vect().Cross(Constants::elBeam.Vect());
        TVector3 n2 = vipho.Vect().Cross(vec2.Vect());
        TVector3 n3 = n1.Cross(n2);
        double angle_norm = ( n1.Dot(n2) ) / ( n1.Mag()*n2.Mag() );
        double sign = n3.Dot(vipho.Vect()); 
        z= (vec2.E())/(       Constants::elBeam.E() - vec1.E() ); // 
        pt2 =  pow(  (n2.Mag() )/ (vipho.Vect().Mag() )  ,    2);
        //if (sign<0 ){
        //    phih =  acos(angle_norm) *180.0/Constants::PI +180;       
        //}
        //if (sign>0 ){   
        //    phih =  acos(angle_norm) *-180.0/Constants::PI+180; 
        //}
        //double acos_angle_norm_deg = acos(angle_norm) * 180.0 / Constants::PI;
        phih = (sign < 0) ? acos(angle_norm) * 180.0 / Constants::PI + 180 : acos(angle_norm) * -180.0 / Constants::PI + 180;
        //std::cout << " phih: " << phih ;
        //std::cout << " sign: " << sign   ;

        //std::cout << " anglenorm: " << angle_norm << std::endl;
        //std::cout << " acos anglenorm: " <<  acos(angle_norm) << std::endl;
        //std::cout <<  n1.Mag() <<"  "   <<n2.Mag() <<"  "<<  vec2.Mag() << std::endl;
        //consider BSA 4 phih TBD
        return 0;
    }

    int Particle::CalcMCHadronKin(TLorentzVector vec1 ,TLorentzVector vec2){
        TLorentzVector vipho = Constants::elBeam - vec1; 
        TVector3 n1 = vipho.Vect().Cross(Constants::elBeam.Vect());
        TVector3 n2 = vipho.Vect().Cross(vec2.Vect());
        TVector3 n3 = n1.Cross(n2);
        double angle_norm = ( n1.Dot(n2) ) / ( n1.Mag()*n2.Mag() );
        double sign = n3.Dot(vipho.Vect()); 
        MCz=  (vec2.E())/(       Constants::elBeam.E() - vec1.E() ); // 
        MCpt2 =   pow(  (n2.Mag() )/ (vipho.Vect().Mag() )  ,    2);
        MCphih =  (sign < 0) ? acos(angle_norm) * 180.0 / Constants::PI + 180 : acos(angle_norm) * -180.0 / Constants::PI + 180;
        //consider BSA 4 phih TBD
        return 0;
    }
