#include <TLorentzVector.h>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "Particle.h"
#include "constants.h"



    Particle::Particle(const TLorentzVector& momentum, int pid, int row)
        : momentum(momentum), pid(pid), particleRow(row) {

        }
    const TLorentzVector& Particle::GetMomentum() const {
        return momentum;
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

    void Particle::SetParticleRow(int row) {
        particleRow = row;
    }

    int Particle::GetParticleRow() const {
        return particleRow  ;
    }

    double Particle::GetTheta() const {
        return theta;
    }

    double Particle::GetPhi() const {
        return phi;
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

    void Particle::SetTheta(double t) {
        theta = t;
    }

    void Particle::SetPhi(double p) {
        phi = p;
    }

    void Particle::SetMomentum( TLorentzVector momentum){

        this-> momentum = momentum;
    }
    //void Particle::SetVertexZ ( double vertz) {
    //    this->
    //}
    double Particle::GetVx() const {
        return vx;
    }
    double Particle::GetVy() const {
        return vy;
    }



    double Particle::GetPx() const {
    return px;
}

double Particle::GetPy() const {
    return py;
}

double Particle::GetPz() const {
    return pz;
}

void Particle::SetPx(double pxValue) {
    px = pxValue;
}

void Particle::SetPy(double pyValue) {
    py = pyValue;
}

void Particle::SetPz(double pzValue) {
    pz = pzValue;
}





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
        if (sign<0 ){
            phih =  acos(angle_norm) *180.0/Constants::PI +180;       
        }
        if (sign>0 ){   
            phih =  acos(angle_norm) *-180.0/Constants::PI+180; 
        }
        //consider BSA 4 phih TBD
        return 0;
    }
