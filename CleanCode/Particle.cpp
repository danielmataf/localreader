#include <TLorentzVector.h>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "Particle.h"
#include "constants.h"



    Particle::Particle(const TLorentzVector& momentum,  int pid) 
        : momentum(momentum),  pid(pid) {

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

    //double Particle::GetQ2() const {
    //    return Q2;
    //}
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
    void Particle::SetMomentum( TLorentzVector momentum){

        this-> momentum = momentum;
    }
    //void Particle::SetVertexZ ( double vertz) {
    //    this->
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
        if (sign<0 ){
            phih =  acos(angle_norm) *180.0/Constants::PI +180;       
        }
        if (sign>0 ){   
            phih =  acos(angle_norm) *-180.0/Constants::PI+180; 
        }
        //consider BSA 4 phih TBD
        return 0;
    }
