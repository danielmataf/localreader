#include <TLorentzVector.h>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#define PI 3.14159265
#include "Particle.h"



    Particle::Particle(const TLorentzVector& momentum,  int pid) 
        : momentum(momentum),  pid(pid) {
    IncidentLepton.SetPxPyPzE(0.0, 0.0, 11.0, 11.0);

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
        Q2bis = q2;
    }

    void Particle::Setnu(double v) {
        nubis = v;
    }
    void Particle::SetKinVariables(double a, double v, double c, double d, double e){
        Q2bis = a;
        nubis = v;
        ybis  = c;
        w2bis = d;
        xbbis = e; 

    }

    double Particle::GetQ2() const {
        return Q2bis;
    }
    double Particle::Getnu() const {
        return nubis;
    }
    double Particle::Gety() const {
        return ybis;
    } 
    double Particle::GetW2() const{
        return w2bis;
    }
    double Particle::Getxb() const{
        return xbbis;
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



                                //vec1=scatterdelec (scal)    vec2=scatteredhadron (scapip)
    int Particle::CalcHadronKin(TLorentzVector vec1 ,TLorentzVector vec2){
        //vec1 = scal ; vec2 = scapip or scah 
        TLorentzVector vipho = IncidentLepton - vec1; 
        TVector3 n1 = vipho.Vect().Cross(IncidentLepton.Vect());
        TVector3 n2 = vipho.Vect().Cross(vec2.Vect());
        TVector3 n3 = n1.Cross(n2);
        double angle_norm = ( n1.Dot(n2) ) / ( n1.Mag()*n2.Mag() );
        double sign = n3.Dot(vipho.Vect());
            //TVector3 p_incil = v_incil->Vect();  //3D
        z= vec1.Px();//(vec2.E())/Getnu();
        pt2 = vec1.Py();// pow(  (n2.Mag() )/ (vipho.Vect().Mag() )  ,    2);
        phih = vec1.Pz();// acos(angle_norm)*180.0/PI;       // check sign TBD!!!!   //consider BSA as well TBD
        return 0;
    }
