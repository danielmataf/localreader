#ifndef PARTICLE_H
#define PARTICLE_H

#include <TLorentzVector.h>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

class Particle {
public:
    Particle(const TLorentzVector& , int );   
    const TLorentzVector& GetMomentum() const; 
    int GetEventIndex() const ;
    int GetPID() const ;  
    void SetQ2(double);
    void Setnu(double );
    void SetKinVariables(double, double, double, double, double);
    double GetQ2() const ;
    double Getnu() const ;
    double Gety() const ;
    double GetW2() const ;
    double Getxb() const ;
    double Getz()    const;
    double Getpt2()  const;
    double Getphih() const;
    void SetMomentum( TLorentzVector);


    int CalcHadronKin(TLorentzVector,TLorentzVector);

private:
    TLorentzVector momentum;
    TLorentzVector IncidentLepton;

    int eventIndex;
    int pid; 
    double z=0, pt2=0, phih=0;
    double Q2bis=0, nubis=0, w2bis=0 ,ybis=0 ,xbbis=0;


};
#endif

