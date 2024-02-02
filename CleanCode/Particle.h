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
    Particle(const TLorentzVector& , int, int );   
    const TLorentzVector& GetMomentum() const; 
    int GetEventIndex() const ;
    int GetPID() const ;  
    void SetMomentum( TLorentzVector);
    void SetTheta(double);
    void SetPhi(double);
    void SetQ2(double);
    void Setnu(double );
    void SetKinVariables(double, double, double, double, double);
    void SetVx(double); 
    void SetVy(double); 
    void SetParticleRow(int);
    void SetParticleIndex(int);

    double GetQ2() const ;
    double Getnu()  const ;
    double Gety() const ;
    double GetW2() const ;
    double Getxb() const ;
    double Getz()    const;
    double Getpt2()  const;
    double Getphih() const;
    double GetTheta() const;
    double GetPhi() const;
    double GetVx() const; 
    double GetVy() const; 
    int GetParticleRow() const;
    int GetParticleIndex() const;




    int CalcHadronKin(TLorentzVector,TLorentzVector);

private:
    TLorentzVector momentum;

    int eventIndex;
    int pid; 
    double z=0, pt2=0, phih=0;
    double Q2=0, nu=0, w2=0 ,y=0 ,xb=0;
    double vx, vy;
    double phi , theta ;
    int particleRow;
    int particleIndex;

};
#endif

