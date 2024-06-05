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
    Particle(const TLorentzVector& , int, int , double);   
    const TLorentzVector& GetMomentum() const; 
    const TLorentzVector& GetMCMomentum() const;
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
    void SetVz(double);
    void SetParticleRow(int);
    void SetParticleIndex(int);
    void SetPx(double );
    void SetPy(double );
    void SetPz(double );

    double GetQ2() const ;
    double Getnu()  const ;
    double Gety() const ;
    double GetW2() const ;
    double Getxb() const ;
        double Getz()    const;
    double GetzMC()    const;
        double Getpt2()  const;
    double Getpt2MC()  const;
        double Getphih() const;
    double GetphihMC() const;
    double GetTheta() const;
    double GetPhi() const;
    double GetVx() const; 
    double GetVy() const;
    double GetVz() const; 
    double GetPx(TLorentzVector) const;
    double GetPy(TLorentzVector) const;
    double GetPz(TLorentzVector) const;
 
    int GetParticleRow() const;
    int GetParticleIndex() const;
    double GetParticleVertexZ() const;




    int CalcHadronKin(TLorentzVector,TLorentzVector);
    int CalcMCHadronKin(TLorentzVector,TLorentzVector);

private:
    TLorentzVector momentum;
    TLorentzVector MCmomentum;

    int eventIndex;
    int pid; 
    double z=0, pt2=0, phih=0;
    double MCz=0, MCpt2=0, MCphih=0;

    double Q2=0, nu=0, w2=0 ,y=0 ,xb=0;
    
    double vx, vy ;
    double vz=0;
    double particleVertexZ = 0;
    double phi , theta ;
    double px, py, pz;
    int particleRow;
    int particleIndex;
    //double particleVertexZ;

    struct CalorimeterInfo {
        //int layer;
        //int sector;
        //double lu;
        //double lv;
        //double lw;
        //double e_pcal = 0, e_ecalin= 0,  e_ecalout=0 ;
    };

    //Not sure if a struct for calorimeter data is needed.
    //int layer; 
    //LAYERS not useful to stor, is mainly useful to determine which layere we're in to recover different energies 
    int sector; //suppose that if a particle crosses all 3 layers, it will be in the same sector CHECK!!!!
    double lu;
    double lv;
    double lw;
    double e_pcal = 0, e_ecalin= 0,  e_ecalout=0 ;
    //differentiationg Edep per layer, this should be associated to the layer nb the particle crossed. 
    // assigning default 0 values. Allows to see wich layer has been crossed. 
};
#endif

