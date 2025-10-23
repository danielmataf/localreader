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
    Particle(const TLorentzVector& , int, int , double, double);   
    const TLorentzVector& GetMomentum() const; 
    const TLorentzVector& GetMCMomentum() const;
    int GetEventIndex() const ;
    int GetPID() const ;  
    double Getchi2() const;
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
    void SetChi2(double );

    //Setters of CALO data (only available for electrons)
    void SetCalSector(int );
    void Setlu(double );
    void Setlv(double );
    void Setlw(double );
    void SetEpcal(double ); //can be grouped in one function for all enrgs and add sum energy
    void SetEcalin(double );    //can be grouped in one function for all enrgs and add sum energy
    void SetEcalout(double );   //can be grouped in one function for all enrgs and add sum energy TBD
    void SetCalX(double );  //can be grouped in one fct XYZ
    void SetCalY(double );  //can be grouped in one fct XYZ
    void SetCalZ(double );  //can be grouped in one fct XYZ TBD
    void Setnphe15(double );
    void Setnphe16(double );



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
    double GetCalSector( )const;
    double Getlu( )const;
    double Getlv( )const;
    double Getlw( )const;
    double GetEpcal( )const;
    double GetEcalin( )const;
    double GetEcalout( )const;

    double GetCalX( ) const;
    double GetCalY( ) const;
    double GetCalZ( ) const;
    //Getters of Cherenkov nphe data (only available for electrons)
    double Getnphe15( ) const;
    double Getnphe16( ) const;

    //double GetChi2() const;

    int GetParticleRow() const;
    int GetParticleIndex() const;
    double GetParticleVertexZ() const;



    int CalcHadronKin(TLorentzVector,TLorentzVector);
    int CalcMCHadronKin(TLorentzVector,TLorentzVector);
    //BETA stuff 
    void   SetBeta(double b) { beta_meas_ = b; }
    double GetBeta() const   { return beta_meas_; }

    // Optional helper: expected β under a mass hypothesis dismiss this 
    double BetaExpected(double mass) const {
    const double p = momentum.P();            // you already store TLorentzVector
    return p / std::sqrt(p*p + mass*mass);
}

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
    double chi2pid = 0;
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
    
    double cal_sector = -1;
    double lu;
    double lv;
    double lw;
    double Epcal=0.0, Ecalin=0, Ecalout = 0;
    double calox=0, caloy=0, caloz=0;
    double nphe15=0, nphe16=0;
    //BETA stuff
    double beta_meas_ = std::numeric_limits<double>::quiet_NaN();

    //differentiationg Edep per layer, this should be associated to the layer nb the particle crossed. 
    // assigning default 0 values. Allows to see wich layer has been crossed. 
};
#endif

