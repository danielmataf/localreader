#ifndef EVENT_H
#define EVENT_H


#include <TLorentzVector.h>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <any>
#include <cstdlib>
#include "Particle.h"




class Event {
public:
    Event ();
    void AddElectron(const TLorentzVector& , int,double); 
    void AddHadron(const TLorentzVector& , int, int,double ); 
    int GetEventIndex() const ;
    void SetVertexZ(double );
    void SetVertexX(double );
    void SetVertexY(double );
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
    void SetHel(int);
    void SetHelRaw(int);

     void CalcPol(Particle&);

    double GetVz()const ;
    double GetVx()const ;
    double GetVy()const ;

    const std::vector<Particle>& GetHadrons() const;
    const Particle& GetElectron() const;
    int CalcKinematics( );
    void Print();
    void calcAll();
    double GetQ2() const ;
    double Getnu()  const ;
    double Gety() const ;
    double GetW2() const ;
    double Getxb() const ;
    
    //Getters of CALO data (only available for electrons)
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
    //Getters for Helicity
    int GetHel( ) const;
    int GetHelRaw( ) const;






    int GetTargetType() const; 
    void SetTargetType(int );

    //double Getz()    const;
    //double Getpt2()  const;
    //double Getphih() const;





private:
    int eventIndex;
    int target_type;
     
    Particle electron;
    std::vector<Particle> hadrons; 
    double Q2=0,nu=0,w2=0,y=0,xb=0;
    double m_D, m_Sn, m_Cu; //GeV/c^2
    double vz=0;  //vertex of the event (?) (should be the same every time)
    double vx;
    double vy;
    double px_el, py_el, pz_el;
    double theta_el, phi_el;
    //CALO data (only available for electrons)
    double cal_sector = -1;
    double lu;
    double lv;
    double lw;
    double Epcal=0, Ecalin=0, Ecalout = 0;
    double calox=0, caloy=0, caloz=0;
    double nphe15=0, nphe16=10;
    int hel =0, helraw=0;



};
//
#endif





