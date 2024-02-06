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
    void AddElectron(const TLorentzVector& , int); 
    void AddHadron(const TLorentzVector& , int, int ); 
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

     void CalcPol(Particle&);

    double GetVz()const ;
    double GetVx()const ;
    double GetVy()const ;

    const std::vector<Particle>& GetHadrons() const;
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
    double vz;  //vertex of the event (?) (should be the same every time)
    double vx;
    double vy;
    //CALO data (only available for electrons)
    double cal_sector;
    double lu;
    double lv;
    double lw;
    double Epcal=0, Ecalin=0, Ecalout = 0;
    double calox=0, caloy=0, caloz=0;


};
//
#endif





