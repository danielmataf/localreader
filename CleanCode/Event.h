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
    void AddElectron(const TLorentzVector& ); 
    void AddHadron(const TLorentzVector& , int ); 
    int GetEventIndex() const ;
    void SetVertexZ(double );
    double GetVz()const ;
    const std::vector<Particle>& GetHadrons() const;
    int CalcKinematics( );
    void Print();
    void calcAll();
    double GetQ2() const ;
    double Getnu()  const ;
    double Gety() const ;
    double GetW2() const ;
    double Getxb() const ;
    //double Getz()    const;
    //double Getpt2()  const;
    //double Getphih() const;





private:
    int eventIndex;
    Particle electron;
    std::vector<Particle> hadrons; 
    double Q2=0,nu=0,w2=0,y=0,xb=0;
    double m_D, m_Sn, m_Cu; //GeV/c^2
    double vz;  //vertex of the event (?) (should be the same every time)

};
//
#endif





