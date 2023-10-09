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
    //const std::vector<Particle>& GetElectrons() const; 
    const std::vector<Particle>& GetHadrons() const;
    //double cordtheta(TLorentzVector* );
    //double cordphi(TLorentzVector* );


    int CalcKinematics( );
    void Print();
    void calcAll();



private:
    int eventIndex;
    Particle electron;
    std::vector<Particle> hadrons; 
    double Q2=0,nu=0,w2=0,y=0,xb=0;
    double m_e =  0.000511, m_n = 0.939565, m_p = 0.938272 ; //GeV/c^2
    double m_D, m_Sn, m_Cu; //GeV/c^2

};
//
#endif





/* 
TO DO
class cut ou struct 
creer un strucs avec un setcut 
avec cuts = a2 low Q2 high (doubles )

class cuts ...
....


struct cutval {
    double Q2 etc //initializing 
    ...
}
CutValues CutVMC(1,5....)
Cut 



class cutSet
..
{
    double Q = 




    bool PassCut(event);
}

CutSet CutMC();
if (CutMC.PassCuts(test))
        do whatever 


cut CutBase, CutTest;
CutTest.set (value )






 */