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
    Particle electron;
    Particle MCelectron;
    std::vector<Particle> hadrons; 
    std::vector<Particle> MChadrons; 

    void AddElectron(const TLorentzVector& , int,double, double); 
    void AddMCElectron(const TLorentzVector& , int,double); 
    void AddHadron(const TLorentzVector& , int, int,double, double  ); 
    void AddMCHadron(const TLorentzVector& , int, int,double );
    //fcts Add[particle] have an extra argument for the chi2pid
    int GetEventIndex() const ;
    void SetVertexZ(double );
    void SetVertexZMC(double );
    void SetVertexX(double );
    void SetVertexXMC(double );
    void SetVertexY(double );
    void SetVertexYMC(double );
    void SetHel(int);
    void SetHelRaw(int);
    
    void Setevtnbr(int);


     void CalcPol(Particle&);

    double GetVz()const ;
    double GetVzMC()const ;
    double GetVx()const ;
    double GetVxMC()const ;
    double GetVy()const ;
    double GetVyMC()const ;

    const std::vector<Particle>& GetHadrons() const;
    const std::vector<Particle>& GetMCHadrons() const;

    const Particle& GetElectron() const;
    const Particle& GetElectronMC() const ;

    int CalcKinematics( );
    int CalcMCKinematics( );

    void Print();
    void PrintMC();
    void calcAll();
    void calcMCAll();
    double GetQ2() const ;
    double GetQ2MC() const ;

    
    double Getnu()  const ;
        double GetnuMC()  const ;
    double Gety() const ;
        double GetyMC() const ;
    double GetW2() const ;
        double GetW2MC() const ;
    double Getxb() const ;
        double GetxbMC() const ;
    //Getters of CALO data (only available for electrons)
    //Getters for Helicity
    int GetHel( ) const;
    int GetHelRaw( ) const;

    double GetThetaElectron() const;
    double GetAcosyada() const;
    double GetThetaElectronMC() const;
    double GetAcosyadaMC() const;



    int GetEvtnbr() const;
    int GetEvtnbrMC() const;

    int GetTargetType() const; 
    void SetTargetType(int );

    //double Getz()    const;
    //double Getpt2()  const;
    //double Getphih() const;





private:
    int eventIndex;
    int target_type;
     

    double Q2=0,nu=0,w2=0,y=0,xb=0;
    double MCQ2=0,MCnu=0,MCw2=0,MCy=0,MCxb=0;

    double m_D, m_Sn, m_Cu; //GeV/c^2
    double vz=0;  //vertex of the event (?) (should be the same every time)
    double MCvz=0;
    double thetaelectron=0;
    double acosyada=0;
    double thetaelectronMC=0;
    double acosyadaMC=0;
    double vx;
    double vy;
    double px_el, py_el, pz_el;
    double theta_el, phi_el;
    //CALO data (only available for electrons)
    //double cal_sector = -1;
    //double lu;
    //double lv;
    //double lw;
    //double Epcal=0.0, Ecalin=0, Ecalout = 0;
    //double calox=0, caloy=0, caloz=0;
    //double nphe15=0, nphe16=10;
    int hel =0, helraw=0;
    int evtnbr=0;



};
//
#endif





