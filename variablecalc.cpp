#include <TArrow.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TLorentzVector.h>
#include "reader.h"


#define _USE_MATH_DEFINES
#include <math.h>



//double cordtheta(TLorentzVector *vec3  = new TLorentzVector)
//{
    //function can work for any particle represented with a vector
    //double theta_el = vec3->Theta();
  //  return theta_el;
//}

class QHistogram {
public:
  // Constructor  initializes the histogram with variable binning
  QHistogram() : h_Q("h_Q", "Q Histogram", 100, 0, 10) {}

  // fill the histogram with  Q values after applying the cuts
  void fill(double Q_value) {
    h_Q.Fill(Q_value);
  }

  // Draw??
  void draw() {
    h_Q.Draw();
  }

private:
  // TH1F histogram with variable binning
  TH1F h_Q;
};

hipo::reader reader;
hipo::dictionary factory;
hipo::event event;

hipo::bank REC_Particle_bank(factory.getSchema("REC::Particle"));  //call REC Bank

TLorentzVector fillVector2(TLorentzVector &vec1,double mass_particle, double evt) {
    // Create a TLorentzVector object with the given components
    //TLorentzVector myVector(px, py, pz, E);
    //reader.open(filelist[filenbr]);
    hipo::reader reader;
hipo::dictionary factory;
hipo::event event;
    hipo::bank REC_Particle_bank(factory.getSchema("REC::Particle"));  //call REC Bank
    reader.readDictionary(factory);
    reader.read(event);
    event.getStructure(REC_Particle_bank);
    //
    double Px1 = REC_Particle_bank.getFloat("px", evt);
    double Py1 = REC_Particle_bank.getFloat("py", evt);
    double Pz1 = REC_Particle_bank.getFloat("pz", evt);
    double E1 = std::sqrt(Px1*Px1 + Py1*Py1 + Pz1*Pz1 + mass_particle*mass_particle);
    vec1.SetPxPyPzE(Px1, Py1, Pz1, E1);
    return vec1;
}
  

    

double cordtheta(TLorentzVector* vec3)
{
    double theta_el = vec3->Theta();
    return theta_el;
}
double cordphi(TLorentzVector* vec3)
{
    //function can work for any particle represented with a vector
    double phi_el = vec3->Phi();
    return phi_el;
}

double calcQ(TLorentzVector* vec1 = nullptr, TLorentzVector* vec2 = nullptr) {
    if (!vec1 || !vec2) {
        // Handle null pointers
        return 0.0;
    }

    TVector3 p_scal = vec2->Vect();
    TVector3 p_incil = vec1->Vect();

    double theta_e = acos(p_incil.Dot(p_scal) / (p_incil.Mag() * p_scal.Mag()));
    double Q2 = 4 * vec1->E() * vec2->E() * pow(sin(theta_e / 2), 2);

    return Q2;
}

double calcgamnu(TLorentzVector* vec1 , TLorentzVector* vec2 )
{
    double gamnu= vec1->E() - vec2->E();
    return gamnu;
}
double calcy(TLorentzVector* vec1 , double nu)
{
    double y = nu / vec1->E();
    return y;
}




//double calcP(TVector3 vec3D = (TLorentzVector *vec3  = new TLorentzVector)->Vect())
double calcP(TLorentzVector* vec3 )
{
    //function maybe can convert directly a 4vect to a 3vect???
    //will get the lagnitude of any 3vect
    TVector3 vec3D = vec3->Vect();
    double p_el = vec3D.Mag();  
    return p_el;
}



double calcx(double Q, double mass, double nu)
{
    double x_b = Q/ (2 * mass * nu);	//CAREFUL, changes with diffeent targets
    return x_b;
}


//double calcz(TLorentzVector* vec4 = nullptr, double nu) 
double calcz(TLorentzVector* vec4, double nu)
{
    double z = (vec4->E())/nu;			//!!!!!!!!!!!!!!!!!!!!!!etc//
    return z;
}

double calcW(double Q, double mass, double nu)
{
    double W = sqrt(mass*mass + 2*mass*nu - Q);
    return W;
}

double calcmissm(TLorentzVector* vec1 ,TLorentzVector* vec2,TLorentzVector* vec3,TLorentzVector* vec4)
{
    TLorentzVector *v_diff  = new TLorentzVector;
    *v_diff= *vec1 + *vec2 -*vec3 - *vec4;  //creation of a missing mass vector
    double missmass = v_diff->M2();         //caalc the squared missing mass of the vector 
    return missmass;
}
double calcM2(TLorentzVector* vec1 )
{
    double masssq = vec1->M2();         //caalc the squared missing mass of the vector 
    return masssq;
}
double calcE(TLorentzVector* vec1 )
{
    double nrgv = vec1->E();         //caalc the squared missing mass of the vector 
    return nrgv;
}
double calcmag(TLorentzVector* vec1 ,TLorentzVector* vec2,TLorentzVector* vec3,TLorentzVector* vec4)
{
    TLorentzVector *v_diff  = new TLorentzVector;
    *v_diff= *vec1 + *vec2 -*vec3 - *vec4;  //creation of a missing mass vector
    double missmass = v_diff->Mag();         //caalc the squared missing mass of the vector 
    return missmass;
}



double calcP_trX(TLorentzVector* vec1 ,TLorentzVector* vec2)
{
    TVector3 vec3D_a = vec1->Vect();
    TVector3 vec3D_b = vec2->Vect();
    TVector3 crossvec = vec3D_a.Cross(vec3D_b);
    double P_tr = (crossvec.Mag() )/ (vec3D_a.Mag() );
    return P_tr;
}



/*
double calcmissP(TLorentzVector* vec1,TLorentzVector* vec2,TLorentzVector* vec3,TLorentzVector* vec4)
{
    TLorentzVector *v_diff  = new TLorentzVector;
    *v_diff= *vec1 + *vec2 -*vec3 - *vec4;  //creation of a missing vector
    TVector3 p_diff= v_diff->Vect();        //conversion to 3D
    double missP = p_diff->Mag();         //calc the  missing momenta of the vector 
    return missP;
}

*/
//pending TO DO getE get Px,y,z


