#ifndef VARIABLECALC_H_INCLUDED
#define VARIABLECALC_H_INCLUDED

#include <vector>
#include "TLorentzVector.h"
#include "Math/Vector3D.h"
#include "reader.h"


using namespace std;

double cordtheta(TLorentzVector* vec3);

double cordphi(TLorentzVector* vec3);

double calcQ(TLorentzVector* vec1 , TLorentzVector* vec2 );

double calcgamnu(TLorentzVector* vec1 , TLorentzVector* vec2 );

double calcy(TLorentzVector* vec1 , double nu);

double calcx(double Q, double mass, double nu);

double calcz(TLorentzVector* vec4, double nu);

double calcW(double Q, double mass, double nu);

double calcmissm(TLorentzVector* vec1 ,TLorentzVector* vec2,TLorentzVector* vec3,TLorentzVector* vec4);

double calcmissP(TLorentzVector* vec1,TLorentzVector* vec2,TLorentzVector* vec3,TLorentzVector* vec4);

double calcM2(TLorentzVector* vec1 );

double calcE(TLorentzVector* vec1 );

double calcmag(TLorentzVector* vec1 ,TLorentzVector* vec2,TLorentzVector* vec3,TLorentzVector* vec4);

TLorentzVector fillVector2(TLorentzVector &vec1,double mass_particle, double evt) ;

//TLorentzVector fillVector2(TLorentzVector *vec1, double mass_particle, double evt);

//class QHistogram;

double calcP_trX(TLorentzVector* vec1 ,TLorentzVector* vec2);



#endif //VARIABLECALC_H
