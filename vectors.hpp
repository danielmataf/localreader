#ifndef VECTORS_H
#define VECTORS_H
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TH1F.h"




class PhysicsVectors {
public:
    TLorentzVector *v_incil;   //incident lepton
	TLorentzVector *v_nucleontarg;   //incident NUCLEON
    TLorentzVector *v_dif4epi;
    TLorentzVector *v_dif4piX;
    TLorentzVector *v_dif4eX;
    TLorentzVector *v_dif4epiX;
    TLorentzVector *v_diffmass; // Reusing this vector for X (hadron)
    TLorentzVector *v_diffmiss;
    TLorentzVector *v_misspho;
    TLorentzVector *v_scal;    // Scattered lepton
    TLorentzVector *v_scaph;   // Scattered REAL photon
    TLorentzVector *v_vipho;   // Virtual photon
    TLorentzVector *v_diff;
    TLorentzVector *v_scapip;
    TLorentzVector *v_scapim;
    TLorentzVector *v_scapit;
    TLorentzVector *v_scaneu;
    TLorentzVector *v_scapro;

    PhysicsVectors() {
        v_incil = new TLorentzVector;
        v_nucleontarg = new TLorentzVector;
        v_dif4epi = new TLorentzVector;
        v_dif4piX = new TLorentzVector;
        v_dif4eX = new TLorentzVector;
        v_dif4epiX = new TLorentzVector;
        v_diffmass = new TLorentzVector;
        v_diffmiss = new TLorentzVector;
        v_misspho = new TLorentzVector;
        v_scal = new TLorentzVector;
        v_scaph = new TLorentzVector;
        v_vipho = new TLorentzVector;
        v_diff = new TLorentzVector;
        v_scapip = new TLorentzVector;
        v_scapim = new TLorentzVector;
        v_scapit = new TLorentzVector;
        v_scaneu = new TLorentzVector;
        v_scapro = new TLorentzVector;
    }

    ~PhysicsVectors() {
        delete v_dif4epi;
        delete v_dif4piX;
        delete v_dif4eX;
        delete v_dif4epiX;
        delete v_diffmass;
        delete v_diffmiss;
        delete v_misspho;
        delete v_scal;
        delete v_scaph;
        delete v_vipho;
        delete v_diff;
        delete v_scapip;
        delete v_scapim;
        delete v_scapit;
        delete v_scaneu;
        delete v_scapro;
    }
};


#endif