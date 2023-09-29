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
    //Particle(const TLorentzVector& momentum,  int pid) : momentum(momentum),  pid(pid) {}
    Particle(const TLorentzVector& , int );   //no need of arg nale !!!
    const TLorentzVector& GetMomentum() const; 
    int GetEventIndex() const ;
    int GetPID() const ;    
private:
    TLorentzVector momentum;
    int eventIndex;
    int pid; // adding Particle ID
};
#endif

