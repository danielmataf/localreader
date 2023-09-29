#include <TLorentzVector.h>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "Particle.h"

    //->GetIndex is currently not being used but could b helpful (?)


    Particle::Particle(const TLorentzVector& momentum,  int pid) 
        : momentum(momentum),  pid(pid) {}

    const TLorentzVector& Particle::GetMomentum() const {
        return momentum;
    }
    int Particle::GetEventIndex() const {
        return eventIndex;
    }
    int Particle::GetPID() const {
        return pid;
    }
