#include <TLorentzVector.h>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include "reader.h"
#include "Event.h"
#include "Particle.h"

Event::Event(){
    
};

void Event::AddElectron(const TLorentzVector& electronMomentum) {
        electrons.push_back(Particle(electronMomentum,  11 ));
        //correct here the 11 PID electron TBD!!
}

void Event::AddHadron(const TLorentzVector& hadronMomentum, int pid) {
    hadrons.push_back(Particle(hadronMomentum,  pid));
}

int Event::GetEventIndex() const {
    return eventIndex;
}

const std::vector<Particle>& Event::GetElectrons() const {
    return electrons;
}

const std::vector<Particle>& Event::GetHadrons() const {
    return hadrons;
}



