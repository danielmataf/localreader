#ifndef EVENT_H
#define EVENT_H


#include <TLorentzVector.h>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <any>
#include "reader.h"
#include <cstdlib>
//#include "EventProcessor.h"



class Particle {
public:
    Particle(const TLorentzVector& momentum, int eventIndex, int pid) 
        : momentum(momentum), eventIndex(eventIndex), pid(pid) {}

    const TLorentzVector& GetMomentum() const {
        return momentum;
    }
    int GetEventIndex() const {
        return eventIndex;
    }
    int GetPID() const {
        return pid;
    }
private:
    TLorentzVector momentum;
    int eventIndex;
    int pid; // adding Particle ID
};

class Event {
public:
    Event(int eventIndex) : eventIndex(eventIndex) {}

    void AddElectron(const TLorentzVector& electronMomentum) {
        electrons.push_back(Particle(electronMomentum, eventIndex, 11 ));
    }
    void AddHadron(const TLorentzVector& hadronMomentum, int pid) {
        hadrons.push_back(Particle(hadronMomentum, eventIndex, pid));
    }
    int GetEventIndex() const {
        return eventIndex;
    }

    const std::vector<Particle>& GetElectrons() const {
        return electrons;
    }
    const std::vector<Particle>& GetHadrons() const {
        return hadrons;
    }
private:
    int eventIndex;
    std::vector<Particle> electrons;
    std::vector<Particle> hadrons; // generalized to include different types of hadrons
    // in the generalized category shouldn't i distinguish between pions.
};
//
#endif