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
//#include "EventProcessor.h"




class Event {
public:
    Event ();
    //Event(int eventIndex) : eventIndex(eventIndex) {}

    void AddElectron(const TLorentzVector& ); 
    void AddHadron(const TLorentzVector& , int ); 
    int GetEventIndex() const ;
    const std::vector<Particle>& GetElectrons() const; 
    const std::vector<Particle>& GetHadrons() const;
    
    


private:
    int eventIndex;
    std::vector<Particle> electrons;
    std::vector<Particle> hadrons; // generalized to include different types of hadrons
    // in the generalized category shouldn't i distinguish between pions.
};
//
#endif
