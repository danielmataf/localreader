#ifndef EVENTPROCESSOR_H
#define EVENTPROCESSOR_H

#include <TLorentzVector.h>
#include <vector>
#include "reader.h"  
#include "Event.h"



static const int ELECTRON_PID = 11;
static const int POSITRON_PID = -11;
static const double mass_pi = 0.134976, mass_e = 0.000511;  //CREATE constantes.h to include csts



class EventProcessor {
public:
    EventProcessor(){
        reader.readDictionary(factory);
        RUNconfig = factory.getSchema("RUN::config");
        RECgen = factory.getSchema("REC::Particle");

    }
    void ProcessParticle(const TLorentzVector& momentum, int pid, Event* currentEvent);
    void PrintEventInfo(int eventIndex);
    Event ProcessEvent(const hipo::event& event, int eventNumber);
    void EvtReader(const std::vector<std::string>& filelist, const double m_e, const double m_pi);
    void ProcessEventsInFile(const std::string& filename);
    void NextFile (const std::vector<std::string>& filelist);


private:
    //const double m_e;
    //const double m_pi;
    //Do i need to call these csts here ? Theyre already defined in the declaration part  TBD  
  //  HipoUtils hipoClass; //calling a previous class 4 databanks
    //std::vector<Event> events;      //DELETE THISS!!!!!
    //hipoclass should be called at the beginning of "public" part of the class so it can be used in every fct of the class 
    hipo::bank RECgen;  
    hipo::bank RUNconfig;  
    hipo::reader reader;    
    hipo::event event;    
    hipo::dictionary factory;

    bool IsHadron(int pid) {
        //KEEP, this fct is OK 
        return pid != ELECTRON_PID && pid != POSITRON_PID;
    }
};







//private:
//    // Member variables
//    std::vector<Event> events;
//    hipo::bank RECgen;
//    hipo::bank RUNconfig;
//    hipo::reader reader;
//    hipo::event event;
//    hipo::dictionary factory;
//
//    // Helper function declaration
//    bool IsHadron(int pid);
//};

#endif // EVENTPROCESSOR_H
