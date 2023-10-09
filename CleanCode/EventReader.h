#ifndef EVENTREADER_H
#define EVENTREADER_H

#include <TLorentzVector.h>
#include <vector>
#include "reader.h"  
#include "dictionary.h"  
#include "Event.h"

class EventReader {
public:
    EventReader(const std::vector<std::string>& );
    void ProcessParticle(const TLorentzVector& , int  );
    bool IsHadron(int );
    double GetMassID(int); 
    void PrintEventInfo(int eventIndex);
    Event ProcessEventsInFile();
private:
    Event ProcessEvent( hipo::event event, int eventNumber);
    hipo::bank RECgen;
    hipo::bank RUNconfig;
    hipo::reader reader;
    hipo::event event;
    hipo::dictionary factory;
    int filenb;
    std::vector<std::string> filelist;
    Event currentEvent;


};
#endif 