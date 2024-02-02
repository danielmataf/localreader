#ifndef EVENTREADER_H
#define EVENTREADER_H

#include <TLorentzVector.h>
#include <vector>
#include "reader.h"  
#include "dictionary.h"  
#include "Event.h"
#include <optional>
class EventReader {
public:
    EventReader(const std::vector<std::string>& );
    void ProcessParticle(const TLorentzVector& , int ,double, double, double, int );
    bool IsHadron(int );
    double GetMassID(int); 
    void PrintEventInfo(int eventIndex);
    std::optional<Event> ProcessEventsInFile();
    int getevtnbr();

private:
    std::optional<Event> ProcessEvent( hipo::event event, int eventNumber);
    std::optional<Event> ProcessEventsWithPositivePions(hipo::event event, int eventNumber) ;


    //Banks
    hipo::bank RECgen;
    hipo::bank RUNconfig;
    hipo::bank RECcalo;
    //if other banks shoul be added, add here then propagate to ProcessEventsInFile in the .cpp

    hipo::reader reader;
    hipo::event event;
    hipo::dictionary factory;
    int filenb;
    std::vector<std::string> filelist;
    Event currentEvent;


};
#endif 