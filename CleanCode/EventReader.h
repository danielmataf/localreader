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
    void ProcessParticleMomentum(int, double , double , double );
    void AddCaloInfo(int,  int, double, double,double,double,double,double); 
    void AddCaloXYZ(int, double , double , double );
    void AddCherInfo(int, double, double);
    void AddHelInfo( int, int );
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
    hipo::bank RECcalo; //for three calorimeters (pcal ecal_in and ecal_out)
    hipo::bank RECcher; //for both cherenkov detectors (15 and 16)
    hipo::bank HELbank;
    hipo::bank RECevt;
    //if other banks shoul be added, add here then propagate to ProcessEventsInFile in the .cpp

    hipo::reader reader;
    hipo::event event;
    hipo::dictionary factory;
    int filenb;
    std::vector<std::string> filelist;
    Event currentEvent;


};
#endif 
