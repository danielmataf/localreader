 #ifndef EVENTREADER_H
#define EVENTREADER_H

#include <TLorentzVector.h>
#include <vector>
#include "reader.h"  
#include "dictionary.h"  
#include "Event.h"
#include <optional>
#include <unordered_map>
#include <string>

class EventReader {
public:
    EventReader(const std::vector<std::string>& );
    void ProcessParticle(const TLorentzVector& , int ,double, double, double, int, double );
    //upper fct has an extra argument for the chi2pid
    void ProcessMCParticle(const TLorentzVector& , int ,double, double, double, int );
    void ProcessParticleMomentum(int, double , double , double );
    void AddCaloInfo(int,  int, double, double,double,double,double,double); 
    void AddCaloXYZ(int, double , double , double );
    void AddCherInfo(int, double, double);
    void AddHelInfo( int, int );
    bool IsHadron(int );
    double GetMassID(int); 
    void PrintEventInfo(int eventIndex);
    std::optional<Event> ProcessEventsInFile();
    std::optional<Event> ProcessEventsInFileMC();
    std::optional<Event> ProcessEventsInFileREC();
    int getevtnbr();
    int getSimulatedCount() const ;
    bool isSimulatedData(hipo::event) ;
    void ReadRunconfig(hipo::event event);

    void SetDuplicateMode(bool on) { duplicateMode = on; }  //this in order to call the reading funvction twice, when using CuSn files taht need to be recalled 
    // wahteve to track opened files and events per file
    const std::unordered_map<std::string, size_t>& fileEventCounts() const { return perFileCounts_; }
    const std::string& currentFile() const { return currentFilePath_; }
    void printFileReadSummary(std::ostream& os = std::cout) const;

private:
    std::optional<Event> ProcessEvent( hipo::event event, int eventNumber, bool isSimulated);
    std::optional<Event> ProcessEventMC(hipo::event event) ;
    std::optional<Event> ProcessEventsWithPositivePions(hipo::event event, int eventNumber) ;
    //std::optional<Event> ReadRunconfigMC(hipo::event event);

    //following statements for WARD  in case of empty or corrupted .hipo files
     bool duplicateMode = false;     // when true: return each event twice (no advance on 2nd)
    bool haveCached = false;        // indicates a cached decoded Event is available
    Event cachedEvent;              // copy of last decoded Event
    bool openNextValidRECFile();
    bool openNextValidMCFile();

    int simulatedEventCount=0;

    //Banks
    hipo::bank RECgen;
    hipo::bank RUNconfig;
    hipo::bank RECcalo; //for three calorimeters (pcal ecal_in and ecal_out)
    hipo::bank RECcher; //for both cherenkov detectors (15 and 16)
    hipo::bank HELbank;
    hipo::bank RECevt;
    hipo::bank MCpart;
    hipo::bank MCevt;
    //if other banks shoul be added, add here then propagate to ProcessEventsInFile in the .cpp

    hipo::reader reader;
    hipo::event event;
    hipo::dictionary factory;
    int filenb ;
    std::vector<std::string> filelist;
    Event currentEvent;
    hipo::event cached_;      //last raw event read by MC
    bool have_cached_ = false;

    //more shit to track files 
    std::string currentFilePath_;
    std::unordered_map<std::string, size_t> perFileCounts_;

    // call when opening a file
    void noteOpen_(const std::string& path) {
        currentFilePath_ = path;
        // ensure key exists so “opened but yielded 0 events” still appears
        perFileCounts_.emplace(currentFilePath_, 0);
    }
    // call when you successfully yield one event to the caller
    void noteEvent_() {
        if (!currentFilePath_.empty()) ++perFileCounts_[currentFilePath_];
    }


};
#endif 