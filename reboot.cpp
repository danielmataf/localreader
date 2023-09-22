#include <TCanvas.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <any>
#include "reader.h"
#include "vectors.hpp"
#include <cstdlib>
#include "EventProcessor.h"
#include "Event.h"


	////////////////	VARIABLES	/////////////////

	//const int ELECTRON_PID = 11;
	//const int PION_PID = 211;
	//const int POSITRON_PID = -11;

	////////////////	CLASSES 	/////////////////


//class HipoUtils {
//public:
//    //USELESS   
//    HipoUtils() {
//        //to be corrected TBD 
//        reader.readDictionary(factory);
//        RECgen = factory.getSchema("REC::Particle");
//        RUNconfig = factory.getSchema("RUN::config");
//        reader.readDictionary(factory);
//    }
//
//    hipo::bank RECgen;  
//    hipo::bank RUNconfig;  
//    hipo::reader reader;    
//    hipo::event event;    
//
//private:
//    hipo::dictionary factory;
//};




/*
class EventProcessor {
public:
    EventProcessor() {}
    //need to create a fct "next event that verifies if there is a next event, else it will open the next file
        // see elements in the EvtReader function 
    //eventReadeer has been integrated to the class. Still needs to be modified 
    //Add PrintEventInfo on the class and modify it so it can be used fo a single evt in a given time

    // feed qnan(?) filelist in argumt



    //purge process particle, probaby tautologic (no need) TBD 
    void ProcessParticle(const TLorentzVector& momentum, int pid, Event* currentEvent) {
        //useless fct (?)   or put it into event class 
        //or create a fct AddParticle inside the evtClass
        //maybe add the bool ishadron fct 
        Particle particle(momentum, currentEvent->GetEventIndex(), pid);
        if (pid == ELECTRON_PID) {
            currentEvent->AddElectron(momentum);
        } else if (IsHadron(pid)) {
            currentEvent->AddHadron(momentum, pid);
        }
    }

     void PrintEventInfo(int eventIndex) {
        const Event& event = events[eventIndex]; // Assuming 'events' is a member variable

        std::cout << "Event Index: " << event.GetEventIndex() << std::endl;

        std::cout << "Electrons:" << std::endl;
        for (const Particle& electron : event.GetElectrons()) {
            std::cout << "  Momentum: (" << electron.GetMomentum().Px() << ", "
                      << electron.GetMomentum().Py() << ", " << electron.GetMomentum().Pz() << ")" << std::endl;
        }

        // Print pions' momentum
        std::cout << "Pions:" << std::endl;
        for (const Particle& pion : event.GetHadrons()) {
            std::cout << "  Momentum: (" << pion.GetMomentum().Px() << ", "
                      << pion.GetMomentum().Py() << ", " << pion.GetMomentum().Pz() << ")" << std::endl;
        }
    }


    //void AddParticle(const TLorentzVector& momentum, int pid, Event* currentEvent) {
    //    //this fct should replace ProcessParticle
    //    Particle particle(momentum, currentEvent->GetEventIndex(), pid);
    //    currentEvent->AddParticle(momentum, pid);
    //}


    Event ProcessEvent(const hipo::event& event, int eventNumber) {
    	//THIS FCT ENTERS AN EVT AND LOOPS IN ITS ROWS.
        //FCT RETURNS AN EVT

        //std::cout << "Processing event: " << eventNumber << std::endl;  //OK!!!
        Event currentEvent(eventNumber);   //!!!!!!!    NEEDS TO BE ADJUSTED
        reader.readDictionary(factory);
        bool el_detect = false;
        double max_energy_electron = 0.0; // reference energy counter
        hipo::bank RUNevt(factory.getSchema("RUN::config"));
        RUNconfig = factory.getSchema("RUN::config");
        RECgen = factory.getSchema("REC::Particle");
        std::cout << "event " << std::endl;
        //std::cout << "testrows: " << RECgen.getRows() << std::endl; //OK!!!  
        //create fct to call RECgen.getrows !!!!!
        for (int i = 0; i < RECgen.getRows(); ++i) {
            int pid = RECgen.getInt("pid", i);
            if (pid == POSITRON_PID) {
                continue;
                //next evt(?) if positron 
            }
            
            double mass = (pid == ELECTRON_PID) ? mass_e : mass_pi; //if/else on masses
                //condition ward can be improved to a non binary statement 
                //also no need bc use of process particle 
    
            //std::cout << "working on particle #" <<i  << std::endl; //OK!!!  
            TLorentzVector momentum;
            momentum.SetPx(RECgen.getFloat("px", i));
            momentum.SetPy(RECgen.getFloat("py", i));
            momentum.SetPz(RECgen.getFloat("pz", i));
            TVector3 momentum3D = momentum.Vect();
            momentum.SetVectM(momentum3D, mass);
            if (pid == ELECTRON_PID) {
                    

                electron_status = RECgen.getInt("status", i);
                if (momentum.E() > max_energy_electron && electron_status < 0) {
                    max_energy_electron = momentum.E();
                    ProcessParticle(momentum, pid, &currentEvent);  // Add electron to the currentEvent
                        //skip this just add particle (?)
                        //no need 4 the process fct 
                    el_detect = true;
                }

                el_detect = true;
            } else if (IsHadron(pid) && el_detect) {
                ProcessParticle(momentum, pid, &currentEvent); // Add hadron to the currentEvent
            }
        }
        
        return currentEvent; // return the evt object w/ particles (TO BE TESTED)
    }

    void EvtReader(const std::vector<std::string>& filelist, const double m_e, const double m_pi) {
        int eventNumber = 0;
        //THIS FUNCTION SHOULD BE PURGED
        //hipo::bank RUNconfig = hipoClass.RUNconfig;

        //hipo::reader reader = hipoClass.reader;

        std::vector<Event> events; // 4 evt storage
        for (const auto& filename : filelist) {
            //separate/handle loop evts & loop files 
            //hipo::reader reader;
            //std::cout << "Processing file: " << filename << std::endl;
            reader.open(filename.c_str());
            //hipo::dictionary factory;
            //reader.readDictionary(factory);
            hipo::event event;
            hipo::bank RUNevt(factory.getSchema("RUN::config"));
            RECgen = factory.getSchema("REC::Particle");

            //  ALL THESE HIPO THINGS SHOULD BE CALLED ONLY ONCE AT THE START OF THE CLASS


            while ( reader.hasNext()) {
                if (reader.next()) {
                    /////////////////////////////
                    //FOLLOWING LINE SHOULD BE MOVED TO PROCESSEVENT FUNCTION
                    reader.read(event);
                    event.getStructure(RECgen);
                    event.getStructure(RUNevt);
                    int numRUNRows  = RUNevt.getRows();
                    int evt_nbr = RUNevt.getInt("event", 0);
                            Event currentEvent = ProcessEvent(event, evt_nbr);
                            //CALLING THE ACTUAL FUNCTION 
                    events.push_back(currentEvent);
                        //ONLY PUSHBACK CURRENTEVT IF IT IS A PROCESSED EVT
                } else {
                    std::cerr << "Error reading event!" << std::endl;
                }
            }
        }

    }

    void ProcessEventsInFiles(const std::vector<std::string>& filelist) {
    EventProcessor eventProcessor;

    //for (const std::string& filename : filelist) {
        // open file
        //hipo::reader reader;
    

            if (!reader.hasNext()) {
                // check if there is a next event in the current file
                // then read the next event; not the current (?)

                reader.open(filename.c_str());  //opening file
                std::cout << "Processing file: " << filename << std::endl;
            } 

            hipo::event hipoEvent;
                    //calling it hipoEvent to diffferentiate read events from saved events
            if (reader.next(hipoEvent)) {
                return ProcessEvent(hipoEvent, actual_evt_nbr);
                        //inside that fct need set a new evt counter of considered evts with belonging particles 
            } else {
                std::cerr << "Error reading event from file: " << filename << std::endl;
            }
   // }
    std::cout << "All files processed." << std::endl;
    }


private:
    //const double m_e;
    //const double m_pi;
    //Do i need to call these csts here ? Theyre already defined in the declaration part  TBD  
  //  HipoUtils hipoClass; //calling a previous class 4 databanks
    std::vector<Event> events;
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
*/
	////////////////	FUNCTIONS   	/////////////////






std::vector<std::string> createFileList(int typeOfFile, int files) {
    //function creating  a list of files to read 
	//type for D or Sn, and files for the nmbr of files, length of list (100)
    //this function is USELESS, can be improved. but is not being used 
    std::string basePath;
    std::string endPath;
    
    if (typeOfFile == 1) {
        basePath = "../../files2read/r_eD-0"; // Path of files in ifarm for type 1
    } else {
        basePath = "../../files2read/r_eSn-0"; // Path of files in ifarm for type 2
    }

    std::vector<std::string> fileList;
    for (int fileNum = 1; fileNum <= files; fileNum++) {
        std::string filePath = basePath + std::to_string(fileNum) + ".hipo";
        fileList.push_back(filePath);
    }
    
    return fileList;
}


    

/*
void EvtReader(const std::vector<std::string>& filelist, const double m_e, const double m_pi) {
    //including masses as arguments may be an issue when dealing w/ other hadrons TBD
    //no need to have masses as  arguments (yuck)
    // !! this fct should be inside the class evtprocessor 
    EventProcessor eventProcessor(m_e, m_pi);
    //importing the class 4 event handling
    int eventNumber = 0;    //not exactly a counter
    HipoUtils hipoClass;
    hipo::bank RUNconfig = hipoClass.RUNconfig;
    hipo::reader reader = hipoClass.reader;
    std::vector<Event> events; // 4 evt storage
    for (const auto& filename : filelist) { //loop through files 
        std::cout << "Processing file: " << filename << std::endl;
        reader.open(filename.c_str());

        while (reader.hasNext()) {  //looping through events (?)
            //verifying we have a next event (? ) 
            if (reader.next()) {
                //entering the evt (?) 
                reader.read(hipoClass.event); // Reading the current event
                std::cout << "before evt structure..." << std::endl;
                //couts to test 4 seg fault 
                hipoClass.event.getStructure(RUNconfig);        //problem is here. Watch the hipo includes!!!
                std::cout << "after evt structure..." << std::endl;

                int evt_nbr = RUNconfig.getInt("event",0);
                //there's no row-loop in this function. see processevent fct
                Event currentEvent = eventProcessor.ProcessEvent(hipoClass.event, evt_nbr);
                //Event currentEvent(evt_nbr);
                events.push_back(currentEvent);
            }else {
                std::cerr << "Error reading event!" << std::endl;
            }
        }
    }

}
*/
	////////////////	MAIN    	/////////////////


int main() {
    std::vector<std::string> filenames =	{"../../files2read/r_eD-01.hipo","../../files2read/r_eD-02.hipo","../../files2read/r_eSn-01.hipo","../../files2read/r_eSn-02.hipo"};
    //std::vector<std::string> filelist2 = createfilelist(1,1);
    //EvtReader(filenames, mass_e, mass_pi);
    std::string filetest =	{"../../files2read/r_eD-01.hipo"};
    EventProcessor eventProcessor;
    //eventProcessor.EvtReader(filenames, mass_e, mass_pi);
    //std::cout << "Processing file: " << filetest << std::endl;
    eventProcessor.ProcessEventsInFile(filetest);
    //eventProcessor.NextFile(filenames);
    //int eventToDisplay = 11; // Change this to the event index you want to display
    
    return 0;
}


/*
void PrintEventInfo(const std::vector<Event>& events) {
    for (const Event& event : events) {
        std::cout << "Event Index: " << event.GetEventIndex() << std::endl;

        // Print electrons' momentum
        std::cout << "Electrons:" << std::endl;
        for (const Particle& electron : event.GetElectrons()) {
            std::cout << "  Momentum: (" << electron.GetMomentum().Px() << ", "
                      << electron.GetMomentum().Py() << ", " << electron.GetMomentum().Pz() << ")" << std::endl;
        }

        // Print pions' momentum
        std::cout << "Pions:" << std::endl;
        for (const Particle& pion : event.GetPions()) {
            std::cout << "  Momentum: (" << pion.GetMomentum().Px() << ", "
                      << pion.GetMomentum().Py() << ", " << pion.GetMomentum().Pz() << ")" << std::endl;
        }
    }
}
*/
/*
void ReadEvents(const double m_e, const double m_pi) {
    EventProcessor eventProcessor(m_e, m_pi);

    // Assuming you have code to open your data source (e.g., a file) and create a reader object.
    hipo::reader reader;
    reader.readDictionary(factory);

    int eventNumber = 0;

    while (reader.next() == true) {
        hipo::event event;
        reader.read(event);

        // Process the event using the EventProcessor
        Event currentEvent = eventProcessor.ProcessEvent(event, eventNumber);

        // Now you have the currentEvent with particles inside
        // You can access its information as needed.
        int eventIndex = currentEvent.GetEventIndex();
        const std::vector<Particle>& electrons = currentEvent.GetElectrons();
        const std::vector<Particle>& hadrons = currentEvent.GetHadrons();

        // Do further processing or analysis with the event data here.

        // Increment the event number
        eventNumber++;
    }

    // Close your data source and perform any necessary cleanup.
}
}
*/

/*      //"working" savestate
void EvtReader(const std::vector<std::string>& filelist, const double m_e, const double m_pi) {
    //including masses as arguments may be an issue when dealing w/ other hadrons TBD
    EventProcessor eventProcessor(m_e, m_pi);
    //importing the class 4 event handling
    hipo::reader reader;
    int eventNumber = 0;
    std::vector<Event> events; // 4 evt storage

    for (const auto& filename : filelist) { //loop through files 
        std::cout << "Processing file: " << filename << std::endl;
        HipoUtils hipoClass;
        hipo::bank RUNconfig = hipoClass.RUNconfig;

        reader.open(filename.c_str());
            
        hipo::dictionary factory;
        reader.readDictionary(factory);
        hipo::event event;
        hipo::bank RUNevt(factory.getSchema("RUN::config"));
        
        //int numRUNRows  = RUNevt.getRows();

        std::vector<Event> events;

        while (reader.hasNext()) {

            if (reader.next()) {
                reader.read(event);
                // process event data
                event.getStructure(RUNevt);
                int evt_nbr = RUNevt.getInt("event",0); //goes inside the event loop
                    //FAILING!  /!\ 
                //there's no row-loop in this function. see processevent fct
                Event currentEvent(evt_nbr);
                //evtprcssr.ProcessEvent(event, &currentEvent);
                //events.push_back(currentEvent);

            } else {
                std::cerr << "Error reading event!" << std::endl;
            }
        }
        PrintEventInfo(events, 11);
        //PrintEventInfo(events);
    }
}
*/
///////////////////////////////

