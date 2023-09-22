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
#include <cstdlib>

#include "Event.h"
#include "EventProcessor.h"


    //EventProcessor(     ) {}
    //need to create a fct "next event that verifies if there is a next event, else it will open the next file
        // see elements in the EvtReader function 
    //eventReadeer has been integrated to the class. Still needs to be modified 
    //Add PrintEventInfo on the class and modify it so it can be used fo a single evt in a given time

    // feed qnan(?) filelist in argumt
    
    void EventProcessor::ProcessParticle(const TLorentzVector& momentum, int pid, Event* currentEvent) {
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

//    void EventProcessor::PrintEventInfo(int eventIndex) {
//        const Event& event = events[eventIndex]; // Assuming 'events' is a member variable
//
//        std::cout << "Event Index: " << event.GetEventIndex() << std::endl;
//
//        std::cout << "Electrons:" << std::endl;
//        for (const Particle& electron : event.GetElectrons()) {
//            std::cout << "  Momentum: (" << electron.GetMomentum().Px() << ", "
//                      << electron.GetMomentum().Py() << ", " << electron.GetMomentum().Pz() << ")" << std::endl;
//        }
//
//        // Print pions' momentum
//        std::cout << "Pions:" << std::endl;
//        for (const Particle& pion : event.GetHadrons()) {
//            std::cout << "  Momentum: (" << pion.GetMomentum().Px() << ", "
//                      << pion.GetMomentum().Py() << ", " << pion.GetMomentum().Pz() << ")" << std::endl;
//        }
//    }


    //void AddParticle(const TLorentzVector& momentum, int pid, Event* currentEvent) {
    //    //this fct should replace ProcessParticle
    //    Particle particle(momentum, currentEvent->GetEventIndex(), pid);
    //    currentEvent->AddParticle(momentum, pid);
    //}

    Event EventProcessor::ProcessEvent(const hipo::event& event, int eventNumber) {
    	//THIS FCT ENTERS AN EVT AND LOOPS IN ITS ROWS.
        //FCT RETURNS AN EVT

        //std::cout << "Processing event: " << eventNumber << std::endl;  //OK!!!
        Event currentEvent(eventNumber);   //!!!!!!!    NEEDS TO BE ADJUSTED
        bool el_detect = false;
        double max_energy_electron = 0.0; // reference energy counter
        hipo::bank RUNevt(factory.getSchema("RUN::config"));
//all, move to constructeur 

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
                    

    	        double  electron_status = RECgen.getInt("status", i);
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

    void EventProcessor::EvtReader(const std::vector<std::string>& filelist, const double m_e, const double m_pi) {
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

                                            //std::vector<std::string>
                                            //const std::string& filename
                                            //const std::vector<std::string>& filename
    void EventProcessor::ProcessEventsInFile(const std::string& filename) { //argument currently is 1 file
        EventProcessor eventProcessor;
        //This fonction should be called  in a loop of files!!! because. it reads an event (and the next one) 


            if (!reader.hasNext()) {
                // check if there is not a next event in the current file
                    //then open a file; 
                    //NEED TO MOFIFY THE FILE LIST FOR IR TO BE LOOPED AUTOMATICALLY (???)
                // then read the next event; not the current (?)

                reader.open(filename.c_str());  //opening file
                std::cout << "Processing file: " << filename << std::endl;
            } 

            hipo::event hipoEvent;
                    //calling it hipoEvent to diffferentiate read events from saved events
            if (reader.next(hipoEvent)) {
                return ProcessEvent(hipoEvent, actual_evt_nbr);
                        //need to  return a given evt ? 
                        //inside that fct need set a new evt counter of considered evts with belonging particles 
            } else {
                //std::cerr << "Error reading event from file: " << filename << std::endl;
            }
   
            //std::cout << "All evts processed?" << std::endl;
    }
    
    void EventProcessor::NextFile (const std::vector<std::string>& filelist){
    //loop files in filelist, the call the ProcessEventsInFiles function 
        for (const std::string& filename : filelist) {
            ProcessEventsInFile(filename);
        }
    }


