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
#include "dictionary.h"  

#include <cstdlib>

#include "Event.h"

#include "EventReader.h"

    EventReader::EventReader(const std::vector<std::string>& Files){
        filelist = Files; 
        filenb = 0;



    }
    //need to create a fct "next event that verifies if there is a next event, else it will open the next file
        // see elements in the EvtReader function 
    //Add PrintEventInfo on the class and modify it so it can be used fo a single evt in a given time
    // feed qnan(?) filelist in argumt

        bool EventReader::IsHadron(int pid) {
        return pid != ELECTRON_PID && pid != POSITRON_PID;
        //taking every other particle that is not an e- nor an e+
        //how to consider only pions without touching the functions ? TBD
    }

    
    void EventReader::ProcessParticle(const TLorentzVector& momentum, int pid) {
        //useless fct (?)   or put it into event class 
        //or create a fct AddParticle inside the evtClass
        //maybe add the bool ishadron fct 
            //previous line is useless? no need of eventIndex
        if (pid == ELECTRON_PID) {
            currentEvent.AddElectron(momentum);
        } else if (IsHadron(pid)) {
            currentEvent.AddHadron(momentum, pid);
        }
    }


    Event EventReader::ProcessEvent( hipo::event event, int eventNumber) {
    	//THIS FCT ENTERS AN EVT AND LOOPS IN ITS ROWS.
        //FCT RETURNS AN EVT

        //std::cout << "Processing event: " << eventNumber << std::endl;  //OK!!!
        bool el_detect = false;
        double max_energy_electron = 0.0; // reference energy counter

        std::cout << "inside the event " << std::endl;
//all, move to constructeur 
        event.getStructure(RECgen);
        std::cout<<"nb of rows: " << RECgen.getRows()<<std::endl ;
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
                    ProcessParticle(momentum, pid);  // Add electron to the currentEvent
                        //skip this just add particle (?)
                        //no need 4 the process fct 
                    el_detect = true;
                }

                el_detect = true;
            } else if (IsHadron(pid) && el_detect) {
                ProcessParticle(momentum, pid ); // Add hadron to the currentEvent
            }
        }
        
        return currentEvent; // return the evt object w/ particles (TO BE TESTED)
    }

    

                                            //std::vector<std::string>
                                            //const std::string& filename
                                            //const std::vector<std::string>& filename
    Event EventReader::ProcessEventsInFile() { //argument currently is 1 file
        //This fonction should be called  in a loop of files!!! because. it reads an event (and the next one) 
            std::string filename;

            if (!reader.hasNext()) {
                // check if there is not aroot-config --cflags --libs next event in the current file
                    //then open a file; 
                    //NEED TO MOFIFY THE FILE LIST FOR IR TO BE LOOPED AUTOMATICALLY (???)
                // then read the next event; not the current (?)
                filenb++;
                filename = filelist.at(filenb);
                reader.open(filename.c_str());  //opening file
                reader.readDictionary(factory);
                RUNconfig = factory.getSchema("RUN::config");
                RECgen = factory.getSchema("REC::Particle");

                std::cout << "Processing file: " << filename << std::endl;
            } 

            hipo::event hipoEvent;
                    //calling it hipoEvent to diffferentiate read events from saved events
            if (reader.next()) {
                reader.read(hipoEvent);

                std::cout<<"starting event"<<std::endl;
                return ProcessEvent(hipoEvent, 0);
                        //need to  return a given evt ? 
                        //inside that fct need set a new evt counter of considered evts with belonging particles 
            } 
            Event empty;
            return empty;
            //std::cout << "All evts processed?" << std::endl;
    }
    


