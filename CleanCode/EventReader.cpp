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
#include "constants.h"

    EventReader::EventReader(const std::vector<std::string>& Files){
        filelist = Files; 
        filenb = 0;

    }

    bool EventReader::IsHadron(int pid) {
        return pid != Constants::ELECTRON_PID && pid != Constants::POSITRON_PID;
    }
    
    void EventReader::ProcessParticle(const TLorentzVector& momentum, int pid) {
        if (pid == Constants::ELECTRON_PID) {
            currentEvent.AddElectron(momentum);
        } else if (IsHadron(pid)) {
            currentEvent.AddHadron(momentum, pid);
        }
        
    }


    Event EventReader::ProcessEvent( hipo::event event, int eventNumber) {
        currentEvent = Event(); 
        bool el_detect = false;
        double max_energy_electron = 0.0; // reference energy counter
        TLorentzVector max_e_momentum;      //reference electron max nrg moment
        event.getStructure(RECgen);

        for (int i = 0; i < RECgen.getRows(); ++i) {
            int pid = RECgen.getInt("pid", i);
            if (pid == Constants::POSITRON_PID) {
                continue;   //skip evt if positron 
            }
            double mass = (pid == Constants::ELECTRON_PID) ? Constants::MASS_ELECTRON : Constants::MASS_PION; //if/else on masses
                //condition ward can be improved to a non binary statement 
                //also no need bc use of process particle 
    
            //std::cout << "working on particle #" <<i  << std::endl; //OK!!!  
            TLorentzVector momentum;
            momentum.SetPx(RECgen.getFloat("px", i));
            momentum.SetPy(RECgen.getFloat("py", i));
            momentum.SetPz(RECgen.getFloat("pz", i));
            TVector3 momentum3D = momentum.Vect();
            momentum.SetVectM(momentum3D, mass);
            if (pid == Constants::ELECTRON_PID) {
    	        double  electron_status = RECgen.getInt("status", i);
                if ( electron_status < 0) {     //momentum.E() > max_energy_electron &&
                    max_energy_electron = momentum.E();
                    max_e_momentum = momentum;
                    //coinsider trigger electron , not most energetic one 
                    ProcessParticle(momentum , Constants::ELECTRON_PID);
                }
                el_detect = true;
            } else if (IsHadron(pid) && el_detect ==true) {
                ProcessParticle(momentum, pid ); // Add hadron to the currentEvent
            }
        }
        return currentEvent; 
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
    





