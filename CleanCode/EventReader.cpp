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
#include <optional>

    EventReader::EventReader(const std::vector<std::string>& Files){
        filelist = Files; 
        filenb = 0;

    }

    bool EventReader::IsHadron(int pid) {
        return pid != Constants::ELECTRON_PID && pid != Constants::POSITRON_PID;
    }
    
    void EventReader::ProcessParticle(const TLorentzVector& momentum, int pid, double vertx,double verty,double vertz ) {
        if (pid == Constants::ELECTRON_PID) {
            currentEvent.AddElectron(momentum);
            currentEvent.SetVertexZ(vertz);
            currentEvent.SetVertexX(vertx);
            currentEvent.SetVertexY(verty);
        } else if (IsHadron(pid)) {
            currentEvent.AddHadron(momentum, pid);
        }
        
    }
    double EventReader::GetMassID(int id) {
            //a function that retunrs the mass according to the particle PID in argument
    switch (id) {
        case Constants::PION_PLUS_PID:
            return Constants::MASS_PION_PLUS;
        case Constants::PION_MINUS_PID:
            return Constants::MASS_PION_MINUS;
        case Constants::PION_ZERO_PID:
            return Constants::MASS_PION_ZERO;
        // others can be added
        default:
            return Constants::default_mass;
            
    }
}


    std::optional<Event> EventReader::ProcessEvent( hipo::event event, int eventNumber) {
        currentEvent = Event();
         
        bool el_detect = false;
        double max_energy_electron = 0.0; 
        event.getStructure(RECgen);
        //RECgen.show();
        bool flag_el = false;
        int counter_el= 0.0;

        for (int i = 0; i < RECgen.getRows(); ++i) {
            
            if ( RECgen.getInt("pid", i) == 11 ){
                flag_el = true;
            }
        }
        if (flag_el == false )return std::nullopt;
        for (int i = 0; i < RECgen.getRows(); ++i) {
            int pid = RECgen.getInt("pid", i);
            if (pid == Constants::POSITRON_PID) {
                return std::nullopt;    
            }
            if ( pid == 0){    //mass == Constants::default_mass ||
                continue;
            }
            //double mass = (pid == Constants::ELECTRON_PID) ? Constants::MASS_ELECTRON : Constants::MASS_PION; //if/else on masses
                //condition ward can be improved to a non binary statement TBD
            double mass = GetMassID(pid);
            double targetvz = RECgen.getFloat("vz", i);
            double targetvx = RECgen.getFloat("vx", i);
            double targetvy = RECgen.getFloat("vy", i);

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
                    //coinsider trigger electron , not most energetic one 
                    ProcessParticle(momentum , Constants::ELECTRON_PID,targetvx,targetvy,targetvz);
                    
                }
                el_detect = true;
            } else if (IsHadron(pid) && el_detect ==true) {
                //if (pid == Constants::PION_PLUS_PID){
                    ProcessParticle(momentum, pid,targetvx,targetvy,targetvz ); 
                //}
            }
        }
        return currentEvent; 
    }
    int EventReader::getevtnbr(){
        std::string filename;
        filename = filelist.at(filenb);
        reader.open(filename.c_str());  
        int evtnbr = reader.getEntries();
        return evtnbr;

    }

////////////////////////CLAScollNOV//////////////////////////


std::optional<Event> EventReader::ProcessEventsWithPositivePions(hipo::event event, int eventNumber) {
    currentEvent = Event();
    bool el_detect = false;
    double max_energy_electron = 0.0; 
    event.getStructure(RECgen);
    //RECgen.show();
    bool flag_el = false;
    int counter_el = 0.0;

    for (int i = 0; i < RECgen.getRows(); ++i) {
        if (RECgen.getInt("pid", i) == 11) {
            flag_el = true;
        }
    }
    if (flag_el == false) return std::nullopt;

    for (int i = 0; i < RECgen.getRows(); ++i) {
        int pid = RECgen.getInt("pid", i);
        if (pid == Constants::POSITRON_PID) {
            return std::nullopt;    
        }
        if (pid == 0) {
            continue;
        }

        // Check for positive pions (pi+)
        if (pid != Constants::PION_PLUS_PID) {
            continue; // Skip this hadron
        }

        double mass = GetMassID(pid);
        double targetvz = RECgen.getFloat("vz", i);
        double targetvx = RECgen.getFloat("vx", i);
        double targetvy = RECgen.getFloat("vy", i);

        TLorentzVector momentum;
        momentum.SetPx(RECgen.getFloat("px", i));
        momentum.SetPy(RECgen.getFloat("py", i));
        momentum.SetPz(RECgen.getFloat("pz", i));
        TVector3 momentum3D = momentum.Vect();
        momentum.SetVectM(momentum3D, mass);

        if (pid == Constants::ELECTRON_PID) {
            double electron_status = RECgen.getInt("status", i);
            if (electron_status < 0) {
                max_energy_electron = momentum.E();
                ProcessParticle(momentum, Constants::ELECTRON_PID, targetvx, targetvy, targetvz);
            }
            el_detect = true;
        } else if (IsHadron(pid) && el_detect == true) {
            ProcessParticle(momentum, pid, targetvx, targetvy, targetvz); 
        }
    }
    return currentEvent; 
}


/////////////////////////////////////////////////////////


    std::optional<Event> EventReader::ProcessEventsInFile() { 
            std::string filename;

            if (!reader.hasNext()) {
                filenb++;
                filename = filelist.at(filenb);
                reader.open(filename.c_str());  
                reader.readDictionary(factory);
                RUNconfig = factory.getSchema("RUN::config");
                RECgen = factory.getSchema("REC::Particle");

                std::cout << "Processing file: " << filename << std::endl;
            } 

            hipo::event hipoEvent;
            if (reader.next()) {
                reader.read(hipoEvent);

                //std::cout<<"starting event"<<std::endl;
                
                
                return ProcessEvent(hipoEvent, 0);    
                //prevoius line was commented and replace by the next one for CLAScollNOV
                //return ProcessEventsWithPositivePions(hipoEvent, 0);
            } 
            Event empty;
            return empty;
    }
