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
    
    void EventReader::ProcessParticle(const TLorentzVector& momentum, int pid, double vertx,double verty,double vertz, int row ) {
        if (pid == Constants::ELECTRON_PID) {
            currentEvent.AddElectron(momentum, row);
            currentEvent.SetVertexZ(vertz);
            currentEvent.SetVertexX(vertx);
            currentEvent.SetVertexY(verty);
            
            //currentEvent.SetParticleRow(row);
            
        } else if (IsHadron(pid)) {
            currentEvent.AddHadron(momentum, pid,row);
            currentEvent.SetVertexZ(vertz);
            //currentEvent.SetParticleRow(row);
        }
        
    }

    void EventReader::ProcessParticleMomentum(int pid, double pxValue, double pyValue, double pzValue) {
        if (pid == Constants::ELECTRON_PID) {
            //currentEvent.SetPx(pxValue);
            //currentEvent.SetPy(pyValue);
            //currentEvent.SetPz(pzValue);
        } else if (IsHadron(pid)) {
            //SetPx(pxValue);
            //SetPy(pyValue);
            //SetPz(pzValue);
        }
    }

    void EventReader::AddCaloInfo(int pid, int sector, double u, double v, double w, double epcal, double ecalin, double ecalout ){
        if (pid == Constants::ELECTRON_PID) {
            currentEvent.SetCalSector(sector);
            currentEvent.Setlu(u);
            currentEvent.Setlv(v);
            currentEvent.Setlw(w);
            currentEvent.SetEpcal(epcal);
            currentEvent.SetEcalin(ecalin);
            currentEvent.SetEcalout(ecalout);

        } else if (IsHadron(pid)) {
            //SetU(u);
            //SetV(v);
            //SetW(w);
            //SetEpcal(epcal);
            //SetEecalin(ecalin);
            //SetEecalout(ecalout);
        }
    }

    void EventReader::AddCaloXYZ(int pid , double xcal , double ycal , double zcal ){
        //testing a first function, only recovering xyz for first layer of calo and only for electrons

        if (pid == Constants::ELECTRON_PID) {
            currentEvent.SetCalX(xcal);
            currentEvent.SetCalY(ycal);
            currentEvent.SetCalZ(zcal); //porbalby useless TBD
        }
        // else if (IsHadron(pid)) {
        //    SetCalX(xcal);
        //    SetCalY(ycal);
        //    SetCalZ(zcal); //porbalby useless TBD
        //}

    }

    void EventReader::AddCherInfo(int pid, double nphe15, double nphe16){
        if (pid == Constants::ELECTRON_PID) {
            currentEvent.Setnphe15(nphe15);
            currentEvent.Setnphe16(nphe16);
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
        event.getStructure(RECcalo);
        event.getStructure(RECcher);
        //RECgen.show();
        bool flag_el = false;
        int counter_el= 0.0;
        
        for (int i = 0; i < RECgen.getRows(); ++i) {
            
            if ( RECgen.getInt("pid", i) == 11 ){
                flag_el = true;
            }
        }
        if (flag_el == false )return std::nullopt;

        for (int i = 0; i <RECgen.getRows(); ++i) {
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
            //std::cout<< "phi : "<< momentum.Phi() << std::endl;
            //std::cout<< "theta : "<< momentum.Theta() << std::endl;
            TVector3 momentum3D = momentum.Vect();
            momentum.SetVectM(momentum3D, mass);
            if (pid == Constants::ELECTRON_PID) {
    	        double  electron_status = RECgen.getInt("status", i);
                if ( electron_status < 0) {     //momentum.E() > max_energy_electron &&
                    max_energy_electron = momentum.E();
                    //coinsider trigger electron , not most energetic one 
                    ProcessParticle(momentum , Constants::ELECTRON_PID,targetvx,targetvy,targetvz, i );
                    
                }
                el_detect = true;
            } else if (IsHadron(pid) && el_detect ==true) {
                //if (pid == Constants::PION_PLUS_PID){
                    ProcessParticle(momentum, pid,targetvx,targetvy,targetvz,i ); 
                //}
            }
            for (int c_row = 0; c_row < RECcalo.getRows(); ++c_row) {
                //int index = RECcalo.getInt("index", c_row);
                int pindex = RECcalo.getInt("pindex", c_row);
                //check if  i matches the pindex
                double e_pcal ;
                double e_ecalin ;
                double e_ecalout ;
                double lu_pcal ;
                double lv_pcal ;
                double lw_pcal ;
                double sector_pcal ;
                double x_cal;
                double y_cal;
                double z_cal;       //porbalby useless TBD
                
                    //once pindex changes  cal values are reset
                if (pindex == i) {
                    int cal_layer = RECcalo.getInt("layer", pindex);
                        //e_pcal = RECcalo.getFloat("energy", pindex);
                        //lu_pcal = RECcalo.getFloat("lu", pindex);
                        //lv_pcal = RECcalo.getFloat("lv", pindex);
                        //lw_pcal = RECcalo.getFloat("lw", pindex);
                        //sector_pcal = RECcalo.getInt("sector", pindex);
                        //x_cal = RECcalo.getFloat("x", pindex);
                        //y_cal = RECcalo.getFloat("y", pindex);
                        //z_cal = RECcalo.getFloat("z", pindex); //porbalby useless TBD

                    if (cal_layer == 1) {   //pcal
                        e_pcal = RECcalo.getFloat("energy", pindex);
                        lu_pcal = RECcalo.getFloat("lu", pindex);
                        lv_pcal = RECcalo.getFloat("lv", pindex);
                        lw_pcal = RECcalo.getFloat("lw", pindex);
                        sector_pcal = RECcalo.getInt("sector", pindex);
                        x_cal = RECcalo.getFloat("x", pindex);
                        y_cal = RECcalo.getFloat("y", pindex);
                        z_cal = RECcalo.getFloat("z", pindex); //porbalby useless TBD

                    }
                    if (cal_layer == 4) {   //ecal_in
                         e_ecalin = RECcalo.getFloat("energy", pindex);
                    }
                    if (cal_layer == 7) {   //ecal_out
                         e_ecalout = RECcalo.getFloat("energy", pindex);
                    } 
                    AddCaloInfo(pid, sector_pcal, lu_pcal, lv_pcal, lw_pcal, e_pcal, e_ecalin, e_ecalout);
                    AddCaloXYZ(pid, x_cal, y_cal, z_cal);
                }
            }
            for (int cher_row = 0; cher_row < RECcher.getRows(); ++cher_row){
                int pindex_cher = RECcher.getInt("pindex", cher_row);
                
                double nphe15;
                double nphe16 ;
                if (pindex_cher == i){
                    int cher_detectornb = RECcher.getInt("detector", pindex_cher);
                    //AddCherInfo(pid, nphe, cher_detectornb);
                    if (cher_detectornb == 16){
                        nphe16 = RECcher.getFloat("nphe", pindex_cher);
                        //std::cout << "nphe16: " << nphe16 << std::endl;
                    }
                    if (cher_detectornb == 15){
                        nphe15 = RECcher.getFloat("nphe", pindex_cher);
                        //std::cout << "nphe15: " << nphe15 << std::endl;
                    }
                    AddCherInfo(pid, nphe15, nphe16);
                }


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
                ProcessParticle(momentum, Constants::ELECTRON_PID, targetvx, targetvy, targetvz, i);
            }
            el_detect = true;
        } else if (IsHadron(pid) && el_detect == true) {
            ProcessParticle(momentum, pid, targetvx, targetvy, targetvz, i); 
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
                RECcalo = factory.getSchema("REC::Calorimeter");
                RECcher = factory.getSchema("REC::Cherenkov");

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
