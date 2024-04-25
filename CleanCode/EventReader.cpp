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
            currentEvent.AddElectron(momentum, row, vertz);
            currentEvent.SetVertexZ(vertz);
            currentEvent.SetVertexX(vertx);
            currentEvent.SetVertexY(verty);
            
            //currentEvent.SetParticleRow(row);
            
        } else if (IsHadron(pid)) {
            currentEvent.AddHadron(momentum, pid,row, vertz);
            //currentEvent.SetVertexZ(vertz);
            //currentEvent.SetParticleRow(row);
            for (const Particle& hadron : currentEvent.GetHadrons()) {
                //proceed from here to set vertex info !!!
                //std::cout << "Hadron momentum: " << hadron.GetMomentum().Px() << std::endl;
                //hadron.SetVz(vertz);
            }
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


    void EventReader::AddHelInfo(int hel, int hel_raw){
        currentEvent.SetHel(hel);
        currentEvent.SetHelRaw(hel_raw);
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


    std::optional<Event> EventReader::ProcessEvent( hipo::event event, int eventNumber, bool isSimulated ) {
        currentEvent = Event();
        
        bool el_detect = false;
        double max_energy_electron = 0.0; 
        event.getStructure(RECgen);
        event.getStructure(RECcalo);
        event.getStructure(RECcher);
        event.getStructure(HELbank);
        event.getStructure(MCpart);
        if (MCpart.getRows()==0){isSimulated = true;  }
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
            if (pid == Constants::ELECTRON_PID) {

                double e_pcal = 0 ;
                double e_ecalin = 0 ;
                double e_ecalout = 0 ;
                double lu_pcal = 0 ;
                double lv_pcal = 0 ;
                double lw_pcal = 0 ;
                double sector_pcal = 0 ;
                double x_cal = 0;
                double y_cal = 0;
                double z_cal = 0;       //porbalby useless TBD
            //RECcalo.show();
            for (int c_row = 0; c_row < RECcalo.getRows(); ++c_row) {
                //int index = RECcalo.getInt("index", c_row);
                int pindex = RECcalo.getInt("pindex", c_row);
                //check if  i matches the pindex
                
                    //once pindex changes  cal values are reset
                if (pindex == i) {
                    int cal_layer = RECcalo.getInt("layer", c_row);
                    if (cal_layer == 1) {   //pcal
                        e_pcal = RECcalo.getFloat("energy", c_row);
                        lu_pcal = RECcalo.getFloat("lu", c_row);
                        lv_pcal = RECcalo.getFloat("lv", c_row);
                        lw_pcal = RECcalo.getFloat("lw", c_row);
                        sector_pcal = RECcalo.getInt("sector", c_row);
                        x_cal = RECcalo.getFloat("x", c_row);
                        y_cal = RECcalo.getFloat("y", c_row);
                        z_cal = RECcalo.getFloat("z", c_row); //porbalby useless TBD

                    }
                    if (cal_layer == 4) {   //ecal_in
                        if (RECcalo.getFloat("energy", c_row) > 0.01){
                            //filetring the 0 energy hits (non existent after first calo layer)
                            e_ecalin = RECcalo.getFloat("energy", c_row);
                        }
                    }
                    if (cal_layer == 7) {   //ecal_out
                        if (RECcalo.getFloat("energy", c_row) > 0.01){
                            e_ecalout = RECcalo.getFloat("energy", c_row);
                        }
                    } 
                }
            }
            //std::cout << "e_pcal: " << e_pcal << std::endl;
            //std::cout << "e_in: " << e_ecalin << std::endl;
            //std::cout << "e_out: " << e_ecalout << std::endl;
                                AddCaloInfo(pid, sector_pcal, lu_pcal, lv_pcal, lw_pcal, e_pcal, e_ecalin, e_ecalout);
                    AddCaloXYZ(pid, x_cal, y_cal, z_cal);
            }
            for (int cher_row = 0; cher_row < RECcher.getRows(); ++cher_row){
                int pindex_cher = RECcher.getInt("pindex", cher_row);
                
                double nphe15;
                double nphe16 ;
                if (pindex_cher == i){
                    int cher_detectornb = RECcher.getInt("detector", cher_row);
                    //AddCherInfo(pid, nphe, cher_detectornb);
                    if (cher_detectornb == 16){
                        nphe16 = RECcher.getFloat("nphe", cher_row);
                        //std::cout << "nphe16: " << nphe16 << std::endl;
                    }
                    if (cher_detectornb == 15){
                        nphe15 = RECcher.getFloat("nphe", cher_row);
                        //std::cout << "nphe15: " << nphe15 << std::endl;
                    }
                    AddCherInfo(pid, nphe15, nphe16);
                }


            }
            int hel_evt = 0;
            int hel_raw = 0;
            //for (int hel_row = 0; hel_row < HELbank.getRows(); ++hel_row){
            hel_evt = HELbank.getInt("helicity", 0);
            hel_raw = HELbank.getInt("helicityRaw", 0);
            AddHelInfo(hel_evt, hel_raw);
            if (isSimulated == true){
                for (int j = 0; j <MCpart.getRows(); ++i) {
                    int pid = MCpart.getInt("pid", j);
                    if (pid == Constants::POSITRON_PID) {
                        return std::nullopt;    
                    }
                    if ( pid == 0){    //mass == Constants::default_mass ||
                        continue;
                    }
                    double MCmass = GetMassID(pid);
                    double MCtargetvz = MCpart.getFloat("vz", j);
                    double MCtargetvx = MCpart.getFloat("vx", j);
                    double MCtargetvy = MCpart.getFloat("vy", j);

                    TLorentzVector MCmomentum;
                    MCmomentum.SetPx(MCpart.getFloat("px", j));
                    MCmomentum.SetPy(MCpart.getFloat("py", j));
                    MCmomentum.SetPz(MCpart.getFloat("pz", j));
                    TVector3 MCmomentum3D = MCmomentum.Vect();
                    MCmomentum.SetVectM(MCmomentum3D, mass);
                    if (pid == Constants::ELECTRON_PID) {
    	                double  electron_status = MCpart.getInt("status", j);
                        if ( electron_status < 0) {     //momentum.E() > max_energy_electron &&
                            max_energy_electron = MCmomentum.E();
                            ProcessParticle(MCmomentum , Constants::ELECTRON_PID,MCtargetvx,MCtargetvy,MCtargetvz, j );

                        }
                        el_detect = true;
                    } else if (IsHadron(pid) && el_detect ==true) {
                            ProcessParticle(momentum, pid,MCtargetvx,MCtargetvy,MCtargetvz,j ); 

                    }
                }
            }
                
            //}
            

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
    bool flag_RECgen = false;
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
                RECevt  = factory.getSchema("REC::Event");  //using bank to recover helicity and helicity raw
                HELbank = factory.getSchema("HEL::online"); //can be HEL::online or HEL::flip or HEL::raw? eventually 
                //HELbank = factory.getSchema("REC::Event"); //can be HEL::online or HEL::flip or HEL::raw? eventually 
                //HELbank = factory.getSchema("HEL::flip"); //can be HEL::online or HEL::flip or HEL::raw? eventually 
                MCpart = factory.getSchema("MC::Particle");

                std::cout << "Processing file: " << filename << std::endl;
            } 

            hipo::event hipoEvent;
            if (reader.next()) {
                reader.read(hipoEvent);

                //std::cout<<"starting event"<<std::endl;
                
                
                return ProcessEvent(hipoEvent, 0, false);    
                //prevoius line was commented and replace by the next one for CLAScollNOV
                //return ProcessEventsWithPositivePions(hipoEvent, 0);
            } 
            Event empty;
            return empty;
    }
