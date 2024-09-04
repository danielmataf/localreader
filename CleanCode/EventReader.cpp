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
        simulatedEventCount =0;

    }

    bool EventReader::IsHadron(int pid) {
        return pid != Constants::ELECTRON_PID && pid != Constants::POSITRON_PID;
    }
    
    void EventReader::ProcessParticle(const TLorentzVector& momentum, int pid, double vertx,double verty,double vertz, int row, double chi2_row) {
        if (pid == Constants::ELECTRON_PID) {
            currentEvent.AddElectron(momentum, row, vertz, chi2_row);
            currentEvent.SetVertexZ(vertz);
            currentEvent.SetVertexX(vertx);
            currentEvent.SetVertexY(verty);
                //Setting VErtex on the EVENT!!!! not in the particle 
            //currentEvent.SetParticleRow(row);
            
        } else if (IsHadron(pid)) {
            currentEvent.AddHadron(momentum, pid,row, vertz, chi2_row);
            //currentEvent.SetVertexZ(vertz);
            //currentEvent.SetParticleRow(row);
            for (const Particle& hadron : currentEvent.GetHadrons()) {
                //proceed from here to set vertex info !!!
                //std::cout << "Hadron momentum: " << hadron.GetMomentum().Px() << std::endl;
                //hadron.SetVz(vertz);      // why is this commented ? ???
            }
        }
        
    }

    void EventReader::ProcessMCParticle(const TLorentzVector& MCmomentum, int MCpid, double MCvertx,double MCverty,double MCvertz, int MCrow ) {
        if (MCpid == Constants::ELECTRON_PID) {
            currentEvent.AddMCElectron(MCmomentum, MCrow, MCvertz);
            currentEvent.SetVertexZMC(MCvertz);       //!!!!!!
            //currentEvent.SetVertexXMC(MCvertx);
            //currentEvent.SetVertexYMC(MCverty);
                //maybe we need ro create SetMCvertex functions
            //currentEvent.SetParticleRow(row);
            
        } else if (IsHadron(MCpid)) {       //no need to change this
            currentEvent.AddMCHadron(MCmomentum, MCpid,MCrow, MCvertz);
            for (const Particle& MChadron : currentEvent.GetMCHadrons()) {
                //nothing here apparently
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
            currentEvent.electron.SetCalSector(sector);
            currentEvent.electron.Setlu(u);
            currentEvent.electron.Setlv(v);
            currentEvent.electron.Setlw(w);
            currentEvent.electron.SetEpcal(epcal);
            currentEvent.electron.SetEcalin(ecalin);
            currentEvent.electron.SetEcalout(ecalout);

            std::cout << "Added calo infoepcal: " << epcal << std::endl;


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
            currentEvent.electron.SetCalX(xcal);
            currentEvent.electron.SetCalY(ycal);
            currentEvent.electron.SetCalZ(zcal); //porbalby useless TBD
        }
        // else if (IsHadron(pid)) {
        //    SetCalX(xcal);
        //    SetCalY(ycal);
        //    SetCalZ(zcal); //porbalby useless TBD
        //}

    }

    void EventReader::AddCherInfo(int pid, double nphe15, double nphe16){
        if (pid == Constants::ELECTRON_PID) {
            currentEvent.electron.Setnphe15(nphe15);
            currentEvent.electron.Setnphe16(nphe16);
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


bool EventReader::isSimulatedData(hipo::event event) {
    bool isitsilulated = false;
    event.getStructure(MCpart);
    if (MCpart.getRows() != 0) {
        isitsilulated = true;
    }
    return isitsilulated;
}


    std::optional<Event> EventReader::ProcessEventMC( hipo::event event ) {
        currentEvent = Event();
        double max_energy_electronMC = 0.0;
        event.getStructure(MCpart);
        int counterel_MC=0;
        bool el_detect = false;
        // no need to worry about returning nullopt in this fct since everyevent in MC bank is supposed to be as it should
        for (int j = 0; j <MCpart.getRows(); ++j) {
            int MCpid = MCpart.getInt("pid", j);
            double MCmass = GetMassID(MCpid);
            double MCtargetvz = MCpart.getFloat("vz", j);
            //std::cout << "MCtargetvz: " << MCtargetvz << std::endl; 
            double MCtargetvx = MCpart.getFloat("vx", j);
            double MCtargetvy = MCpart.getFloat("vy", j);
            TLorentzVector MCmomentum;
            MCmomentum.SetPx(MCpart.getFloat("px", j));
            MCmomentum.SetPy(MCpart.getFloat("py", j));
            MCmomentum.SetPz(MCpart.getFloat("pz", j));
            TVector3 MCmomentum3D = MCmomentum.Vect();
            MCmomentum.SetVectM(MCmomentum3D, MCmass);
            if (MCpid == Constants::ELECTRON_PID && MCmomentum.E() > max_energy_electronMC){
    	        //double  electron_status = MCpart.getInt("status", j);
                //if ( electron_status < 0) {     //momentum.E() > max_energy_electron &&
                    counterel_MC++;
                    max_energy_electronMC = MCmomentum.E();
                    ProcessMCParticle(MCmomentum , Constants::ELECTRON_PID,MCtargetvx,MCtargetvy,MCtargetvz, j );
                    ReadRunconfig(event);
               //}
                el_detect = true;
            } else if (IsHadron(MCpid) && el_detect ==true) {
                    ProcessMCParticle(MCmomentum, MCpid,MCtargetvx,MCtargetvy,MCtargetvz,j ); 
            }
        }
        return currentEvent;
    }



    //std::optional<Event> EventReader::ReadRunconfig(hipo::event event){
    //    currentEvent = Event();
    //    event.getStructure(RUNconfig);
    //    int evtnbr = RUNconfig.getInt("event", 0);
    //    currentEvent.Setevtnbr(evtnbr);
    //    return currentEvent;
//
//
    //}
    void EventReader::ReadRunconfig(hipo::event event){
    event.getStructure(RUNconfig);
    int evtnbr = RUNconfig.getInt("event", 0);
    currentEvent.Setevtnbr(evtnbr);
}


    //std::optional<Event> EventReader::ReadRunconfigMC(hipo::event event);



    std::optional<Event> EventReader::ProcessEvent( hipo::event event, int eventNumber, bool isSimulated ) {
        currentEvent = Event();
        //std::cout << "Processing event: " << eventNumber << std::endl;
        bool el_detect = false;
        int electron_count = 0;
        double max_energy_electron = 0.0; 
        event.getStructure(RECgen);
        event.getStructure(RECcalo);
        event.getStructure(RECcher);
        event.getStructure(HELbank);
        event.getStructure(MCpart);

        
        //RECgen.show();
        bool flag_el = false;
        int counter_el= 0.0;
        
        for (int i = 0; i < RECgen.getRows(); ++i) {
            if ( RECgen.getInt("pid", i) == 11 ){
                flag_el = true;
                electron_count++;
            }
        }
        if (electron_count > 1) { return std::nullopt;} //DISCARDS EVENTS WITH MORE THAN ONE ELECTRON
        if (flag_el == false )return std::nullopt;
        //std::cout << "RECgen rows : " << RECgen.getRows() << std::endl;
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
            double chi2_row = RECgen.getFloat("chi2pid", i);
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
                    //considering trigger electron , not most energetic one
                    //but trigger el should be the most energetic one 
                    ProcessParticle(momentum , Constants::ELECTRON_PID,targetvx,targetvy,targetvz, i, chi2_row);
                    ReadRunconfig(event);
                }
                el_detect = true;
            } else if (IsHadron(pid) && el_detect ==true) {
                //if (pid == Constants::PION_PLUS_PID){
                    ProcessParticle(momentum, pid,targetvx,targetvy,targetvz,i, chi2_row ); 

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
                    //layers (1,4,7) for each calo part, explains why many pindex are the same

                    if (cal_layer == 1) {   //pcal
                        e_pcal = RECcalo.getFloat("energy", c_row);
                        std :: cout << "source e_pcal: " << e_pcal << std::endl;
                        //we interested in saving this value for a given particle. 
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
                    std::cout<< "Adding to calo e_pcal: " << e_pcal << std::endl;
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
        }
        
                
            //}
            
        
        return currentEvent;

    }
    int EventReader::getevtnbr(){
        std::string filename;
        filename = filelist.at(filenb);
        reader.open(filename.c_str());  
        int evtnbr = reader.getEntries();
        return evtnbr;

    }

    int EventReader::getSimulatedCount() const {
    return simulatedEventCount;
    }




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
            //RUNconfig = factory.getSchema("RUN::config");
            std::cout << "Processing file: " << filename << std::endl;
        } 
        hipo::event hipoEvent;
        if (reader.next()) {
            reader.read(hipoEvent);
            //std::cout<<"starting event"<<std::endl;
            Event event ;
            if (isSimulatedData(hipoEvent)==true) {
                ProcessEventMC(hipoEvent);
                //std::cout << "Simulated data" << std::endl;
            }
            return ProcessEvent(hipoEvent, 0, false);    
            //prevoius line was commented and replace by the next one for CLAScollNOV
            //return ProcessEventsWithPositivePions(hipoEvent, 0);
        } 
        Event empty;
        return empty;
}


std::optional<Event> EventReader::ProcessEventsInFileMC() { 
        std::string filename;
        //change reader has next to 
        if (!reader.hasNext()) {
            filenb++;
            filename = filelist.at(filenb);
            reader.open(filename.c_str());  
            reader.readDictionary(factory);
            RUNconfig = factory.getSchema("RUN::config");
            MCpart = factory.getSchema("MC::Particle");
            std::cout << "Processing file MC: " << filename << std::endl;   //I guess we can delete this!
        } 
        hipo::event hipoEvent;
        //if (reader.next()) {
            reader.read(hipoEvent);
            Event event ;
            if (isSimulatedData(hipoEvent)==true) {
                ProcessEventMC(hipoEvent);
                return ProcessEventMC(hipoEvent);    
            }
        //} 
        Event empty;
        return empty;
}