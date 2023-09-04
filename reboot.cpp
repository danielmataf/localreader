
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


	////////////////	VARIABLES	/////////////////

	const int ELECTRON_PID = 11;
	const int PION_PID = 211;
	const int POSITRON_PID = -11;
	const double m_pi = 0.134976, m_e = 0.000511;
	double  electron_status;



	////////////////	CLASSES 	/////////////////


class Particle {
public:
    Particle(const TLorentzVector& momentum, int eventIndex) : momentum(momentum), eventIndex(eventIndex) {}

    const TLorentzVector& GetMomentum() const {
        return momentum;
        //useful to get momentum vectors
    }

    int GetEventIndex() const {
        return eventIndex;
        //useful to get the event number of the cosidered particle 
    }

private:
    TLorentzVector momentum;
    int eventIndex;
};

class Event {
public:
    Event(int eventIndex) : eventIndex(eventIndex) {}

    void AddElectron(const TLorentzVector& electronMomentum) {
        electrons.push_back(Particle(electronMomentum, eventIndex));
    }
	//both functions used to add a pion or an elet=ctron to a given vector called electron/pionMomentum
    void AddPion(const TLorentzVector& pionMomentum) {
        pions.push_back(Particle(pionMomentum, eventIndex));
    }

    int GetEventIndex() const {
        return eventIndex;
    }

    const std::vector<Particle>& GetElectrons() const {
        return electrons;
    }

    const std::vector<Particle>& GetPions() const {
        return pions;
    }

private:
    int eventIndex;
    std::vector<Particle> electrons;
    std::vector<Particle> pions;
    // electrons and pions are lists (c++ vector) 
    // each element of the list contains a given pion/electron with its momentum and index    
};


	////////////////	FUNCTIONS   	/////////////////


std::vector<std::string> createfilelist(int typeoffile, int files){ 
	//function creating  a list of files to read 
	//type for D or Sn, and files for the nmbr of files, length of list (100)
	    std::string basePath;
	    std::string EndPath;
	    if (typeoffile==1)
	    {
		basePath = "../../files2read/r_eD-0";       //path of files in ifarm
		EndPath = ".hipo";

	    }
	    else
	    {
		basePath = "../../files2read/r_eSn-0";       //path of files in ifarm
		EndPath = ".hipo";
		
	    }
	    std::vector<std::string> filelist;
	    for (int filenb = 1; filenb <= files; filenb++) {
		std::string filePath = basePath + std::to_string(filenb) + EndPath;
		filelist.push_back(filePath);
	    }
	    return filelist; 
	}


void ProcessParticle(const TLorentzVector& momentum, int pid, Event* currentEvent) {
    Particle particle(momentum, currentEvent->GetEventIndex());
    
    if (pid == 11) {
        currentEvent->AddElectron(momentum);
    } else if (pid == 211) {
        currentEvent->AddPion(momentum);
    }
    // Other particle types are not considered.
}

void ProcessEvent(const hipo::event& event, Event* currentEvent) {
    //hipo::bank RECgen(event.getBank("REC::Particle"));
    //std::cout << "  event:test  "<< &event <<std::endl;

    hipo::reader reader;
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::bank RECgen(factory.getSchema("REC::Particle"));  //call REC Bank 
    int numRECrows = RECgen.getRows();
    bool el_detect = false;
    std::cout << "  rowss->  " << numRECrows<<std::endl;

    for (int i = 0; i < 1; ++i) {
        int pid = RECgen.getInt("pid", i);
        //std::cout << "  row:test  " << std::endl;
        //add here POSITRON  condition TBD
        if (pid == ELECTRON_PID) {
            TLorentzVector momentum;
            //TLorentzVector *momentum  = new TLorentzVector;
            momentum.SetPx(RECgen.getFloat("px", i));
            momentum.SetPy(RECgen.getFloat("py", i));
            momentum.SetPz(RECgen.getFloat("pz", i));
            std::cout << "  Momentum:test  " << RECgen.getFloat("pz", i) << std::endl;
            double mass = (pid == 11) ? m_e : m_pi;
            TVector3 momentum3D = momentum.Vect();
            momentum.SetVectM(momentum3D, mass);
            electron_status = RECgen.getFloat("status", i);
            //Add electron status condition here TBD
                //processPart inside the condition
            ProcessParticle(momentum, pid, currentEvent);
            
            el_detect = true;
        } else if (pid == PION_PID && el_detect) {
            TLorentzVector momentum;
            momentum.SetPx(RECgen.getFloat("px", i));
            momentum.SetPy(RECgen.getFloat("py", i));
            momentum.SetPz(RECgen.getFloat("pz", i));   
            double mass = m_pi;
            TVector3 momentum3D = momentum.Vect();
            momentum.SetVectM(momentum3D, mass);

            //ProcessParticle(momentum, pid, currentEvent);
        }
    }
}

void PrintEventInfo(const std::vector<Event>& events, int eventIndex) {
    const Event& event = events[eventIndex]; // Get the specific event

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

void EvtReader(const std::vector<std::string>& filelist) {
    hipo::reader reader;

    for (const auto& filename : filelist) {
        std::cout << "Processing file: " << filename << std::endl;

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
                ProcessEvent(event, &currentEvent);
                events.push_back(currentEvent);

            } else {
                std::cerr << "Error reading event!" << std::endl;
            }
        }
        PrintEventInfo(events, 11);
        //PrintEventInfo(events);
    }
}
///////////////////////////////
	////////////////	MAIN    	/////////////////


int main() {
    std::vector<std::string> filenames =	{"../../files2read/r_eD-01.hipo","../../files2read/r_eD-02.hipo","../../files2read/r_eSn-01.hipo","../../files2read/r_eSn-02.hipo"};
    //std::vector<std::string> filelist2 = createfilelist(1,1);
    EvtReader(filenames);
    int eventToDisplay = 11; // Change this to the event index you want to display
    
    return 0;
}
