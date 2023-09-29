#include <iostream>
#include <vector>
#include "reader.h"
#include "Event.h"
#include "EventReader.h"
//#include "EventReader.cpp"


//compile :
//g++ -o analyse main.cpp `root-config --cflags --libs` -I/home/matamoros/Desktop/localreader/hipo4/



void PrintEventInfo(const Event& event) {
    std::cout << "Event Index: " << event.GetEventIndex() << std::endl;

    // Print electrons' momentum
    std::cout << "Electrons:" << std::endl;
    for (const Particle& electron : event.GetElectrons()) {
        std::cout << "  Particle ID: " << electron.GetPID() << std::endl;
        std::cout << "  Momentum: (" << electron.GetMomentum().Px() << ", "
                  << electron.GetMomentum().Py() << ", " << electron.GetMomentum().Pz() << ")" << std::endl;
    }

    // Print hadrons' momentum
    std::cout << "Hadrons:" << std::endl;
    for (const Particle& hadron : event.GetHadrons()) {
        std::cout << "  Particle ID: " << hadron.GetPID() << std::endl;
        std::cout << "  Momentum: (" << hadron.GetMomentum().Px() << ", "
                  << hadron.GetMomentum().Py() << ", " << hadron.GetMomentum().Pz() << ")" << std::endl;
    }
}

int main() {
    std::vector<std::string> filenames =	{"../../../files2read/r_eD-01.hipo","../../../files2read/r_eSn-01.hipo","../../../files2read/r_eSn-01.hipo","../../../files2read/r_eSn-02.hipo"};
    std::cout<< "Hello world \n";
    EventReader MC(filenames);
    Event test;
    for (int i=0; i<10; i++){
            test = MC.ProcessEventsInFile();

            PrintEventInfo(test);

    }    


    return 0;
}
