#include <iostream>
#include <vector>
#include "reader.h"
#include "Event.h"
#include "EventReader.h"
#include "CutSet.h"
#include "Monitoring.h"

//compile :
//g++ -o analyse main.cpp `root-config --cflags --libs` -I/home/matamoros/Desktop/localreader/hipo4/
// or use makefile 
int main() {
    std::vector<std::string> filenames =	{"../../../files2read/r_eD-01.hipo","../../../files2read/r_ttest3S.hipo","../../../files2read/r_ttest3S.hipo","../../../files2read/r_ttest3S.hipo"};
    std::cout<< "Hello world \n";
    EventReader MC(filenames);
    std::optional<Event> test;
    CutSet ccuts;
    CutSet bcuts;
    bcuts.SetCutQ(1.5,10);
    bcuts.SetCutY(0.25, 0.85);
    bcuts.SetCutW(0,4);
    bcuts.SetCutZ(0.3,0.7);
    ccuts.SetCutQ(0, 10);
    ccuts.SetCutY(0, 10);
    ccuts.SetCutW(0, 100);
    ccuts.SetCutZ(0, 10);

    Monitoring mon(ccuts);
    int sumevts = 0;
    //Monitoring monb(bcuts);
//    for (auto & f : filenames ){
//        std::cout<< MC.getevtnbr()<< std::endl;
//        sumevts += MC.getevtnbr();
//
//    }

            int counter_el= 0.0;

    for (int i=0; i<30000; i++){
           //std::cout<< "Event number  " << i << std::endl;

            test = MC.ProcessEventsInFile();
           if (test.has_value()==false) continue;
            counter_el ++;
           Event eventtest =  test.value();
            eventtest.calcAll();
            //eventtest.Print();
            //monitoring(test,)
                mon.FillHistograms(eventtest);
                //monb.FillHistograms(test);
                //loop in hadron "list"
    }    
    //writing outside of the loop
    std::cout<<counter_el<<std::endl;  
    mon.WriteHistogramsToFile("output_ccuts.root");
    //monb.WriteHistogramsToFile("output_bcuts.root");


    return 0;
}
