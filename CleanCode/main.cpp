#include <iostream>
#include <vector>
#include "reader.h"
#include "Event.h"
#include "EventReader.h"
#include "CutSet.h"
#include "Monitoring.h"
#include "Ratio.h"


//compile :
//g++ -o analyse main.cpp `root-config --cflags --libs` -I/home/matamoros/Desktop/localreader/hipo4/
// or use makefile 



int main() {
    std::vector<std::string> filenamesD = {"../../../files2read/r_eD-01.hipo", "../../../files2read/r_eD-01.hipo", "../../../files2read/r_eD-02.hipo"};
    std::vector<std::string> filenamesSn = {"../../../files2read/r_eSn-01.hipo", "../../../files2read/r_eSn-01.hipo", "../../../files2read/r_eSn-02.hipo"};
    
    EventReader MC_D(filenamesD);
    EventReader MC_Sn(filenamesSn);
    //std::optional<Event> testD;
    //std::optional<Event> testSn;


    //call class cuts 4 both nuclei
    CutSet cutD;
    CutSet cutSn;

    cutD.SetCutQ(1.5, 10);
    cutD.SetCutY(0.25, 0.85);
    cutD.SetCutW(0, 30);
    cutD.SetCutZ(0.3, 0.7);
    //both cuts are thesame 
    cutSn.SetCutQ(1.5, 10); 
    cutSn.SetCutY(0.25, 0.85);
    cutSn.SetCutW(0, 30);
    cutSn.SetCutZ(0.3, 0.7);
    int counter_elD= 0;
    int counter_elSn= 0;
    //call class ratio 
    Ratio ratiotest(cutD, cutSn);

    for (int i = 0; i < 20000; i++) {
        std::optional<Event> eventD;
        std::optional<Event> eventSn;

        //process events for deuterium and tin separately
        eventD = MC_D.ProcessEventsInFile();
        eventSn = MC_Sn.ProcessEventsInFile();

        if (eventD.has_value()) {
            //not using here the if hasvalue== false continue 
            counter_elD++;
            Event eventtestD =  eventD.value();
            eventtestD.calcAll();

            //ratiotest.FillHistograms(eventD.value(), "D");
            ratiotest.FillHistograms(eventtestD, "D");
            //ratiotest.calcR();

        }

        if (eventSn.has_value()) {
            //not using here the if hasvalue== false continue 
            counter_elSn++;
            Event eventtestSn =  eventSn.value();
            eventtestSn.calcAll();
            ratiotest.FillHistograms(eventtestSn, "Sn");
            //ratiotest.calcR();
        }
        ratiotest.calcR();
    }

    //ratiotest.WriteHistos("output_ratio_SnetD.root");

    return 0;
}
/*
int main() {
    std::vector<std::string> filenames =	{"../../../files2read/r_eD-01.hipo","../../../files2read/r_ttest3S.hipo","../../../files2read/r_ttest3S.hipo","../../../files2read/r_ttest3S.hipo"};
    std::cout<< "Hello world \n";
    EventReader MC(filenames);
    EventReader MC_D(filenamesD);
    EventReader MC_Sn(filenamesSn);
    std::optional<Event> test;
    std::optional<Event> testD;
    std::optional<Event> testSn;
    CutSet ccuts;
    CutSet bcuts;
    bcuts.SetCutQ(1.5,10);
    bcuts.SetCutY(0.25, 0.85);
    bcuts.SetCutW(0,30);
    bcuts.SetCutZ(0.3,0.7);
    ccuts.SetCutQ(0, 10);
    ccuts.SetCutY(0, 10);
    ccuts.SetCutW(0, 30);
    ccuts.SetCutZ(0, 10);

    Monitoring mon(ccuts);
    int sumevts = 0;
    Monitoring monb(bcuts);
//    for (auto & f : filenames ){
//        std::cout<< MC.getevtnbr()<< std::endl;
//        sumevts += MC.getevtnbr();
//
//    }

    int counter_el= 0.0;
    int counter_elD= 0.0;
    int counter_elSn= 0.0;


    for (int i=0; i<20000; i++){
            test = MC.ProcessEventsInFile();
           if (test.has_value()==false) continue;
            counter_el ++;
           Event eventtest =  test.value();
            eventtest.calcAll();
            //eventtest.Print();
            //monitoring(test,)
                mon.FillHistograms(eventtest);
                monb.FillHistograms(eventtest);
                //loop in hadron "list"
    }
    std::cout<<counter_el<<std::endl;  
    mon.WriteHistogramsToFile("output_ccuts.root");
    monb.WriteHistogramsToFile("output_bcuts.root");
    monb.DrawHistograms("after_bcuts.root");
    bcuts.Chop("chop_bcuts.root");
    bcuts.DrawChop("chopped_bcuts");


    return 0;
}
*/