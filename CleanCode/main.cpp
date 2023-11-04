#include <iostream>
#include <vector>
#include <string>
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


/*
int main() {
    std::vector<std::string> filenamesD = {"../../../files2read/r_eD-01.hipo", "../../../files2read/r_eD-01.hipo", "../../../files2read/r_eD-02.hipo"};
    std::vector<std::string> filenamesSn = {"../../../files2read/r_eSn-01.hipo", "../../../files2read/r_eSn-01.hipo", "../../../files2read/r_eSn-02.hipo"};
    
    EventReader MC_D(filenamesD);
    EventReader MC_Sn(filenamesSn);
    //std::optional<Event> testD;
    //std::optional<Event> testSn;
        std::optional<Event> eventD;
        std::optional<Event> eventSn;


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
    ratiotest.PlotRatio("test3D");

    return 0;
}
*/

int main() {
        std::vector<std::string> filenames = {
    "/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/rec_clas_018428.evio.00010-00014.hipo",
    "/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/rec_clas_018428.evio.00015-00019.hipo",
    "/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/rec_clas_018428.evio.00020-00024.hipo",
    "/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/rec_clas_018428.evio.00025-00029.hipo"};
    //std::vector<std::string> filenames = {"../../../files2read/r_eD-01.hipo", "../../../files2read/r_eD-01.hipo", "../../../files2read/r_eD-02.hipo"};
    //std::vector<std::string> filenamesSn = {"../../../files2read/r_eSn-01.hipo", "../../../files2read/r_eSn-01.hipo", "../../../files2read/r_eSn-02.hipo"};

    std::cout<< "Hello world \n";
    EventReader MC(filenames);
    std::optional<Event> test;
    std::optional<Event> testD;
    std::optional<Event> testSn;
    CutSet ccuts;
    CutSet bcuts;
    bcuts.SetCutQ(1.5,10);
    bcuts.SetCutY(0.25, 0.85);
    bcuts.SetCutW(6,30);
    bcuts.SetCutZ(0.3,0.7);
    bcuts.SetCutVz(-7.5,-2.5);
    ccuts.SetCutQ(0, 10);
    ccuts.SetCutY(0, 10);
    ccuts.SetCutW(10, 0);
    ccuts.SetCutZ(0, 10);
    Monitoring mon(ccuts);
    int sumevts = 0;
    Monitoring monb(bcuts);
    for (auto & f : filenames ){
        try{ 
        std::cout << "Processing file: " << f << std::endl; 
      
        std::cout<< MC.getevtnbr()<< std::endl;
        sumevts += MC.getevtnbr();
        } catch (const std::exception &e) { 
            std::cerr << "Error opening or processing file: " << f << std::endl;
            std::cerr << "[ERROR] something went wrong with opening file: " << f << std::endl;
            std::cerr << e.what() << std::endl; 
            continue; // Skip to the next file in case of an error
        }
    }
    //std::cout << "Total number of events: " << sumevts << std::endl;
    //int counter_el= 0.0;
    //int counter_elD= 0.0;
    //int counter_elSn= 0.0;
    //for (int i=0; i<30000; i++){
    //        test = MC.ProcessEventsInFile();
    //       if (test.has_value()==false) continue;
    //        counter_el ++;
    //       Event eventtest =  test.value();
    //        eventtest.calcAll();
    //        //eventtest.Print();
    //        //monitoring(test,)
    //            mon.FillHistograms(eventtest);
    //            monb.FillHistograms(eventtest);
    //            //loop in hadron "list"
    //}
    //std::cout<<counter_el<<std::endl;  
    //mon.WriteHistogramsToFile("output_testVz.root");
    //monb.WriteHistogramsToFile("output_testVz.root");
    //monb.DrawHistograms("after_bcutsVz.root");
    //bcuts.Chop("chop_bcutsVz.root");
    //bcuts.DrawChop("chopped_bcutsVz");

    return 0;
}