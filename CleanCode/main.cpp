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
"/work/clas12/dmat/simus/Deut/01/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/02/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/03/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/04/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/05/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/06/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/07/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/08/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/09/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/10/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/11/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/12/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/13/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/14/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/15/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/16/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/17/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/18/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/19/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/20/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/21/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/22/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/23/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/24/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/25/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/26/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/27/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/28/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/29/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/30/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/31/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/32/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/33/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/34/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/35/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/36/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/37/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/38/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/39/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/40/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/41/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/42/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/43/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/44/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/45/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/46/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/47/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/48/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/49/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/50/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/51/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/52/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/53/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/54/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/55/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/56/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/57/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/58/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/59/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/60/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/61/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/62/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/63/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/64/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/65/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/66/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/67/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/68/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/69/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/70/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/71/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/72/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/73/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/74/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/75/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/76/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/77/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/78/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/79/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/80/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/81/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/82/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/83/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/84/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/85/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/86/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/87/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/88/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/89/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/90/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/91/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/92/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/93/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/94/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/95/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/96/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/97/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/98/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/99/sidis_mc-master/r_ttest3D.hipo",
"/work/clas12/dmat/simus/Deut/100/sidis_mc-master/r_ttest3D.hipo"};
    //"/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/rec_clas_018325.evio.00010-00014.hipo",
    //"/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/rec_clas_018325.evio.00015-00019.hipo",
    //"/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/rec_clas_018325.evio.00020-00024.hipo",
    //"/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/rec_clas_018325.evio.00025-00029 .hipo"};
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
    ccuts.SetCutQ(0, 10);
    ccuts.SetCutY(0, 10);
    ccuts.SetCutW(10, 0);
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
    for (int i=0; i<30000; i++){
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