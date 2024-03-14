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
#include "Dpt.h"
#include "cratio.h"
#include "FilePathGenerator.h"


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
    FilePathGenerator files;
    std::vector<std::string> filenamesLD2;
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/", filenamesLD2);
    std::vector<std::string> filenamesCuSn;
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesCuSn);


    std::cout<< "Hello world \n";
    //EventReader MC(filenames);
    EventReader MC_LD2(filenamesLD2);
    EventReader MC_CuSn(filenamesCuSn);
    //std::optional<Event> test;
    std::optional<Event> testLD2;
    std::optional<Event> testCuSn;
    CutSet Sncuts;   //Sn
    CutSet Cucuts;   //Cu
    CutSet LDcuts;   //LD2

    Sncuts.SetCutGen4Rat();
    Sncuts.SetCutVz(-3.5,-1.5);     //vz cut for Sn filter in the double target
    Cucuts.SetCutGen4Rat();
    Cucuts.SetCutVz(-8.5,-6.5);     //vz cut for Cu filter in the double target
    LDcuts.SetCutGen4Rat();
    LDcuts.SetCutVz(-7.5,-2.5);     //vz cut for precision in LD2 target


    //bcuts.SetCutPt2(3,10);        //this cut has not been added yet to passcuts
    int sumevts = 0;
    Monitoring monSn(Sncuts, "Sn");
    Monitoring monLD(LDcuts, "LD2");    // Modify the class so it doesnt have to take target name. (get targetevent) TBD!
    Monitoring monCu(Cucuts, "Cu");
    Ratio rat(LDcuts, Sncuts, "Sn"); //calling the class with the corresponding cuts
    //deltaptsq dpt(LDcuts, Sncuts);  
    //cratio crat(LDcuts, Sncuts);   
    Ratio rat2(LDcuts, Cucuts, "Cu"); // RE calling the class with the corresponding cuts for study with Cu
    //deltaptsq dpt2(LDcuts, Cucuts);
    //cratio crat2(LDcuts, Cucuts); 


    

    //add only 1 cut, contains all targets
//    for (auto & f : filenames ){
//        try{ 
//        std::cout << "Processing file: " << f << std::endl; 
//      
//        std::cout<< MC.getevtnbr()<< std::endl;
//        sumevts += MC.getevtnbr();
//        } catch (const std::exception &e) { 
//            std::cerr << "Error opening or processing file: " << f << std::endl;
//            std::cerr << "[ERROR] something went wrong with opening file: " << f << std::endl;
//            std::cerr << e.what() << std::endl; 
//            continue; // skips to next file if error
//        }
//    }
//    std::cout << "Total number of events: " << sumevts << std::endl;        //this line counts events if we run the commented (up) lines first


    int counter_el= 0.0;
    int counter_elLD2 = 0;
    int counter_elSn= 0.0;
    for (int i=0; i<900000; i++){
            //std::optional<Event> 
            testLD2 = MC_LD2.ProcessEventsInFile();
            //std::optional<Event> 
            testCuSn = MC_CuSn.ProcessEventsInFile();
            //test = MC.ProcessEventsInFile();
            if (testLD2.has_value()) {
                counter_elLD2++;
                Event eventtestLD2 = testLD2.value();
                eventtestLD2.SetTargetType(0);
                eventtestLD2.calcAll();
                //LDcuts.SetCutGentest(eventtestLD2);
                //monLD.FillHistograms(eventtestLD2);       //WHAT is this again ? 
                monLD.FillHistogramsNoCuts(eventtestLD2);
                //monLD.FillHistogramswCuts(eventtestLD2);
                rat.FillHistograms(eventtestLD2);  
                rat2.FillHistograms(eventtestLD2);
                //dpt.FillHistograms(eventtestLD2);
                //crat.FillHistograms(eventtestLD2);
            }
            if (testCuSn.has_value()) {
                counter_elSn++;
                Event eventtestCuSn = testCuSn.value();
                eventtestCuSn.SetTargetType(1);
                eventtestCuSn.calcAll();

                //Sncuts.SetCutGentest(eventtestCuSn);
                //monSn.FillHistograms(eventtestCuSn);
                monSn.FillHistogramsNoCuts(eventtestCuSn);
                monCu.FillHistogramsNoCuts(eventtestCuSn);
                //monSn.FillHistogramswCuts(eventtestCuSn);
                rat.FillHistograms(eventtestCuSn);
                rat2.FillHistograms(eventtestCuSn);

                //dpt.FillHistograms(eventtestCuSn);
                //crat.FillHistograms(eventtestCuSn);
            }

           
           //if (test.has_value()==false) continue;
           // counter_el ++;
           //Event eventtest =  test.value();
           // eventtest.calcAll();
           // //eventtest.Print();
           // //monitoring(test,)
           //     monSn.FillHistograms(eventtest);
           //     //rat.FillHistograms( )
           //     //loop in hadron "list"
    
    
    
    }

    std::cout << "Events processed for LD2: " << counter_elLD2 << std::endl;
    std::cout << "Events processed for Sn: " << counter_elSn << std::endl;

    //monLD.WriteHistogramsToFile("output_LD2.root");
    //monSn.WriteHistogramsToFile("output_CuSn.root");
    monLD.DrawHistograms("noCutsVarLD2");
    monLD.DrawCaloHistograms("REwCutscaloLD2");
    //monLD.DrawHelicityHistograms("REwCutshelicityLD2");
    //monLD.DrawCherenkovHistograms("REwCutscherenkovLD2");
    //monLD.DrawMomentumHistograms("noCutsmomentumLD2");
    //monSn.DrawMomentumElectronHistograms("noCutsmomentumElectronLD2Sn");
    //monLD.DrawMomentumElectronHistograms("noCutsmomentumElectronLD2");
    monLD.DrawVertexHistograms("noCutsVertexLD2");
    //monLD.DrawMomentumHadronHistograms("noCutsmomentumHadronLD2");
    //monLD.DrawEnergyHistograms("noCutsEnergyLD2");

    //monSn.DrawHistograms("after_cuts_CuSn");
    //monLD.DrawHistogramsPos("comp2D");  
    //monLD.DrawR_Histograms("RmonitoringLD2");

    rat.calcR();
    //rat.writeMatrixToFile("matrix_test.txt");
    rat2.calcR();
    rat2.writeMatrixToFile("matrix2_test.txt");
    //rat.multiplotR();
    //rat2.multiplotR();
    rat.multiplotR(rat2);

    //dpt.calcDpt();
        //dpt.writeMatrixToFile("matrix_Dpt.txt");
    //dpt.multiplotDpt();
    //crat.calcCratio();
        //crat.writeMatrixToFile("matrix_Cratio.txt");
    //crat.multiplotCratio();
    //


    return 0;
}