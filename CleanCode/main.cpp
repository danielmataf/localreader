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
#include "sratio.h"
#include "c2ratio.h"
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
    std::vector<std::string> simufilesLD2;
    std::vector<std::string> simufilesSn;
    std::vector<std::string> filenamesCuSn;
    std::vector<std::string> filenamesCxC;

    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/", filenamesLD2);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesCuSn);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesCxC);
    files.ParDir2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/", simufilesLD2);
    
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder", simufilesLD2);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder", simufilesSn);  //uncomment for sim



    //Uncomment 4 test on ifarm, comment all above
    //files.Files2Vector("/volatile/clas12/rg-d/production/prod/v3_ob_LD2/dst/recon/019033/", filenamesLD2);
    //files.Files2Vector("/volatile/clas12/rg-d/production/prod/v3_ob_CuSn/dst/recon/019130/", filenamesCuSn);
    //files.Files2Vector("/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018835/", filenamesCxC);
    //uncommment also below for sim  on ifarm
    //files.SnDir2Vector("/volatile/clas12/dmat/gen/Deut/", simufilesLD2);  //uncomment for sim
    //files.SnDir2Vector("/volatile/clas12/dmat/gen/Sn/", simufilesSn);  //uncomment for sim 


    //files.SnDir2Vector("/volatile/clas12/dmat/gen/Sn/", simufilesSn);  //uncomment for sim 


    //std::vector<std::string> filenamesSimLD2;
    //files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/", filenamesSimLD2);
    //std::vector<std::string> filenamesSimCuSn;
    //files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesSimCuSn);
    //std::vector<std::string> filenamesSimCxC;
    //files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesSimCxC);


    std::cout<< "Hello world \n";
    //EventReader MC(filenames);
    EventReader MC_LD2(filenamesLD2);
    EventReader MC_CuSn(filenamesCuSn);
    EventReader MC_CxC(filenamesCxC);
    EventReader Sim_LD2(simufilesLD2);
    EventReader Sim_CuSn(simufilesSn);
    //EventReader Sim_CxC(filenamesSimCxC);


    //std::optional<Event> test;
    std::optional<Event> testLD2;
    std::optional<Event> testCuSn;
    std::optional<Event> testCxC;
    std::optional<Event> simLD2;
    std::optional<Event> simCuSn;
    //std::optional<Event> simCxC;

    CutSet Sncuts;   //Sn
    CutSet Cucuts;   //Cu
    CutSet LDcuts;   //LD2
    CutSet CCcuts;   //CxC
    CutSet simSncuts;   //Sn
    CutSet simCucuts;   //Cu    
    CutSet simLDcuts;   //LD2
    //CutSet simCCcuts;   //CxC

    Sncuts.SetCutGen4Rat();
    Sncuts.SetCutVz(-3.5,-1.5);     //vz cut for Sn filter in the double target
    Cucuts.SetCutGen4Rat();
    Cucuts.SetCutVz(-8.5,-6.5);     //vz cut for Cu filter in the double target
    LDcuts.SetCutGen4Rat();
    LDcuts.SetCutVz(-7.5,-2.5);     //vz cut for precision in LD2 target
    CCcuts.SetCutGen4Rat();
    CCcuts.SetCutVz(-8.5,-6.5);     //vz cut for CxC target
    simSncuts.SetCutGen4Rat();
    simSncuts.SetCutVz(-2.0,2.0);     //vz cut for Sn filter in the double target
    //simCucuts.SetCutGen4Rat();
    //simCucuts.SetCutVz(-8.5,-6.5);     //vz cut for Cu filter in the double target
    simLDcuts.SetCutGen4Rat();
    simLDcuts.SetCutVz(-2.0,2.0);     //vz cut for precision in LD2 target
    //simCCcuts.SetCutGen4Rat();
    //simCCcuts.SetCutVz(-8.5,-6.5);     //vz cut for CxC target
    //Fix target borders in constants. remember is not the definite vz cut, but a cut for target separation 



    //bcuts.SetCutPt2(3,10);        //this cut has not been added yet to passcuts
    int sumevts = 0;
    Monitoring monSn(Sncuts, "Sn");
    Monitoring monLD(LDcuts, "LD2");    // Modify the class so it doesnt have to take target name. (get targetevent) TBD!
    Monitoring monCu(Cucuts, "Cu");
    Monitoring monC1(CCcuts, "C1");
    Monitoring monC2(CCcuts, "C2");
    Monitoring monSimLD(simLDcuts, "LD2_sim");
    Monitoring monSimSn(simSncuts, "Sn_sim");
    //Monitoring monSimCu(simCucuts, "Cu");
    //Monitoring monSimCC(simCCcuts, "CxC");
    Ratio rat(LDcuts, Sncuts, "Sn"); //calling the class with the corresponding cuts
    Ratio rat2(LDcuts, Cucuts, "Cu"); // RE calling the class with the corresponding cuts for study with Cu
    Ratio rat3(LDcuts, CCcuts, "C1"); // RE calling the class with the corresponding cuts for study with C1
    Ratio rat4(LDcuts, CCcuts, "C2"); // RE calling the class with the corresponding cuts for study with C2
    Ratio simrat( simLDcuts, simSncuts, "Sn_sim"); //calling the class with the corresponding cuts
    Ratio simrat2(simLDcuts, simCucuts, "Cu_sim"); // RE calling the class with the corresponding cuts for study with Cu

    //Ratio simrat3(simLDcuts, simCCcuts, "CxC"); // RE calling the class with the corresponding cuts for study with CxC
    deltaptsq dpt(LDcuts, Sncuts, "Sn");  
    deltaptsq dpt2(LDcuts, Cucuts, "Cu");
    deltaptsq dpt3(LDcuts, CCcuts, "CxC");
    //deltaptsq simdpt( simLDcuts, simSncuts, "Sn");
    //deltaptsq simdpt2(simLDcuts, simCucuts, "Cus");
    //deltaptsq simdpt3(simLDcuts, simCCcuts, "CxC");

    cratio crat(LDcuts, Sncuts, "Sn");   
    cratio crat2(LDcuts, Cucuts, "Cu");
    cratio crat3(LDcuts, CCcuts, "CxC"); 
    //cratio simcrat( simLDcuts, simSncuts, "Sn");
    //cratio simcrat2(simLDcuts, simCucuts, "Cu");
    //cratio simcrat3(simLDcuts, simCCcuts, "CxC");
    sratio srat (LDcuts, Sncuts, "Sn");
    sratio srat2 (LDcuts, Cucuts, "Cu");
    sratio srat3 (LDcuts, CCcuts, "CxC");
    //sratio simsrat ( simLDcuts, simSncuts, "Sn");
    //sratio simsrat2 (simLDcuts, simCucuts, "Cu");
    //sratio simsrat3 (simLDcuts, simCCcuts, "CxC");
    c2ratio c2rat (LDcuts, Sncuts, "Sn");
    c2ratio c2rat2 (LDcuts, Cucuts, "Cu");
    c2ratio c2rat3 (LDcuts, CCcuts, "CxC");
    //c2ratio simc2rat ( simLDcuts, simSncuts, "Sn");
    //c2ratio simc2rat2 (simLDcuts, simCucuts, "Cu");
    //c2ratio simc2rat3 (simLDcuts, simCCcuts, "CxC");



    

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
    int counter_elCxC= 0.0;
    for (int i=0; i<90000; i++){
            //std::optional<Event> 
            testLD2 = MC_LD2.ProcessEventsInFile();
            testCuSn = MC_CuSn.ProcessEventsInFile();
            testCxC = MC_CxC.ProcessEventsInFile();
            simLD2 = Sim_LD2.ProcessEventsInFile();
            simCuSn = Sim_CuSn.ProcessEventsInFile();
            //simCxC = Sim_CxC.ProcessEventsInFile();

            //test = MC.ProcessEventsInFile();
            if (testLD2.has_value()) {
                counter_elLD2++;
                Event eventtestLD2 = testLD2.value();
                
                eventtestLD2.SetTargetType(0);
                //eventsimLD2.SetTargetType(0);
                
                eventtestLD2.calcAll();
                //eventsimLD2.calcAll();

                //LDcuts.SetCutGentest(eventtestLD2);
                monLD.FillHistograms(eventtestLD2);       //Re comment   this again ? 
                //monLD.FillHistogramsNoCuts(eventtestLD2);
                //monLD.FillHistogramswCuts(eventtestLD2);
                rat.FillHistograms(eventtestLD2);  
                rat2.FillHistograms(eventtestLD2);
                rat3.FillHistograms(eventtestLD2);
                //simrat.FillHistograms(eventtestLD2);
                //simrat2.FillHistograms(eventtestLD2);
                dpt.FillHistograms(eventtestLD2);
                //dpt2.FillHistograms(eventtestLD2);
                //dpt3.FillHistograms(eventtestLD2);
                //crat.FillHistograms(eventtestLD2);
                //crat2.FillHistograms(eventtestLD2);
                //crat3.FillHistograms(eventtestLD2);
                //srat.FillHistograms(eventtestLD2);
                //srat2.FillHistograms(eventtestLD2);
                //srat3.FillHistograms(eventtestLD2);
                //c2rat.FillHistograms(eventtestLD2);
                //c2rat2.FillHistograms(eventtestLD2);
                //c2rat3.FillHistograms(eventtestLD2);
                //rat, dpt , crat, srat and c2rat are all observables 
            }
            if (testCuSn.has_value()) {
                counter_elSn++;
                Event eventtestCuSn = testCuSn.value();
                //Event eventsimCuSn = simCuSn.value();
                
                eventtestCuSn.SetTargetType(1);
                //eventsimCuSn.SetTargetType(1);
                
                eventtestCuSn.calcAll();
                //eventsimCuSn.calcAll();

                //Sncuts.SetCutGentest(eventtestCuSn);
                monSn.FillHistograms(eventtestCuSn);    //recomment this 
                monCu.FillHistograms(eventtestCuSn);    //recomment this 
                //monSn.FillHistogramsNoCuts(eventtestCuSn);
                //monCu.FillHistogramsNoCuts(eventtestCuSn);
                //monSn.FillHistogramswCuts(eventtestCuSn);
                rat.FillHistograms(eventtestCuSn);
                rat2.FillHistograms(eventtestCuSn);
                //simrat.FillHistograms(eventtestCuSn);
                //simrat2.FillHistograms(eventtestCuSn);

                //dpt.FillHistograms(eventtestCuSn);
                //dpt2.FillHistograms(eventtestCuSn);
//
                //crat.FillHistograms(eventtestCuSn);
                //crat2.FillHistograms(eventtestCuSn);
//
                //srat.FillHistograms(eventtestCuSn);
                //srat2.FillHistograms(eventtestCuSn);
//
                //c2rat.FillHistograms(eventtestCuSn);
                //c2rat2.FillHistograms(eventtestCuSn);
            }
            if (testCxC.has_value()) {
                counter_elCxC++;
                Event eventtestCxC = testCxC.value();
                //Event eventsimCxC = simCxC.value();
                eventtestCxC.SetTargetType(1);              //is it targettype 1 or targettype 2? needs to be checked and modified accordingly

                eventtestCxC.calcAll();
                //eventsimCxC.calcAll();
                monC1.FillHistograms(eventtestCxC);     //recomment this
                rat3.FillHistograms(eventtestCxC);
                
            }
            //Simulated data handler
            if (simLD2.has_value()) {
                Event eventsimLD2 = simLD2.value();
                eventsimLD2.SetTargetType(0);
                eventsimLD2.calcAll();
                //monSimLD.FillHistograms(eventsimLD2);
                monSimLD.FillHistogramsNoCuts(eventsimLD2);
                //monSimLD.FillHistogramswCuts(eventsimLD2);
                simrat.FillHistograms(eventsimLD2);  
                simrat2.FillHistograms(eventsimLD2);

            }
            if (simCuSn.has_value()) {
                Event eventsimCuSn = simCuSn.value();
                eventsimCuSn.SetTargetType(1);
                eventsimCuSn.calcAll();
                //monSimSn.FillHistograms(eventsimCuSn);
                monSimSn.FillHistogramsNoCuts(eventsimCuSn);
                //monSimSn.FillHistogramswCuts(eventsimCuSn);
                simrat.FillHistograms(eventsimCuSn);
                simrat2.FillHistograms(eventsimCuSn);
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
    std::cout << "Events processed for CxC: " << counter_elCxC << std::endl;

    //monLD.WriteHistogramsToFile("output_LD2.root");
    //monSn.WriteHistogramsToFile("output_CuSn.root");
    monLD.DrawHistograms("noCutsVarLD2");
    monSimLD.DrawHistTrueandSIM(monLD,"testsimLD2");
    
    
    //monC1.DrawHistograms("dddd");   //use this to check Vz for CxC 
    //monLD.DrawCaloHistograms("REwCutscaloLD2");
    //monLD.DrawHelicityHistograms("REwCutshelicityLD2");
    //monSimLD.DrawHistograms("noCutsSimVarLD2");
    //monSimSn.DrawHistograms("noCutsSimVarSn");    
//monLD.DrawCherenkovHistograms("REwCutscherenkovLD2");
    //monLD.DrawMomentumHistograms("noCutsmomentumLD2");
    //monSn.DrawMomentumElectronHistograms("noCutsmomentumElectronLD2Sn");
    //monLD.DrawMomentumElectronHistograms("noCutsmomentumElectronLD2");
    //monLD.DrawVertexHistograms("noCutsVertexLD2");
    //monC1.DrawVertexHistograms("noCutsVertexCxC");
    //monSimLD.DrawVertexHistograms("noCutsVertexSimLD2");
    //monSimSn.DrawVertexHistograms("noCutsVertexSimSn");
    //monLD.DrawMomentumHadronHistograms("noCutsmomentumHadronLD2");
    //monLD.DrawEnergyHistograms("noCutsEnergyLD2");

    //monSn.DrawHistograms("after_cuts_CuSn");
    //monLD.DrawHistogramsPos("LD2","comp2D_LD");    //  plots for RD
    //monSn.DrawHistogramsPos("Sn","comp2D_Sn"); //  plots for RD
    //monCu.DrawHistogramsPos("Cu","comp2D_Cu"); //  plots for RD
    //monCC.DrawHistogramsPos("CxC","comp2DCxC");    //  plots for RD  
    //monLD.DrawR_Histograms("RmonitoringLD2");

    rat.calcR();
    //rat.writeMatrixToFile("matrix_test.txt");
    rat2.calcR();
    //rat2.writeMatrixToFile("matrix2_test.txt");
    rat3.calcR();
    simrat2.calcR();
    simrat.calcR(); 



    //simrat.multiplotR(simrat2,simrat);    
    //rat.multiplotR(rat2, rat3, simrat, simrat2);  //uncomment for sim 

    //dpt.calcDpt();
    //dpt2.calcDpt();
    //dpt3.calcDpt();
    ////dpt.writeMatrixToFile("matrix_Dpt.txt");
    //dpt.multiplotDpt(dpt2, dpt3);
    
    //crat.calcCratio();
    //crat2.calcCratio();
    //crat3.calcCratio();
    ////crat.writeMatrixToFile("matrix_Cratio.txt");
    //crat.multiplotCratio( crat2, crat3);
    
    //srat.calcSratio();
    //srat2.calcSratio();
    //srat3.calcSratio();
    //srat.multiplotSratio(srat2, srat3);

    //c2rat.calcC2ratio();
    //c2rat2.calcC2ratio();
    //c2rat3.calcC2ratio();
    //c2rat.multiplotC2ratio(c2rat2, c2rat3);
    //


    return 0;
}
