#include <iostream>
#include <vector>
#include <string>
#include <vector>

#include "reader.h"
#include "Event.h"
#include "EventReader.h"
#include "CutSet.h"
#include "Monitoring.h"
#include "Monunfold.h"
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
    std::vector<std::string> simufilesCxC;

    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/", filenamesLD2);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesCuSn);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesCxC);
    ////files.ParDir2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/", simufilesLD2);   //do not uncomment this line
    
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/Deut", simufilesLD2);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/C", simufilesCxC);  //uncomment for sim
    
    //test on ifarm 4 self on LD2 true
    //files.Files2Vector("/volatile/clas12/rg-d/production/prod/v4_ob_LD2/dst/recon/018869/", simufilesCxC);
    
    //files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder", simufilesSn);  //uncomment for sim



    //Uncomment 4 test on ifarm, comment all above
    //files.Files2Vector("/volatile/clas12/rg-d/production/prod/v3_ob_LD2/dst/recon/019033/", filenamesLD2);
    //files.Files2Vector("/volatile/clas12/rg-d/production/prod/v3_ob_CuSn/dst/recon/019130/", filenamesCuSn);
    //files.Files2Vector("/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018835/", filenamesCxC);
    //uncommment also below for sim  on ifarm
    //files.SnDir2Vector("/volatile/clas12/dmat/gen/Deut/", simufilesLD2);  //uncomment for sim
    //files.SnDir2Vector("/volatile/clas12/dmat/gen/C/", simufilesCxC);  //uncomment for sim 


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
    //EventReader Sim_CuSn(simufilesSn);
    EventReader Sim_CxC(simufilesCxC);
    


    //std::optional<Event> test;
    std::optional<Event> testLD2;
    std::optional<Event> testCuSn;
    std::optional<Event> testCxC;
    std::optional<Event> simLD2;
    //std::optional<Event> simCuSn;
    std::optional<Event> simCxC;

    CutSet Sncuts;   //Sn
    CutSet Cucuts;   //Cu
    CutSet LDcuts;   //LD2
    CutSet CCcuts;   //CxC
    CutSet C1cuts;   //C1
    CutSet C2cuts;   //C2
    //CutSet simSncuts;   //Sn
    //CutSet simCucuts;   //Cu    
    CutSet simLDcuts;   //LD2
    //CutSet simCCcuts;   //CxC
    CutSet simC1cuts;   //C1
    CutSet simC2cuts;   //C2

    Sncuts.SetCutGen4Rat();
    Sncuts.SetCutVz(-3.5,-1.5);     //vz cut for Sn filter in the double target
    Cucuts.SetCutGen4Rat();
    Cucuts.SetCutVz(-8.5,-6.5);     //vz cut for Cu filter in the double target
    LDcuts.SetCutGen4Rat();
    LDcuts.SetCutVz(-7.5,-2.5);     //vz cut for precision in LD2 target
    CCcuts.SetCutGen4Rat();
    CCcuts.SetCutVz(-10.0,0.0);     //general cut, in order to include both Carbon positions in Vz (initial)
    C1cuts.SetCutGen4Rat();
    C1cuts.SetCutVz(-3.5,-1.5);     //vz cut for C1 target
    C2cuts.SetCutGen4Rat();
    C2cuts.SetCutVz(-8.5,-6.5);     //vz cut for C2 target
    //simSncuts.SetCutGen4Rat();
    //simSncuts.SetCutVz(-2.0,2.0);     //vz cut for Sn filter in the double target
    //simCucuts.SetCutGen4Rat();
    //simCucuts.SetCutVz(-8.5,-6.5);     //vz cut for Cu filter in the double target
    simLDcuts.SetCutGen4Rat();
    simLDcuts.SetCutVz(-2.0,2.0);     //vz cut for precision in LD2 target
    //simCCcuts.SetCutGen4Rat();
    // simCCcuts.SetCutVz(-8.5,-6.5);     //vz cut for CxC target
    simC1cuts.SetCutGen4Rat();
    simC1cuts.SetCutVz(-3,-2);     //vz cut for C1 target
    simC2cuts.SetCutGen4Rat();
    simC2cuts.SetCutVz(-8,-7);     //vz cut for C2 target
    //Fix target borders in constants. remember is not the definite vz cut, but a cut for target separation 



    //bcuts.SetCutPt2(3,10);        //this cut has not been added yet to passcuts
    int sumevts = 0;
    Monitoring monSn(Sncuts, "Sn");
    Monitoring monLD(LDcuts, "LD2");    // Modify the class so it doesnt have to take target name. (get targetevent) TBD!
    Monitoring monCu(Cucuts, "Cu");
    Monitoring monC1(C1cuts, "C1");
    Monitoring monC2(C2cuts, "C2");
    Monunfold munfLD(LDcuts, "LD2b");
    Monunfold munfC1(C1cuts, "C1");
    Monunfold munfC2(C2cuts, "C2");

    Monitoring monSimLD(simLDcuts, "LD2_sim");
    Monunfold munfSimCC(CCcuts, "CxC_sim");
    //Monitoring monSimSn(simSncuts, "Sn_sim");
    //Monitoring monSimCu(simCucuts, "Cu");
    Monitoring monSimC1(simC1cuts, "C1_sim");     //This needs to be figured out ASAP
    //Monitoring monSimC2(simC2cuts, "C2");     //This needs to be figured out ASAP
    Ratio rat(LDcuts, Sncuts, "Sn"); //calling the class with the corresponding cuts
    Ratio rat2(LDcuts, Cucuts, "Cu"); // RE calling the class with the corresponding cuts for study with Cu
    Ratio rat3(LDcuts, C1cuts, "C1"); // RE calling the class with the corresponding cuts for study with C1
    Ratio rat4(LDcuts, C2cuts, "C2"); // RE calling the class with the corresponding cuts for study with C2
    Ratio simrat3( simLDcuts, simC1cuts, "C1_sim"); //calling the class with the corresponding cuts
    Ratio simrat4(simLDcuts, simC2cuts, "C2_sim"); // RE calling the class with the corresponding cuts for study with Cu

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
    int counter_elsimLD2= 0.0;
    int totalevts = 10000;
    for (int i=0; i<totalevts; i++){
            //std::optional<Event> 
            //testLD2 = MC_LD2.ProcessEventsInFile();
            //testCuSn = MC_CuSn.ProcessEventsInFile();
            //testCxC = MC_CxC.ProcessEventsInFile();
            //simLD2 = Sim_LD2.ProcessEventsInFile();
            //simCuSn = Sim_CuSn.ProcessEventsInFile();
            simCxC = Sim_CxC.ProcessEventsInFile();

            if ( simCxC.has_value()) {
                ////std::cout << "entering Simulated data " << std::endl;
                Event eventsimCxC = simCxC.value();
                eventsimCxC.SetTargetType(1);
                eventsimCxC.calcAll();
                //eventsimCxC.calcMCAll();
                ////monSimC1.FillHistogramsNoCuts(eventsimCxC);
                //munfSimCC.Fill(eventsimCxC);
                //monSimC1.FillHistogramsNoCutsMC(eventsimCxC);
                ////monSimCC.FillHistogramswCuts(eventsimCxC);
                simrat3.FillHistograms(eventsimCxC);
                simrat4.FillHistograms(eventsimCxC);
            }
            files.displayProgress(i + 1, totalevts);
            

           
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
    std::cout << "\nProcessing completed \n";

    std::cout << "Events processed for LD2: " << counter_elLD2 << std::endl;
    std::cout << "Events processed for Sn: " << counter_elSn << std::endl;
    std::cout << "Events processed for CxC: " << counter_elCxC << std::endl;
    std::cout << "Events processed for simLD2: " << counter_elsimLD2 << std::endl;
    std::cout << "Events processed for simCxC: " << counter_elCxC << std::endl;

    //monLD.WriteHistogramsToFile("output_LD2.root");
    //monSn.WriteHistogramsToFile("output_CuSn.root");
    munfSimCC.DrawHistograms("testsimCxC");
    munfSimCC.DrawCompRECMC("compsimCxC");
    monSimC1.DrawCompRECMC("compsimC1_plot2D");
    munfC1.DrawOnlyVz(munfC2, "VzC1C2");
    munfC2.DrawOnlyVz(munfC1, "VzC2C1");
    munfLD.DrawOnlyVz(munfC1, "VzLD2C1");
    simrat3.calcRcarbon(simrat4);
    simrat3.multiplotRbis();

    return 0;
}
