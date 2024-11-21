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
#include "constants.h"  



int main() {
    FilePathGenerator files;
    std::vector<std::string> filenamesLD2;
    std::vector<std::string> simufilesLD2;
    std::vector<std::string> simufilesSn;
    std::vector<std::string> simufilesCu;
    std::vector<std::string> filenamesCuSn;
    std::vector<std::string> filenamesCxC;
    std::vector<std::string> simufilesCxC;

    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/0v6/", filenamesLD2);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/0v6/", filenamesCuSn);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/C_v0/0v6/", filenamesCxC);
    ////files.ParDir2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/", simufilesLD2);   //do not uncomment this line
    
    //files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/LD2", simufilesLD2);  //uncomment for sim
files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novLD2", simufilesLD2);  //uncomment for sim
    
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novC", simufilesCxC);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novCu", simufilesCu);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novSn", simufilesSn);  //uncomment for sim

    //Uncomment 4 test on ifarm, comment all above
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/skim_pass0v6/LD2/", filenamesLD2);
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/skim_pass0v6/CuSn/", filenamesCuSn);
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/skim_pass0v6/CxC/", filenamesCxC);
    ////uncommment also below for sim  on ifarm
    //files.SnDir2Vector("/volatile/clas12/dmat/test/cv_newrgd/", simufilesLD2);  //uncomment for sim
    //files.SnDir2Vector("/volatile/clas12/dmat/test/cv_newrgdCC/", simufilesCxC);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/test/cv_newrgdCu/", simufilesCu);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/test/cv_newrgdSn/", simufilesSn);  //uncomment for sim 


    std::cout<< "Analysis in progress... \n";
    EventReader Sim_CxC(simufilesCxC);
    EventReader RGD_CxC(filenamesCxC);
    EventReader Sim_LD2(simufilesLD2);
    EventReader RGD_LD2(filenamesLD2);   
    EventReader Sim_Cu(simufilesCu);
    EventReader Sim_Sn(simufilesSn);
    EventReader RGD_CuSn(filenamesCuSn);   



    std::optional<Event> simCxC;
    std::optional<Event> simCxC_MC; //testing this if we can loop twice top read MC events//
    std::optional<Event> testCxC;
    std::optional<Event> simLD2;
    std::optional<Event> simLD2_MC;
    std::optional<Event> testLD2;
    std::optional<Event> simCu;
    std::optional<Event> simCu_MC;
    std::optional<Event> testCu;
    std::optional<Event> simSn;
    std::optional<Event> simSn_MC;
    std::optional<Event> testSn;



    CutSet testC1cuts;  //C1 RGD
    CutSet simC1cuts;   //C1
    CutSet testC2cuts;  //C2 RGD  
    CutSet simC2cuts;   //C2
    CutSet testLD2cuts;
    CutSet simLD2cuts;
    CutSet testCucuts;
    CutSet simCucuts;
    CutSet testSncuts;
    CutSet simSncuts;



    //Defining cuts here. Vz cuts for SIM and RGD may differ -> check Mon Vz
    //simC1cuts.SetCutVz(Constants::RcutminVzC1sim,Constants::RcutminVzC1sim);     //vz cut for C1 target
    simC1cuts.SetCutVz(Constants::RcutminVzC1sim,Constants::RcutmaxVzC1sim);     //vz cut for C1 target
    simC1cuts.SetCutGen4Rat();
    //testC1cuts.SetCutVz(Constants::RcutminVzC1data,Constants::RcutminVzC1data);     //vz cut for C1 target
    testC1cuts.SetCutVz(Constants::RcutminVzC1data,Constants::RcutmaxVzC1data);     //vz cut for C1 target
    testC1cuts.SetCutGen4Rat();

    simC2cuts.SetCutVz(Constants::RcutminVzC2sim    , Constants::RcutmaxVzC2sim);     //vz cut for C2 target
    simC2cuts.SetCutGen4Rat();
    testC2cuts.SetCutVz(Constants::RcutminVzC2data, Constants::RcutmaxVzC2data);    //vz cut for C2 target
    testC2cuts.SetCutGen4Rat();

    simLD2cuts.SetCutVz(Constants::RcutminVzLD2,Constants::RcutmaxVzLD2);     //vz cut for LD2 target
    simLD2cuts.SetCutGen4Rat();
    testLD2cuts.SetCutVz(Constants::RcutminVzLD2,Constants::RcutmaxVzLD2);     //vz cut for LD2 target
    testLD2cuts.SetCutGen4Rat();

    simCucuts.SetCutVz(Constants::RcutminVzCu,Constants::RcutmaxVzCu);     //vz cut for Cu target
    simCucuts.SetCutGen4Rat();
    testCucuts.SetCutVz(Constants::RcutminVzCu,Constants::RcutmaxVzCu);     //vz cut for Cu target
    testCucuts.SetCutGen4Rat();
    
    simSncuts.SetCutVz(Constants::RcutminVzSn,Constants::RcutmaxVzSn);     //vz cut for Sn targeta
    simSncuts.SetCutGen4Rat();
    testSncuts.SetCutVz(Constants::RcutminVzSndata,Constants::RcutmaxVzSndata);     //vz cut for Sn target
    testSncuts.SetCutGen4Rat();

    int sumevts = 0;
    Monunfold munfSimC2(simC2cuts, "C2_simf");
    Monunfold munfSimLD2(simLD2cuts, "LD2_simf");
    Monunfold munfSimCu(simCucuts, "Cu_simf");
    Monunfold munfSimSn(simSncuts, "Sn_simf");


    Monitoring monSimC1(simC1cuts, "C1_sim");     //This needs to be figured out ASAP
    Monitoring monSimC2(simC2cuts, "C2_sim");     //This needs to be figured out ASAP
    Monitoring monTestC1(testC1cuts, "C1_RGD");  //This needs to be figured out ASAP
    Monitoring monTestC2(testC2cuts, "C2_RGD");  //This needs to be figured out ASAP
    Monitoring monSimLD2(simLD2cuts, "LD2_sim");  
    Monitoring monTestLD2(testLD2cuts, "LD2_RGD");
    Monitoring monSimCu(simCucuts, "Cu_sim");
    Monitoring monTestCu(testCucuts, "Cu_RGD");
    Monitoring monSimSn(simSncuts, "Sn_sim");
    Monitoring monTestSn(testSncuts, "Sn_RGD");

    Ratio ratC1(testLD2cuts, testC1cuts, "C1_RGD");
    Ratio ratC2(testLD2cuts, testC2cuts, "C2_RGD");
    Ratio ratSn(testLD2cuts, testSncuts, "Sn_RGD");
    Ratio ratCu(testLD2cuts, testCucuts, "Cu_RGD");
    Ratio simratC1(simLD2cuts, simC1cuts, "C1_sim");
    Ratio simratC2(simLD2cuts, simC2cuts, "C2_sim");
    Ratio simratSn(simLD2cuts, simSncuts, "Sn_sim");
    Ratio simratCu(simLD2cuts, simCucuts, "Cu_sim");

    sratio sratC2(testLD2cuts, testC2cuts, "C2_RGD");
    sratio sratSn(testLD2cuts, testSncuts, "Sn_RGD");
    sratio sratCu(testLD2cuts, testCucuts, "Cu_RGD");
    sratio simsratC2(simLD2cuts, simC2cuts, "C2_sim");
    sratio simsratSn(simLD2cuts, simSncuts, "Sn_sim");
    sratio simsratCu(simLD2cuts, simCucuts, "Cu_sim");

    int totalevts = 19000;

    for (int i=0; i<totalevts; i++){
        simCxC  = Sim_CxC.ProcessEventsInFile();
        simCxC_MC  = Sim_CxC.ProcessEventsInFileMC();
        testCxC = RGD_CxC.ProcessEventsInFile(); 

        simSn = Sim_Sn.ProcessEventsInFile();
        simSn_MC = Sim_Sn.ProcessEventsInFileMC();
        testSn = RGD_CuSn.ProcessEventsInFile();

        simLD2 = Sim_LD2.ProcessEventsInFile();
        simLD2_MC = Sim_LD2.ProcessEventsInFileMC();
        testLD2 = RGD_LD2.ProcessEventsInFile();

//
//        simCu = Sim_Cu.ProcessEventsInFile();
//        simCu_MC = Sim_Cu.ProcessEventsInFileMC();
//        testCu = RGD_CuSn.ProcessEventsInFile();
//
//        simSn = Sim_Sn.ProcessEventsInFile();
//        simSn_MC = Sim_Sn.ProcessEventsInFileMC();
//        testSn = RGD_CuSn.ProcessEventsInFile();
        if (testLD2.has_value()) {
            Event eventtestLD2 = testLD2.value();
            eventtestLD2.SetTargetType(0);
            eventtestLD2.calcAll();
            monTestLD2.FillHistogramswCuts(eventtestLD2);
            ratC1.FillHistograms(eventtestLD2);
            ratC2.FillHistograms(eventtestLD2);
            ratSn.FillHistograms(eventtestLD2);
            sratCu.FillHistograms(eventtestLD2);
            sratC2.FillHistograms(eventtestLD2);

        }
        if (simLD2_MC.has_value()) {
            Event eventsimLD2_MC = simLD2_MC.value();
            eventsimLD2_MC.SetTargetType(0);
            eventsimLD2_MC.calcMCAll();
            munfSimLD2.FillHistogramswCutsMC(eventsimLD2_MC);
            if (simLD2.has_value()) {
                Event eventsimLD2 = simLD2.value();
                eventsimLD2.SetTargetType(0);
                eventsimLD2.calcAll();
                munfSimLD2.FillHistCompwCuts(eventsimLD2, eventsimLD2_MC);
                monSimLD2.FillHistogramswCuts(eventsimLD2);
                simratC1.FillHistograms(eventsimLD2);
                simratC2.FillHistograms(eventsimLD2);
                simratSn.FillHistograms(eventsimLD2);
                //simsratC2.FillHistograms(eventsimLD2);

            }
        }
        if (testCxC.has_value()) {
            Event eventtestCxC = testCxC.value();
            eventtestCxC.SetTargetType(1);
            eventtestCxC.calcAll();
            //munfTestC2.FillHistComp(eventtestCxC);
            monTestC2.FillHistogramswCuts(eventtestCxC);
            monTestC1.FillHistogramswCuts(eventtestCxC);
            ratC2.FillHistograms(eventtestCxC);
            ratC1.FillHistograms(eventtestCxC);
            //sratC2.FillHistograms(eventtestCxC);

        }
        if (simCxC_MC.has_value())                          {//CxC sim
            Event eventsimCxC_MC = simCxC_MC.value();
            eventsimCxC_MC.SetTargetType(1);
            eventsimCxC_MC.calcMCAll();
            munfSimC2.FillHistogramswCutsMC(eventsimCxC_MC);
                if ( simCxC.has_value()) {
                    Event eventsimCxC = simCxC.value();
                    eventsimCxC.SetTargetType(1);
                    eventsimCxC.calcAll();
                    munfSimC2.FillHistComp(eventsimCxC, eventsimCxC_MC);
                    monSimC1.FillHistogramswCuts(eventsimCxC);
                    monSimC2.FillHistogramswCuts(eventsimCxC);
                    //simratC1.FillHistograms(eventsimCxC);
                    simratC2.FillHistograms(eventsimCxC);
                    //simsratC2.FillHistograms(eventsimCxC);
                }
        }
        if (testSn.has_value()) {
            Event eventtestSn = testSn.value();
            eventtestSn.SetTargetType(1);
            eventtestSn.calcAll();
            monTestSn.FillHistogramswCuts(eventtestSn);
            ratSn.FillHistograms(eventtestSn);
            sratSn.FillHistograms(eventtestSn);

        }
        if (simSn_MC.has_value()) {
            Event eventsimSn_MC = simSn_MC.value();
            eventsimSn_MC.SetTargetType(1);
            eventsimSn_MC.calcMCAll();
            if (simSn.has_value()) {
                Event eventsimSn = simSn.value();
                eventsimSn.SetTargetType(1);
                eventsimSn.calcAll();
                munfSimSn.FillHistCompwCuts(eventsimSn, eventsimSn_MC);
                monSimSn.FillHistogramswCuts(eventsimSn);
                simratSn.FillHistograms(eventsimSn);
                simsratSn.FillHistograms(eventsimSn);
            }
        }
            //else{ counter_restCxC++;}
        files.displayProgress(i + 1, totalevts);
    }
    std::cout << "\nProcessing completed \n";
    std::cout << "//========= RGD data CxC ==========//  \n";
    monTestC2.SaveHistRoot("v6C2_test");
    monTestC1.SaveHistRoot("v6C1_test");
    std::cout << "//========= Simulation C2 ==========//  \n";
    monSimC2.SaveHistRoot("v6C2_sim");
    monSimC1.SaveHistRoot("v6C1_sim");
    monTestC2.DrawHelicityHistograms("HELIXc2");
    //ratC2.calcRatios();
    //ratC1.calcRatios();
    ratC2.calcR();  
    simratC2.calcR();
    ratSn.calcR();
    simratSn.calcR();

    ratC2.Rtargetsimcomp(simratC2);
    ratSn.Rtargetsimcomp(simratSn);
    simratSn.DrawSelfHistos(simratSn);
    simratSn.multiplotR();
    munfSimC2.DrawCompRECMC("MUNFv6C2");
    //munfSimC2.SaveHistRoot("munfC2_sim");

    //simsratC2.calcSratio();
    //simsratC2.calcAsymmetries();
    //simsratC2.writeAsymmToFile("testasymmetriesC2");
    //simsratC2.multiplotSratio();
    //simsratC2.DrawMonSinrat("augmonSratioC2");

    return 0;
}
