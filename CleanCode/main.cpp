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

    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/", filenamesLD2);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesCuSn);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/C_v0/018451/", filenamesCxC);
    ////files.ParDir2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/", simufilesLD2);   //do not uncomment this line
    
    //files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/LD2", simufilesLD2);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/Deut", simufilesLD2);  //uncomment for sim
    
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/C", simufilesCxC);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/C", simufilesCu);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/Sn", simufilesSn);  //uncomment for sim

    //Uncomment 4 test on ifarm, comment all above
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/prod/v4ob_aideLD2/dst/recon/018528/", filenamesLD2);
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/prod/v4ob_aideCuSn/dst/recon/018624/", filenamesCuSn);
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/prod/v4ob_aideCxC/dst/recon/018451/", filenamesCxC);
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




    CutSet simC1cuts;   //C1
    CutSet simC2cuts;   //C2
    CutSet testC2cuts;  //C2 RGD  
    CutSet simLD2cuts;
    CutSet testLD2cuts;
    CutSet simCucuts;
    CutSet testCucuts;
    CutSet simSncuts;
    CutSet testSncuts;



    //Defining cuts here. Vz cuts for SIM and RGD may differ -> check Mon Vz
    simC1cuts.SetCutVz(Constants::RcutminVzC1,Constants::RcutminVzC1);     //vz cut for C1 target
    simC1cuts.SetCutGen4Rat();

    simC2cuts.SetCutVz(Constants::RcutminVzC2, Constants::RcutmaxVzC2);     //vz cut for C2 target
    simC2cuts.SetCutGen4Rat();
    testC2cuts.SetCutVz(Constants::RcutminVzC2, Constants::RcutmaxVzC2);    //vz cut for C2 target
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
    testSncuts.SetCutVz(Constants::RcutminVzSn,Constants::RcutmaxVzSn);     //vz cut for Sn target
    testSncuts.SetCutGen4Rat();

    int sumevts = 0;
    Monunfold munfSimC2(simC2cuts, "C2_simf");
    Monunfold munfSimLD2(simLD2cuts, "LD2_simf");
    Monunfold munfSimCu(simCucuts, "Cu_simf");
    Monunfold munfSimSn(simSncuts, "Sn_simf");


    Monitoring monSimC1(simC1cuts, "C1_sim");     //This needs to be figured out ASAP
    Monitoring monSimC2(simC2cuts, "C2_sim");     //This needs to be figured out ASAP
    Monitoring monTestC2(testC2cuts, "C2_RGD");  //This needs to be figured out ASAP
    Monitoring monSimLD2(simLD2cuts, "LD2_sim");  
    Monitoring monTestLD2(testLD2cuts, "LD2_RGD");
    Monitoring monSimCu(simCucuts, "Cu_sim");
    Monitoring monTestCu(testCucuts, "Cu_RGD");
    Monitoring monSimSn(simSncuts, "Sn_sim");
    Monitoring monTestSn(testSncuts, "Sn_RGD");


    Ratio ratC2(testLD2cuts, testC2cuts, "C2_RGD");
    Ratio ratSn(testLD2cuts, testSncuts, "Sn_RGD");
    Ratio ratCu(testLD2cuts, testCucuts, "Cu_RGD");
    Ratio simratC2(simLD2cuts, simC2cuts, "C2_sim");
    Ratio simratSn(simLD2cuts, simSncuts, "Sn_sim");
    Ratio simratCu(simLD2cuts, simCucuts, "Cu_sim");

    sratio sratC2(testLD2cuts, testC2cuts, "C2_RGD");
    sratio sratSn(testLD2cuts, testSncuts, "Sn_RGD");
    sratio sratCu(testLD2cuts, testCucuts, "Cu_RGD");
    sratio simsratC2(simLD2cuts, simC2cuts, "C2_sim");
    sratio simsratSn(simLD2cuts, simSncuts, "Sn_sim");
    sratio simsratCu(simLD2cuts, simCucuts, "Cu_sim");

    int totalevts = 9000;

    for (int i=0; i<totalevts; i++){
        simCxC  = Sim_CxC.ProcessEventsInFile();
        simCxC_MC  = Sim_CxC.ProcessEventsInFileMC();
        testCxC = RGD_CxC.ProcessEventsInFile(); 

        simLD2 = Sim_LD2.ProcessEventsInFile();
        simLD2_MC = Sim_LD2.ProcessEventsInFileMC();
        testLD2 = RGD_LD2.ProcessEventsInFile();

        simCu = Sim_Cu.ProcessEventsInFile();
        simCu_MC = Sim_Cu.ProcessEventsInFileMC();
        testCu = RGD_CuSn.ProcessEventsInFile();

        simSn = Sim_Sn.ProcessEventsInFile();
        simSn_MC = Sim_Sn.ProcessEventsInFileMC();
        testSn = RGD_CuSn.ProcessEventsInFile();

        if (testLD2.has_value()) {
            Event eventtestLD2 = testLD2.value();
            eventtestLD2.SetTargetType(0);
            eventtestLD2.calcAll();
            monTestLD2.FillHistogramswCuts(eventtestLD2);
            ratC2.FillHistograms(eventtestLD2);
            ratSn.FillHistograms(eventtestLD2);
            ratCu.FillHistograms(eventtestLD2);
            sratC2.FillHistograms(eventtestLD2);
            sratSn.FillHistograms(eventtestLD2);
            sratCu.FillHistograms(eventtestLD2);
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
                simratC2.FillHistograms(eventsimLD2);
                simratSn.FillHistograms(eventsimLD2);
                simratCu.FillHistograms(eventsimLD2);
                simsratC2.FillHistograms(eventsimLD2);
                simsratSn.FillHistograms(eventsimLD2);
                simsratCu.FillHistograms(eventsimLD2);

            }
        }


        if (testCu.has_value()) {
            Event eventtestCu = testCu.value();
            eventtestCu.SetTargetType(1);
            eventtestCu.calcAll();
            monTestCu.FillHistogramswCuts(eventtestCu);
            ratCu.FillHistograms(eventtestCu);
            sratCu.FillHistograms(eventtestCu);

        }
        if (simCu_MC.has_value()) {
            Event eventsimCu_MC = simCu_MC.value();
            eventsimCu_MC.SetTargetType(1);
            eventsimCu_MC.calcMCAll();
            if (simCu.has_value()) {
                Event eventsimCu = simCu.value();
                eventsimCu.SetTargetType(1);
                eventsimCu.calcAll();
                munfSimCu.FillHistCompwCuts(eventsimCu, eventsimCu_MC);
                monSimCu.FillHistogramswCuts(eventsimCu);
                simratCu.FillHistograms(eventsimCu);
                simsratCu.FillHistograms(eventsimCu);
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
        if (testCxC.has_value()) {
            Event eventtestCxC = testCxC.value();
            eventtestCxC.SetTargetType(1);
            eventtestCxC.calcAll();
            //munfTestC2.FillHistComp(eventtestCxC);
            monTestC2.FillHistogramsNoCuts(eventtestCxC);
            ratC2.FillHistograms(eventtestCxC);
            sratC2.FillHistograms(eventtestCxC);

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
                    munfSimC2.FillHistCompwCuts(eventsimCxC, eventsimCxC_MC);
                    monSimC2.FillHistogramsNoCuts(eventsimCxC);
                    simratC2.FillHistograms(eventsimCxC);
                    simsratC2.FillHistograms(eventsimCxC);
                }
        }
            //else{ counter_restCxC++;}
        files.displayProgress(i + 1, totalevts);
    }
    std::cout << "\nProcessing completed \n";
    std::cout << "//========= RGD data C2 ==========//  \n";
    monTestC2.SaveHistRoot("novC2_test");
    std::cout << "//========= Simulation C2 ==========//  \n";
    monSimC2.SaveHistRoot("novC2_sim");
    std::cout << "//========= RGD data LD2 ==========//  \n";
    monTestLD2.SaveHistRoot("novLD2_test");
    std::cout << "//========= Simulation LD2 ==========//  \n";
    monSimLD2.SaveHistRoot("novLD2_sim");
    std::cout << "//========= RGD data Sn ==========//  \n";
    monTestSn.SaveHistRoot("novSn_test");
    std::cout << "//========= Simulation Sn ==========//  \n";
    monSimSn.SaveHistRoot("novSn_sim");
    std::cout << "//========= RGD data Cu ==========//  \n";
    monTestCu.SaveHistRoot("novCu_test");
    std::cout << "//========= Simulation Cu ==========//  \n";
    monSimCu.SaveHistRoot("novCu_sim");

    
    //monTestCu.SaveHistRoot("octCu_test");
    //monSimCu.SaveHistRoot("octCu_sim");

    //ratC2.DrawHistos(simratC2);
    //simratC2.writeMatrixToFile("Rmatrix");
    //ratC2.writeMatrixToFile("RmatrixRGD");
    //simratC2.multiplotR();
    //munfSimC1.DrawCompRECMC("augcompC1_sim");
    //////!!!ratC2.calcR();
    //////!!!ratSn.calcR();
    //////!!!ratCu.calcR();
    //simratC2.calcR();
    //simratSn.calcR();
    //simratCu.calcR();
    sratC2.calcSratio();
    sratSn.calcSratio();
    sratCu.calcSratio();
    simsratC2.calcSratio();
    simsratSn.calcSratio();
    simsratCu.calcSratio();
    //////!!!sratC2.calcSratio1D();
    //////!!!sratC2.writeMatrixToFile("sRmatrixC2");
    //////!!!sratSn.writeMatrixToFile("sRmatrixSn");
    sratC2.multiplotSratio();
    sratC2.calcAsymmetries();
    sratC2.DrawMonSinrat("augmonSratioC2");
    //sratSn.multiplotR();    

    //munfSimC2.DrawCompRECMC("seventcompC2_sim");
    //munfSimLD2.DrawCompRECMC("seventcompLD2_sim");
    //munfSimCu.DrawCompRECMC("seventcompCu_sim");
    //munfSimSn.DrawCompRECMC("seventcompSn_sim");
    ////munfSimC2.DrawMomentainSim("augmomentumC2_sim");



    return 0;
}
