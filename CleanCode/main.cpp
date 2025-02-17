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
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/skim_pass0v7/LD2/", filenamesLD2);
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/skim_pass0v7/CuSn/", filenamesCuSn);
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/skim_pass0v7/CxC/", filenamesCxC);
    ////uncommment also below for sim  on ifarm
    //files.SnDir2Vector("/volatile/clas12/dmat/test/fulltorus_aicv_newrgdLD2/", simufilesLD2);  //uncomment for sim
    //files.SnDir2Vector("/volatile/clas12/dmat/test/fullTorus_aicv_newrgdCC/", simufilesCxC);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/test/cv_newrgdCu/", simufilesCu);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/test/cv_newrgdSn/", simufilesSn);  //uncomment for sim 

    ////uncommment also below for sim BIS and/or THIRD on ifarm
    //files.SnDir2VectorBis("/volatile/clas12/dmat/test/fullTorus_aicv_newrgdCC/", simufilesCxC);  //uncomment for sim
    //files.SnDir2VectorThird("/volatile/clas12/dmat/test/fullTorus_aicv_newrgdCC/", simufilesCxC);  //uncomment for sim 
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
    //testC1cuts.SetCutVz(Constants::RcutminVzC1data,Constants::RcutmaxVzC1data);     //previous cut, b4jan2025. We change for uncalib data (?) peaks shifted of -1cm
    testC1cuts.SetCutVz(Constants::RcutminVzC1data,Constants::RcutmaxVzC1data);     //Changing values in constants, also widening the cut
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

    deltaptsq dptC1(testLD2cuts, testC1cuts, "C1_RGD");
    deltaptsq dptC2(testLD2cuts, testC2cuts, "C2_RGD");
    deltaptsq dptSn(testLD2cuts, testSncuts, "Sn_RGD");
    deltaptsq dptCu(testLD2cuts, testCucuts, "Cu_RGD");
    deltaptsq simdptC1(simLD2cuts, simC1cuts, "C1_sim");
    deltaptsq simdptC2(simLD2cuts, simC2cuts, "C2_sim");
    deltaptsq simdptSn(simLD2cuts, simSncuts, "Sn_sim");
    deltaptsq simdptCu(simLD2cuts, simCucuts, "Cu_sim");

    sratio sratC2(testLD2cuts, testC2cuts, "C2_RGD");
    sratio sratSn(testLD2cuts, testSncuts, "Sn_RGD");
    sratio sratCu(testLD2cuts, testCucuts, "Cu_RGD");
    sratio simsratC2(simLD2cuts, simC2cuts, "C2_sim");
    sratio simsratSn(simLD2cuts, simSncuts, "Sn_sim");
    sratio simsratCu(simLD2cuts, simCucuts, "Cu_sim");

    int totalevts = 19000;
    // int totalevts =3500000;   //uncomment this on farm for fullT on C2 

    for (int i=0; i<totalevts; i++){
        simCxC  = Sim_CxC.ProcessEventsInFile();
        simCxC_MC  = Sim_CxC.ProcessEventsInFileMC();
        testCxC = RGD_CxC.ProcessEventsInFile(); 

        simLD2 = Sim_LD2.ProcessEventsInFile();
        simLD2_MC = Sim_LD2.ProcessEventsInFileMC();
        testLD2 = RGD_LD2.ProcessEventsInFile();

        testSn = RGD_CuSn.ProcessEventsInFile();
        testCu = RGD_CuSn.ProcessEventsInFile();

        if (testLD2.has_value()) {
            Event eventtestLD2 = testLD2.value();
            eventtestLD2.SetTargetType(0);
            eventtestLD2.calcAll();
            monTestLD2.FillHistogramswCuts(eventtestLD2);
            ratC1.FillHistograms(eventtestLD2);
            ratC2.FillHistograms(eventtestLD2);
            ratSn.FillHistograms(eventtestLD2);
            ratCu.FillHistograms(eventtestLD2);
            dptC1.FillHistograms(eventtestLD2);
            dptC2.FillHistograms(eventtestLD2);
            dptSn.FillHistograms(eventtestLD2);
            dptCu.FillHistograms(eventtestLD2);
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
                munfSimLD2.FillHistComp(eventsimLD2, eventsimLD2_MC);
                monSimLD2.FillHistogramswCuts(eventsimLD2);
                simratC1.FillHistograms(eventsimLD2);
                simratC2.FillHistograms(eventsimLD2);
                simdptC1.FillHistograms(eventsimLD2);
                simdptC2.FillHistograms(eventsimLD2);
                simdptCu.FillHistograms(eventsimLD2);
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
            dptC1.FillHistograms(eventtestCxC);
            dptC2.FillHistograms(eventtestCxC);
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
                    munfSimC2.FillHistComp(eventsimCxC, eventsimCxC_MC);
                    monSimC1.FillHistogramswCuts(eventsimCxC);
                    monSimC2.FillHistogramswCuts(eventsimCxC);
                    simratC1.FillHistograms(eventsimCxC);
                    simratC2.FillHistograms(eventsimCxC);
                    simsratC2.FillHistograms(eventsimCxC);
                    simdptC1.FillHistograms(eventsimCxC);
                    simdptC2.FillHistograms(eventsimCxC);
                }
        }
        if (testSn.has_value()) {
            Event eventtestSn = testSn.value();
            eventtestSn.SetTargetType(1);
            eventtestSn.calcAll();

            monTestSn.FillHistogramswCuts(eventtestSn);
            //monTestSn.FillHistogramsNoCuts(eventtestSn);
            ratSn.FillHistograms(eventtestSn);
            dptSn.FillHistograms(eventtestSn);
            sratSn.FillHistograms(eventtestSn);
        }
        if (testCu.has_value()) {
            Event eventtestCu = testCu.value();
            eventtestCu.SetTargetType(1);
            eventtestCu.calcAll();
            monTestCu.FillHistogramswCuts(eventtestCu);
            ratCu.FillHistograms(eventtestCu);
            dptCu.FillHistograms(eventtestCu);
            sratCu.FillHistograms(eventtestCu);
        }
            //else{ counter_restCxC++;}
        files.displayProgress(i + 1, totalevts);
    }
    std::cout << "\nProcessing completed \n";
    std::cout << "//========= RGD data CxC ==========//  \n";
    monTestC2.SaveHistRoot("janC2_test");
    monTestC2.DrawHistograms("monC2_test");
    monTestC2.SaveKeyHistograms();
    munfSimC2.DrawCompRECMC("DeltasC2_sim");
    monTestSn.SaveHistRoot("janSn_test");
    monTestSn.DrawHistograms("SnUrgent");
    monTestCu.SaveHistRoot("janCu_test");
    monTestCu.DrawHistograms("CuUrgent");
    std::cout << "//========= Simulation C2 ==========//  \n";
    monSimC2.SaveHistRoot("janC2_sim");
    monSimC1.SaveHistRoot("janC1_sim");
    //monTestLD2.SaveHistRoot("janLD2_test");
    //monSimLD2.SaveHistRoot("janLD2_sim");
    //monTestC2.DrawHelicityHistograms("HELIXc2");
    //ratC2.calcRatios();
    simratC1.DrawHistos(ratC1);
    simratC1.DrawSelfHistos(ratC1);

    ratC2.calcR();  
    ratSn.calcR();
    ratCu.calcR();
    dptC2.calcDpt();
    dptSn.calcDpt();
    dptCu.calcDpt();
    simdptC2.calcDpt();
    sratC2.calcSratio();
    simsratC2.calcSratio();

    ratC2.multiplotR(ratSn);
    ratSn.multiplotR(ratCu, ratC2);
    dptSn.multiplotDpt(dptCu, dptC2);
    
    //ratC2.Rtargetsimcomp(simratC2);
    //dptC2.Dpttargetsimcomp(simdptC2);
    //simratC2.DrawHistos(simratC1);
    //simratC2.DrawSelfHistos(simratC1);
    std::cout << "//========= Validations ==========//  \n";
    ratC2.ValidateHistograms();
    ratC2.LogBinContent();
    //std::cout << "//===================//  \n";
//
    //dptC2.ValidateHistograms();
    //dptC2.LogBinContent();
    //std::cout << "//===================//  \n";
//
    //sratC2.ValidateHistograms();
    //sratC2.LogBinContent();




    return 0;
}
