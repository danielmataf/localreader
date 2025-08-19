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
    std::vector<std::string> simufilesLD2full;
    std::vector<std::string> simufilesLD2symm;

    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/0v6/", filenamesLD2);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/0v6/", filenamesCuSn);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/C_v0/0v6/", filenamesCxC);
    ////files.ParDir2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/", simufilesLD2);   //do not uncomment this line
    
    //files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/LD2", simufilesLD2);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novLD2", simufilesLD2);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novC", simufilesCxC);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novCu", simufilesCu);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novSn", simufilesSn);  //uncomment for sim

    //clean maybe //uncomment for test on ifarm
    //files.SnDir2Vector("/volatile/clas12/dmat/clean/LD2fullob/", simufilesLD2);  //uncomment for sim
    //files.SnDir2Vector("/volatile/clas12/dmat/clean/CCfullob/", simufilesCxC);  //uncomment for sim
    //files.SnDir2Vector("/volatile/clas12/dmat/clean/Cufullob/", simufilesCu);  //uncomment for sim
    //files.SnDir2Vector("/volatile/clas12/dmat/clean/Snfullob/", simufilesSn);  //uncomment for sim

    //check this for pass1
    //files.pass1search("/cache/clas12/rg-d/production/pass1/recon/LD2/", filenamesLD2);
    //files.pass1search("/cache/clas12/rg-d/production/pass1/recon/CuSn/", filenamesCuSn);
    //files.pass1search("/cache/clas12/rg-d/production/pass1/recon/CxC/", filenamesCxC);


    //Uncomment 4 test on ifarm, comment all above
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/skim_pass0v7/LD2/", filenamesLD2);
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/skim_pass0v7/CuSn/", filenamesCuSn);
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/skim_pass0v7/CxC/", filenamesCxC);
    ////uncommment also below for sim  on ifarm
    //files.SnDir2Vector("/volatile/clas12/dmat/test/fulltorus_aicv_newrgdLD2/", simufilesLD2);  //uncomment for sim
    //files.SnDir2Vector("/volatile/clas12/dmat/test/fullTorus_aicv_newrgdCC/", simufilesCxC);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/test/cv_newrgdCu/", simufilesCu);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/test/cv_newrgdSn/", simufilesSn);  //uncomment for sim 
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novLD2/", simufilesLD2full);  //uncomment for sim 
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novLD2", simufilesLD2symm);  //uncomment for sim 


    ////uncommment also below for sim BIS and/or THIRD on ifarm
    //files.SnDir2VectorBis("/volatile/clas12/dmat/test/fullTorus_aicv_newrgdCC/", simufilesCxC);  //uncomment for sim
    //files.SnDir2VectorThird("/volatile/clas12/dmat/test/fullTorus_aicv_newrgdCC/", simufilesCxC);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/test/cv_newrgdCu/", simufilesCu);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/test/cv_newrgdSn/", simufilesSn);  //uncomment for sim 


    std::cout<< "Analysis in progress... \n";
    EventReader RGD_CxC(filenamesCxC);
    EventReader RGD_LD2(filenamesLD2);   
    EventReader RGD_CuSn(filenamesCuSn);   
    EventReader Sim_CxC(simufilesCxC);
    EventReader Sim_LD2symm(simufilesLD2);
    EventReader Sim_LD2full(simufilesLD2full);
    EventReader Sim_Cu(simufilesCu);


    std::optional<Event> simCxC;
    std::optional<Event> simCxC_MC; //testing this if we can loop twice top read MC events//
    std::optional<Event> simCu;
    std::optional<Event> simCu_MC;
    std::optional<Event> testCxC;
    std::optional<Event> testLD2;
    std::optional<Event> testCu;
    std::optional<Event> testSn;
    std::optional<Event> testLD2symm;
    std::optional<Event> testLD2full;
    std::optional<Event> testLD2symmMC;
    std::optional<Event> testLD2fullMC;

    CutSet testC1cuts;  //C1 RGD
    CutSet testC2cuts;  //C2 RGD  
    CutSet simC2cuts;  
    CutSet simCucuts;  //Cu SIM

    CutSet testLD2cuts;
    CutSet testCucuts;
    CutSet testSncuts;



    //Defining cuts here. Vz cuts for SIM and RGD may differ -> check Mon Vz
    //simC1cuts.SetCutVz(Constants::RcutminVzC1sim,Constants::RcutminVzC1sim);     //vz cut for C1 target
    //testC1cuts.SetCutVz(Constants::RcutminVzC1data,Constants::RcutmaxVzC1data);     //previous cut, b4jan2025. We change for uncalib data (?) peaks shifted of -1cm
    testC1cuts.SetCutVz(Constants::v11cutminVzC1data,Constants::v11cutmaxVzC1data);     //Changing values in constants, also widening the cut
    testC1cuts.SetCutGen4Rat();

    testC2cuts.SetCutVz(Constants::v11cutminVzC2data, Constants::v11cutmaxVzC2data);    //vz cut for C2 target
    testC2cuts.SetCutGen4Rat();

    simC2cuts.SetCutVz(Constants::RcutminVzC2sim    , Constants::RcutmaxVzC2sim);     //vz cut for C2 target
    simC2cuts.SetCutGen4Rat();  
    simCucuts.SetCutVz(Constants::RcutminVzC2sim, Constants::RcutmaxVzC2sim);     //vz cut for Cu target
    simCucuts.SetCutGen4Rat();  //Cu target

    testLD2cuts.SetCutVz(Constants::v11cutminVzLD2data,Constants::v11cutmaxVzLD2data);     //vz cut for LD2 target
    testLD2cuts.SetCutGen4Rat();

    testCucuts.SetCutVz(Constants::v11cutminVzCudata,Constants::v11cutmaxVzCudata);     //vz cut for Cu target
    testCucuts.SetCutGen4Rat();
    
    testSncuts.SetCutVz(Constants::v11cutminVzSndata,Constants::v11cutmaxVzSndata);     //vz cut for Sn target
    testSncuts.SetCutGen4Rat();

    int sumevts = 0;

    Monitoring monSimC2(simC2cuts, "C2_sim");     //This needs to be figured out ASAP
    Monitoring monMCC2(simC2cuts, "C2_MC");
    Monitoring monSimCu(simCucuts, "Cu_sim");     //This needs to be figured out ASAP
    Monitoring monMCCu(simCucuts, "Cu_MC");
    Monunfold munftrueC2(simC2cuts, "C2_truef");    
    Monunfold munftrueCu(testCucuts, "Cu_truef"); 
    Monitoring monTestC1(testC1cuts, "C1_RGD");  //This needs to be figured out ASAP
    Monitoring monTestC2(testC2cuts, "C2_RGD");  //This needs to be figured out ASAP
    Monitoring monTestLD2(testLD2cuts, "LD2_RGD");
    Monunfold munfTestC2(testC2cuts, "C2_RGD");
    Monunfold munfTestCu(simCucuts, "C1_RGD");
    Monitoring monTestCu(testCucuts, "Cu_RGD");
    Monitoring monTestSn(testSncuts, "Sn_RGD");
    Monitoring monLD2symm(testLD2cuts, "LD2_RGD_symm");
    Monitoring monLD2full(testLD2cuts, "LD2_RGD_full");

    Ratio ratC1(testLD2cuts, testC1cuts, "C1_RGD");
    Ratio ratC2(testLD2cuts, testC2cuts, "C2_RGD");
    Ratio ratSn(testLD2cuts, testSncuts, "Sn_RGD");
    Ratio ratCu(testLD2cuts, testCucuts, "Cu_RGD");

    deltaptsq dptC1(testLD2cuts, testC1cuts, "C1_RGD");
    deltaptsq dptC2(testLD2cuts, testC2cuts, "C2_RGD");
    deltaptsq dptSn(testLD2cuts, testSncuts, "Sn_RGD");
    deltaptsq dptCu(testLD2cuts, testCucuts, "Cu_RGD");

    sratio sratC2(testLD2cuts, testC2cuts, "C2_RGD");
    sratio sratSn(testLD2cuts, testSncuts, "Sn_RGD");
    sratio sratCu(testLD2cuts, testCucuts, "Cu_RGD");

    cratio cratC1(testLD2cuts, testC1cuts, "C1_RGD");
    cratio cratC2(testLD2cuts, testC2cuts, "C2_RGD");
    cratio cratSn(testLD2cuts, testSncuts, "Sn_RGD");
    cratio cratCu(testLD2cuts, testCucuts, "Cu_RGD");


    int totalevts = 10000;
    // int totalevts =3500000;   //uncomment this on farm for fullT on C2 

    for (int i=0; i<totalevts; i++){
        testCxC = RGD_CxC.ProcessEventsInFile(); 
        testLD2 = RGD_LD2.ProcessEventsInFile();
        testSn = RGD_CuSn.ProcessEventsInFile();
        testCu = RGD_CuSn.ProcessEventsInFile();
        simCxC = Sim_CxC.ProcessEventsInFile();
        simCxC_MC = Sim_CxC.ProcessEventsInFileMC(); //testing this if we can loop twice top read MC events//
        simCu = Sim_Cu.ProcessEventsInFile();
        simCu_MC = Sim_Cu.ProcessEventsInFileMC(); //testing this if we can loop twice top read MC events//
        testLD2symm = Sim_LD2symm.ProcessEventsInFile();
        testLD2full = Sim_LD2full.ProcessEventsInFile();
        testLD2symmMC = Sim_LD2symm.ProcessEventsInFileMC();
        testLD2fullMC = Sim_LD2full.ProcessEventsInFileMC();

        if (testLD2.has_value()) {
            Event eventtestLD2 = testLD2.value();
            eventtestLD2.SetTargetType(0);
            eventtestLD2.calcAll();
            monTestLD2.FillHistogramswCuts(eventtestLD2);
            monTestLD2.CheckLargeBins(eventtestLD2); //check large bins in order to fill them with the correct values
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
            cratC2.FillHistograms(eventtestLD2);
            cratC1.FillHistograms(eventtestLD2);
            cratSn.FillHistograms(eventtestLD2);
            cratCu.FillHistograms(eventtestLD2);


        }

        if (testCxC.has_value()) {
            Event eventtestCxC = testCxC.value();
            eventtestCxC.SetTargetType(1);
            eventtestCxC.calcAll();

            //munfTestC2.FillHistComp(eventtestCxC);
            monTestC2.FillHistogramswCuts(eventtestCxC);
            monTestC2.thetabinning(eventtestCxC);    //for theta pending
            monTestC1.FillHistogramswCuts(eventtestCxC);
            munfTestC2.CheckLargeBins(eventtestCxC); //check large bins in order to fill them with the correct values
            ratC2.FillHistograms(eventtestCxC);
            ratC1.FillHistograms(eventtestCxC);
            dptC1.FillHistograms(eventtestCxC);
            dptC2.FillHistograms(eventtestCxC);
            sratC2.FillHistograms(eventtestCxC);
            cratC1.FillHistograms(eventtestCxC);
            cratC2.FillHistograms(eventtestCxC);
            cratC2.FillDebug(eventtestCxC);

        }
        if (simCxC_MC.has_value() )                          {//CxC sim
            Event eventsimCxC_MC = simCxC_MC.value();
            eventsimCxC_MC.SetTargetType(1);
            eventsimCxC_MC.calcMCAll();
            munftrueC2.FillHistogramswCutsMC(eventsimCxC_MC);

                if ( simCxC.has_value()) {
                    Event eventsimCxC = simCxC.value();
                    eventsimCxC.SetTargetType(1);
                    eventsimCxC.calcAll();
                    monSimC2.FillHistogramswCuts(eventsimCxC);
                    munftrueC2.FillHistogramswCuts(eventsimCxC);
                    munftrueC2.CheckLargeBins(eventsimCxC); //check large bins in order to fill them with the correct values

                    //munftrueC2.FillHistogramswCutsMC(eventsimCxC_MC);
                }
        }
        if (simCu_MC.has_value()) { //Cu sim
            Event eventsimCu_MC = simCu_MC.value();
            eventsimCu_MC.SetTargetType(1);
            eventsimCu_MC.calcMCAll();
            munftrueCu.FillHistogramswCutsMC(eventsimCu_MC);

                if (simCu.has_value()) {
                    Event eventsimCu = simCu.value();
                    eventsimCu.SetTargetType(1);
                    eventsimCu.calcAll();
                    monSimCu.FillHistogramswCuts(eventsimCu);
                    munftrueCu.FillHistogramswCuts(eventsimCu);
                    munftrueCu.CheckLargeBins(eventsimCu); //check large bins in order to fill them with the correct values
                }
        }
        if (testLD2symmMC.has_value()) { //LD2 sim
            Event eventtestLD2symmMC = testLD2symmMC.value();
            eventtestLD2symmMC.SetTargetType(0);
            eventtestLD2symmMC.calcMCAll();
            //monLD2symm.FillHistogramswCuts(eventtestLD2symmMC);
                if (testLD2symm.has_value()) {
                    Event eventtestLD2symm = testLD2symm.value();
                    eventtestLD2symm.SetTargetType(0);
                    eventtestLD2symm.calcAll();
                    monLD2symm.FillHistogramswCuts(eventtestLD2symm);
                }
        }
        if (testLD2fullMC.has_value()){
            Event eventtestLD2fullMC = testLD2fullMC.value();
            eventtestLD2fullMC.SetTargetType(0);
            eventtestLD2fullMC.calcMCAll();
                if (testLD2full.has_value()){
                    Event eventtestLD2full = testLD2full.value();
                    eventtestLD2full.SetTargetType(0);
                    eventtestLD2full.calcAll();
                    monLD2full.FillHistogramswCuts(eventtestLD2full);
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
            cratSn.FillHistograms(eventtestSn);

        }
        if (testCu.has_value()) {
            Event eventtestCu = testCu.value();
            eventtestCu.SetTargetType(1);
            eventtestCu.calcAll();
            monTestCu.FillHistogramswCuts(eventtestCu);
            ratCu.FillHistograms(eventtestCu);
            dptCu.FillHistograms(eventtestCu);
            sratCu.FillHistograms(eventtestCu);
            cratCu.FillHistograms(eventtestCu);
        }

            //else{ counter_restCxC++;}
        files.displayProgress(i + 1, totalevts);
    }
    munftrueC2.PrintRegionCounters(); //check large bins in order to fill them with the correct values
    monTestLD2.PrintRegionCounters(); //check large bins in order to fill them with the correct values
    std::cout << "\nProcessing completed \n";
    std::cout << "//========= RGD data CxC ==========//  \n";
    monTestC2.SaveHistRoot("janC2_test");
    monTestC2.DrawHistograms("monC2_test");
    monTestC2.SaveKeyHistograms();
    monTestSn.SaveHistRoot("janSn_test");
    monTestSn.DrawHistograms("SnUrgent");
    //monTestCu.SaveHistRoot("janCu_test");
    //monTestCu.DrawHistograms("CuUrgent");
    //monSimC2.SaveHistRoot("marC2_sim");
    //monSimC2.DrawHistograms("monC2_sim");
    cratC2.WriteDebugHistos("cosdebug.root");
    std::cout << "//========= Simulation C2 ==========//  \n";
    ratC2.saveRhistos();
    ratSn.saveRhistos();
    dptC2.saveDptHistos();
    dptSn.saveDptHistos();
    dptCu.saveDptHistos();
    cratC2.saveCratioHistos();
    cratSn.saveCratioHistos();
    cratCu.saveCratioHistos();
    monLD2full.SaveHistRoot("AprLD2full_test");
    monLD2symm.SaveHistRoot("AprLD2symm_test");
    monTestLD2.SaveHistRoot("testLD2_test");
    munftrueC2.SaveHistRoot("junC2_truef_test");
    munftrueC2.SaveHistMCRoot("junC2_truef_testMC");
    
    monLD2full.DrawEnergyHistograms("nrg_angLD2full");
    monLD2symm.DrawEnergyHistograms("nrg_angLD2symm");
    monTestLD2.DrawEnergyHistograms("nrg_angLD2RGDd");
    //ratCu.saveRhistos();
    //dptC2.saveDptHistos();
    //dptSn.saveDptHistos();
    //dptCu.saveDptHistos();
    //cratC2.saveCratioHistos();
    //cratSn.saveCratioHistos();
    //cratCu.saveCratioHistos();

    //ratC2.calcR();  
    //ratSn.calcR();
    //ratCu.calcR();
    //ratC2.calcR_xB_Q2_z();
    //dptC2.calcDpt();
    //dptSn.calcDpt();
    //dptCu.calcDpt();
    //cratC2.calcCratio();
    //cratSn.calcCratio();
    //cratCu.calcCratio();

    //sratC2.calcSratio();
    

//    ratC2.multiplotR(ratSn);
//    ratSn.multiplotR(ratCu, ratC2);
//    dptSn.multiplotDpt(dptCu, dptC2);
//    cratSn.multiplotCratio(cratCu, cratC2);
    //cratC2.calcCratioCarbon(cratC1);
    //cratC2.multiplotCratioBis();
    
    std::cout << "//========= Validations ==========//  \n";
//    ratC2.ValidateHistograms();
//    ratC2.LogBinContent();

    return 0;
}
