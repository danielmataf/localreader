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
    std::vector<std::string> filenamesCuSn;
    std::vector<std::string> filenamesCxC;

    std::vector<std::string> simufilesLD2;
    std::vector<std::string> simufilesCu;
    std::vector<std::string> simufilesSn;
    std::vector<std::string> simufilesCxC;    



    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/0v6/", filenamesLD2);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/0v6/", filenamesCuSn);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/C_v0/0v6/", filenamesCxC);
    ////files.ParDir2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/", simufilesLD2);   //do not uncomment this line
    //files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/LD2", simufilesLD2);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novLD2", simufilesLD2);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novC", simufilesCxC);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novCu", simufilesCu);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/novSn", simufilesSn);

    //check this for pass1
    //files.ParDir2Vector("/cache/clas12/rg-d/production/pass1/recon/LD2/dst/recon/", filenamesLD2);
    //files.ParDir2Vector("/cache/clas12/rg-d/production/pass1/recon/CuSn/dst/recon/", filenamesCuSn);
    //files.ParDir2Vector("/cache/clas12/rg-d/production/pass1/recon/CxC/dst/recon/", filenamesCxC);


    //Uncomment 4 test on ifarm, comment all above
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/skim_pass0v7/LD2/", filenamesLD2);
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/skim_pass0v7/CuSn/", filenamesCuSn);
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/skim_pass0v7/CxC/", filenamesCxC);
    ////uncommment also below for sim  on ifarm
    //files.SnDir2Vector("/volatile/clas12/dmat/test/fulltorus_aicv_newrgdLD2/", simufilesLD2);  //uncomment for sim
 //   files.SnDir2Vector("/volatile/clas12/dmat/clean/Cufullob", simufilesCu);  //uncomment for sim 
    //files.pass1search("/cache/clas12/rg-d/production/pass1/recon/CxC/dst/recon/", simufilesCu);  //uncomment for sim 
//    files.SnDir2Vector("/volatile/clas12/dmat/clean/Cufullob", simufilesCxC);  //uncomment for sim 

//files.SnDir2Vector("/volatile/clas12/dmat/test/cv_newrgdSn/", simufilesSn);  //uncomment for sim 


    ////uncommment also below for sim BIS and/or THIRD on ifarm
    //files.SnDir2VectorBis("/volatile/clas12/dmat/test/fullTorus_aicv_newrgdCC/", simufilesCxC);  //uncomment for sim
    //files.SnDir2VectorThird("/volatile/clas12/dmat/test/fullTorus_aicv_newrgdCC/", simufilesCxC);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/test/cv_newrgdCu/", simufilesCu);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/test/cv_newrgdSn/", simufilesSn);  //uncomment for sim 


    std::cout<< "Analysis in progress... \n";
    EventReader RGD_CxC(filenamesCxC);
    EventReader RGD_LD2(filenamesLD2);   
    EventReader RGD_CuSn(filenamesCuSn);   

    EventReader SIM_LD2(simufilesLD2);
    EventReader SIM_Cu(simufilesCu);
    EventReader SIM_Sn(simufilesSn);
    EventReader SIM_CxC(simufilesCxC);

    std::optional<Event> testCxC;
    std::optional<Event> testLD2;
    std::optional<Event> testCu;
    std::optional<Event> testSn;
    std::optional<Event> testLD2full;
    std::optional<Event> testLD2fullMC;
    
    std::optional<Event> simuLD2;
    std::optional<Event> simuCu;
    std::optional<Event> simuSn;
    std::optional<Event> simuCxC;

    std::optional<Event> simuLD2_MC;
    std::optional<Event> simuCu_MC;
    std::optional<Event> simuSn_MC;
    std::optional<Event> simuCxC_MC;

    CutSet testC1cuts;  //C1 RGD
    CutSet testC2cuts;  //C2 RGD  
    CutSet testCCcuts;  //CxC RGD
    CutSet testLD2cuts;
    CutSet testCucuts;
    CutSet testSncuts;
    //We'll use the same cuts for simu since they seem wide enough.
    


    //Defining cuts here. Vz cuts for SIM and RGD may differ -> check Mon Vz
    //simC1cuts.SetCutVz(Constants::RcutminVzC1sim,Constants::RcutminVzC1sim);     //vz cut for C1 target
    //testC1cuts.SetCutVz(Constants::RcutminVzC1data,Constants::RcutmaxVzC1data);     //previous cut, b4jan2025. We change for uncalib data (?) peaks shifted of -1cm

    //testC2cuts.SetCutVz(Constants::v11cutminVzC2data, Constants::v11cutmaxVzC2data);    //vz cut for C2 target
    testCCcuts.SetCutGen4Rat();
    testCCcuts.SetCutVz(-10.56, 5);    //vz cut for CxC target, using C2 values for now
    testC1cuts.SetCutGen4Rat();
    testC1cuts.SetCutVz(-10.55,-6.62);     //vz cut for Cu target
    testC2cuts.SetCutGen4Rat();
    testC2cuts.SetCutVz(-6.21,-1.26);


    //testLD2cuts.SetCutVz(Constants::v11cutminVzLD2data,Constants::v11cutmaxVzLD2data);     //vz cut for LD2 target
    testLD2cuts.SetCutGen4Rat();
    testLD2cuts.SetCutVz(-15,5);
    //testCucuts.SetCutVz(Constants::v11cutminVzCudata,Constants::v11cutmaxVzCudata);     //vz cut for Cu target
    testCucuts.SetCutGen4Rat();
    testCucuts.SetCutVz(-10.55,-6.62);     //vz cut for Cu target

    //testSncuts.SetCutVz(Constants::v11cutminVzSndata,Constants::v11cutmaxVzSndata);     //vz cut for Sn target
    testSncuts.SetCutGen4Rat();
    testSncuts.SetCutVz(-6.21,-1.26);

    int sumevts = 0;

    Monitoring monTestC1(testC1cuts, "C1_RGD");  //This needs to be figured out ASAP
    Monitoring monTestC2(testC2cuts, "C2_RGD");  //This needs to be figured out ASAP
    Monitoring monTestCC(testCCcuts, "CxC_RGD");
    Monitoring monTestLD2(testLD2cuts, "LD2_RGD");
    Monitoring monTestCu(testCucuts, "Cu_RGD");
    Monitoring monTestSn(testSncuts, "Sn_RGD");
    Monitoring monLD2full(testLD2cuts, "LD2_RGD_full");

    Monitoring monSimLD2(testLD2cuts, "LD2_SIMU");
    Monitoring monSimCu(testCucuts, "Cu_SIMU");
    Monitoring monSimSn(testSncuts, "Sn_SIMU");
    Monitoring monSimCxC(testCCcuts, "CxC_SIMU");
    Monitoring monSimC2(testC2cuts, "C2_SIMU");
    Monitoring monSimC1(testC1cuts, "C1_SIMU");

    Monunfold munfSimLD2(testLD2cuts, "LD2_SIMU");
    Monunfold munfSimCu(testCucuts, "Cu_SIMU");
    Monunfold munfSimSn(testSncuts, "Sn_SIMU");
    Monunfold munfSimCxC(testCCcuts, "CxC_SIMU");

    Ratio ratC1(testLD2cuts, testC1cuts, "C1_RGD");
    Ratio ratC2(testLD2cuts, testC2cuts, "C2_RGD");
    Ratio ratCC(testLD2cuts, testCCcuts, "CxC_RGD");
    Ratio ratSn(testLD2cuts, testSncuts, "Sn_RGD");
    Ratio ratCu(testLD2cuts, testCucuts, "Cu_RGD");


    deltaptsq dptC1(testLD2cuts, testC1cuts, "C1_RGD");
    deltaptsq dptC2(testLD2cuts, testC2cuts, "C2_RGD");
    deltaptsq dptCC(testLD2cuts, testCCcuts, "CxC_RGD");
    deltaptsq dptSn(testLD2cuts, testSncuts, "Sn_RGD");
    deltaptsq dptCu(testLD2cuts, testCucuts, "Cu_RGD");

    sratio sratC2(testLD2cuts, testC2cuts, "C2_RGD");
    sratio sratSn(testLD2cuts, testSncuts, "Sn_RGD");
    sratio sratCu(testLD2cuts, testCucuts, "Cu_RGD");
    sratio sratCC(testLD2cuts, testCCcuts, "CxC_RGD");

    cratio cratC1(testLD2cuts, testC1cuts, "C1_RGD");
    cratio cratC2(testLD2cuts, testC2cuts, "C2_RGD");
    cratio cratCC(testLD2cuts, testCCcuts, "CxC_RGD");
    cratio cratSn(testLD2cuts, testSncuts, "Sn_RGD");
    cratio cratCu(testLD2cuts, testCucuts, "Cu_RGD");


    int totalevts = 10000;       //local
//    int totalevts =10000000;   //farm

    for (int i=0; i<totalevts; i++){
        testCxC = RGD_CxC.ProcessEventsInFile(); 
        testLD2 = RGD_LD2.ProcessEventsInFile();
        testSn = RGD_CuSn.ProcessEventsInFile();
        testCu = RGD_CuSn.ProcessEventsInFile();

        simuCxC = SIM_CxC.ProcessEventsInFile();
        simuLD2 = SIM_LD2.ProcessEventsInFile();
        simuSn = SIM_Sn.ProcessEventsInFile();
        simuCu = SIM_Cu.ProcessEventsInFile();
        simuCxC_MC = SIM_CxC.ProcessEventsInFileMC();
        simuLD2_MC = SIM_LD2.ProcessEventsInFileMC();
        simuSn_MC = SIM_Sn.ProcessEventsInFileMC();
        simuCu_MC = SIM_Cu.ProcessEventsInFileMC();

        if (testLD2.has_value()) {
            Event eventtestLD2 = testLD2.value();
            eventtestLD2.SetTargetType(0);
            eventtestLD2.calcAll();
            monTestLD2.FillHistogramswCuts( eventtestLD2);
            monTestLD2.CheckLargeBins(eventtestLD2); //check large bins in order to fill them with the correct values
	        monTestLD2.CheckFewBins(eventtestLD2);   
            monTestLD2.FillPrePostCutHistograms(eventtestLD2); 
            
            ratC1.FillHistograms(eventtestLD2);
            ratC2.FillHistograms(eventtestLD2);
            ratCC.FillHistograms(eventtestLD2);
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
            monTestCC.FillPrePostCutHistograms(eventtestCxC);
            monTestC2.thetabinning(eventtestCxC);    //for theta pending
            monTestC1.FillHistogramswCuts(eventtestCxC);
            monTestC2.CheckLargeBins(eventtestCxC); //check large bins in order to fill them with the correct values
            monTestC2.CheckFewBins(eventtestCxC); //check large bins in order to fill them with the correct values

	    ratC2.FillHistograms(eventtestCxC);
            ratC1.FillHistograms(eventtestCxC);
            dptC1.FillHistograms(eventtestCxC);
            dptC2.FillHistograms(eventtestCxC);
            sratC2.FillHistograms(eventtestCxC);
            cratC1.FillHistograms(eventtestCxC);
            cratC2.FillHistograms(eventtestCxC);
            cratC2.FillDebug(eventtestCxC);

        }

        if (testSn.has_value()) {
            Event eventtestSn = testSn.value();
            eventtestSn.SetTargetType(1);
            eventtestSn.calcAll();

            monTestSn.FillHistogramswCuts(eventtestSn);
            monTestSn.FillPrePostCutHistograms(eventtestSn);
            //monTestSn.FillHistogramsNoCuts(eventtestSn);
            ratSn.FillHistograms(eventtestSn);
            dptSn.FillHistograms(eventtestSn);
            sratSn.FillHistograms(eventtestSn);
            cratSn.FillHistograms(eventtestSn);
	        monTestSn.CheckFewBins(eventtestSn);
            monTestSn.CheckLargeBins(eventtestSn); //check large bins in order to fill them with the correct values

        }
        if (testCu.has_value()) {
            Event eventtestCu = testCu.value();
            eventtestCu.SetTargetType(1);
            eventtestCu.calcAll();
            monTestCu.FillHistogramswCuts(eventtestCu);
            monTestCu.FillPrePostCutHistograms(eventtestCu);
            ratCu.FillHistograms(eventtestCu);
            dptCu.FillHistograms(eventtestCu);
            sratCu.FillHistograms(eventtestCu);
            cratCu.FillHistograms(eventtestCu);
            monTestCu.CheckLargeBins(eventtestCu); //check large bins in order to fill them with the correct values
            monTestCu.CheckFewBins(eventtestCu); //check large bins in order to fill them with the correct values

        }
        //now the simus 
        if (simuLD2_MC.has_value()){
            Event eventsimuLD2_MC = simuLD2_MC.value();
            eventsimuLD2_MC.SetTargetType(0);
            eventsimuLD2_MC.calcMCAll();
            if (simuLD2.has_value()){
                Event eventsimuLD2 = simuLD2.value();
                eventsimuLD2.SetTargetType(0);
                eventsimuLD2.calcAll();
                monSimLD2.FillHistogramswCuts(eventsimuLD2);
            }
        }
        if (simuCxC_MC.has_value()){
            Event eventsimuCxC_MC = simuCxC_MC.value();
            eventsimuCxC_MC.SetTargetType(1);
            eventsimuCxC_MC.calcMCAll();
            if (simuCxC.has_value()){
                Event eventsimuCxC = simuCxC.value();
                eventsimuCxC.SetTargetType(1);
                eventsimuCxC.calcAll();
                monSimCxC.FillHistogramswCuts(eventsimuCxC);
            }
        }
        if (simuSn_MC.has_value()){
            Event eventsimuSn_MC = simuSn_MC.value();
            eventsimuSn_MC.SetTargetType(1);
            eventsimuSn_MC.calcMCAll();
            if (simuSn.has_value()){
                Event eventsimuSn = simuSn.value();
                eventsimuSn.SetTargetType(1);
                eventsimuSn.calcAll();
                monSimSn.FillHistogramswCuts(eventsimuSn);
            }
        }
        if (simuCu_MC.has_value()){
            Event eventsimuCu_MC = simuCu_MC.value();
            eventsimuCu_MC.SetTargetType(1);
            eventsimuCu_MC.calcMCAll();
            if (simuCu.has_value()){
                Event eventsimuCu = simuCu.value();
                eventsimuCu.SetTargetType(1);
                eventsimuCu.calcAll();
                monSimCu.FillHistogramswCuts(eventsimuCu);
            }
        }


            //else{ counter_restCxC++;}
        files.displayProgress(i + 1, totalevts);
    }
//    monTestLD2.CalcElectronRatio();
//
//  
    std::cout << "region large counters LD2 \n";
    monTestLD2.PrintRegionCounters(); //check large bins in order to fill them with the correct values
    std::cout << "region large counters Sn \n";
    monTestSn.PrintRegionCounters(); //check large bins in order to fill them with the correct values
    std::cout << "region few counters LD2 \n";
    monTestLD2.PrintFewRegionCounters(); //check large bins in order to fill them with the correct values
    std::cout << "region few counters Sn \n";
    monTestSn.PrintFewRegionCounters(); //check large bins in order to fill them with the correct values
    std::cout << " EMC electrons Sn  \n";
    monTestSn.CalcElectronRatio(monTestLD2);
    monTestSn.CalcLARGEElectronRatio(monTestLD2);
    std::cout << " EMC electrons Cu  \n";
    monTestCu.CalcElectronRatio(monTestLD2);
    monTestCu.CalcLARGEElectronRatio(monTestLD2);
    std::cout << " EMC electrons CxC  \n";
    monTestC2.CalcElectronRatio(monTestLD2);
    monTestC2.CalcLARGEElectronRatio(monTestLD2);


std::cout << "\nProcessing completed \n";
    std::cout << "//========= RGD data CxC ==========//  \n";
    monTestC2.SaveHistRoot("janC2_test");
    monTestC2.DrawHistograms("monC2_test");
    monTestC2.SaveKeyHistograms();
    monTestSn.SaveHistRoot("pass1Sn");
    monTestC2.SaveHistRoot("pass1C2");
    monTestLD2.SaveHistRoot("pass1LD2");
    monTestCu.SaveHistRoot("pass1Cu");

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
    monTestCC.SaveFINRoot("fin_CxC_test");
    monTestSn.SaveFINRoot("fin_Sn_test");
    monTestCu.SaveFINRoot("fin_Cu_test");
    monTestLD2.SaveFINRoot("fin_LD2_test");
    monLD2full.SaveHistRoot("AprLD2full_test");
    monTestLD2.SaveHistRoot("testLD2_test");
//    munftrueC2.SaveHistRoot("junC2_truef_test");
  //  munftrueC2.SaveHistMCRoot("junC2_truef_testMC");
    
    monLD2full.DrawEnergyHistograms("nrg_angLD2full");
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
