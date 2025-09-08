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

    

    //check this for pass1
    files.ParDir2Vector("/cache/clas12/rg-d/production/pass1/recon/LD2/dst/recon/", filenamesLD2);
    files.ParDir2Vector("/cache/clas12/rg-d/production/pass1/recon/CuSn/dst/recon/", filenamesCuSn);
    files.ParDir2Vector("/cache/clas12/rg-d/production/pass1/recon/CxC/dst/recon/", filenamesCxC);


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


    std::optional<Event> testCxC;
    std::optional<Event> testLD2;
    std::optional<Event> testCu;
    std::optional<Event> testSn;
    std::optional<Event> testLD2full;
    std::optional<Event> testLD2fullMC;

    CutSet testC1cuts;  //C1 RGD
    CutSet testC2cuts;  //C2 RGD  

    CutSet testLD2cuts;
    CutSet testCucuts;
    CutSet testSncuts;



    //Defining cuts here. Vz cuts for SIM and RGD may differ -> check Mon Vz
    //simC1cuts.SetCutVz(Constants::RcutminVzC1sim,Constants::RcutminVzC1sim);     //vz cut for C1 target
    //testC1cuts.SetCutVz(Constants::RcutminVzC1data,Constants::RcutmaxVzC1data);     //previous cut, b4jan2025. We change for uncalib data (?) peaks shifted of -1cm
    testC1cuts.SetCutVz(Constants::v11cutminVzC1data,Constants::v11cutmaxVzC1data);     //Changing values in constants, also widening the cut
    testC1cuts.SetCutGen4Rat();

    //testC2cuts.SetCutVz(Constants::v11cutminVzC2data, Constants::v11cutmaxVzC2data);    //vz cut for C2 target
    testC2cuts.SetCutGen4Rat();
    testC2cuts.SetCutVz(-10.56, 5);    //vz cut for C2 target


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
    Monitoring monTestLD2(testLD2cuts, "LD2_RGD");
    Monitoring monTestCu(testCucuts, "Cu_RGD");
    Monitoring monTestSn(testSncuts, "Sn_RGD");
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


//    int totalevts = 10;
    int totalevts =10000000;   //uncomment this on farm for fullT on C2 

    for (int i=0; i<totalevts; i++){
        testCxC = RGD_CxC.ProcessEventsInFile(); 
        testLD2 = RGD_LD2.ProcessEventsInFile();
        testSn = RGD_CuSn.ProcessEventsInFile();
        testCu = RGD_CuSn.ProcessEventsInFile();

        if (testLD2.has_value()) {
            Event eventtestLD2 = testLD2.value();
            eventtestLD2.SetTargetType(0);
            eventtestLD2.calcAll();
            monTestLD2.FillHistogramswCuts(eventtestLD2);
            monTestLD2.CheckLargeBins(eventtestLD2); //check large bins in order to fill them with the correct values
	    monTestLD2.CheckFewBins(eventtestLD2);    
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
	    monTestSn.CheckFewBins(eventtestSn);
            monTestSn.CheckLargeBins(eventtestSn); //check large bins in order to fill them with the correct values

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
            monTestCu.CheckLargeBins(eventtestCu); //check large bins in order to fill them with the correct values
            monTestCu.CheckFewBins(eventtestCu); //check large bins in order to fill them with the correct values

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
