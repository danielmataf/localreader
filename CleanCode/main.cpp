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
    //files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesCuSn);
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
    //files.SnDir2Vector("/volatile/clas12/dmat/gen/LD2/", simufilesLD2);  //uncomment for sim
    //files.SnDir2Vector("/volatile/clas12/dmat/gen/CC/", simufilesCxC);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/gen/Cu/", simufilesCu);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/gen/Sn/", simufilesSn);  //uncomment for sim 


    std::cout<< "Hello world \n";
    EventReader Sim_CxC(simufilesCxC);
    EventReader RGD_CxC(filenamesCxC);
    EventReader Sim_LD2(simufilesLD2);
    EventReader RGD_LD2(filenamesLD2);   
    EventReader Sim_Cu(simufilesCu);   



    std::optional<Event> simCxC;
    std::optional<Event> simCxC_MC; //testing this if we can loop twice top read MC events//
    std::optional<Event> testCxC;
    std::optional<Event> simLD2;
    std::optional<Event> simLD2_MC;
    std::optional<Event> testLD2;
    std::optional<Event> simCu;
    std::optional<Event> simCu_MC;



    CutSet simC1cuts;   //C1
    CutSet simC2cuts;   //C2
    CutSet testC2cuts;  //C2 RGD  
    CutSet simLD2cuts;
    CutSet testLD2cuts;
    CutSet simCucuts;


    simC1cuts.SetCutVz(Constants::RcutminVzC1,Constants::RcutminVzC2);     //vz cut for C1 target
    simC2cuts.SetCutVz(Constants::RcutminVzC2, Constants::RcutmaxVzC2);     //vz cut for C2 target
    testC2cuts.SetCutVz(Constants::RcutminVzC2, Constants::RcutmaxVzC2);    //vz cut for C2 target
    simLD2cuts.SetCutGen4Rat();

    simLD2cuts.SetCutVz(Constants::RcutminVzLD2,Constants::RcutmaxVzLD2);     //vz cut for LD2 target
    testLD2cuts.SetCutGen4Rat();
    testLD2cuts.SetCutVz(Constants::RcutminVzLD2,Constants::RcutmaxVzLD2);     //vz cut for LD2 target

    simCucuts.SetCutVz(Constants::RcutminVzCu,Constants::RcutmaxVzCu);     //vz cut for Cu target
    simC1cuts.SetCutGen4Rat();
    simC2cuts.SetCutGen4Rat();
    testC2cuts.SetCutGen4Rat();
    
    simCucuts.SetCutGen4Rat();

    int sumevts = 0;
    Monunfold munfSimC2(simC2cuts, "C2_simf");
    Monunfold munfSimLD2(simLD2cuts, "LD2_simf");

    Monitoring monSimC1(simC1cuts, "C1_sim");     //This needs to be figured out ASAP
    Monitoring monSimC2(simC2cuts, "C2_sim");     //This needs to be figured out ASAP
    Monitoring monTestC2(testC2cuts, "C2_true");  //This needs to be figured out ASAP
    Monitoring monSimLD2(simLD2cuts, "LD2_sim");  
    Monitoring monTestLD2(testLD2cuts, "LD2_true");
    Monitoring monSimCu(simCucuts, "Cu_sim");

    Ratio simratC2(simLD2cuts, simC2cuts, "C2_sim");
    Ratio ratC2(testLD2cuts, testC2cuts, "C2_true");

    int totalevts = 9000;

    for (int i=0; i<totalevts; i++){
        simCxC  = Sim_CxC.ProcessEventsInFile();
        simCxC_MC  = Sim_CxC.ProcessEventsInFileMC();
        testCxC = RGD_CxC.ProcessEventsInFile(); 

        simLD2 = Sim_LD2.ProcessEventsInFile();
        simLD2_MC = Sim_LD2.ProcessEventsInFileMC();
        testLD2 = RGD_LD2.ProcessEventsInFile();

        //simCu = Sim_Cu.ProcessEventsInFile();
        //simCu_MC = Sim_Cu.ProcessEventsInFileMC();

        if (testCxC.has_value()) {
            Event eventtestCxC = testCxC.value();
            eventtestCxC.SetTargetType(1);
            eventtestCxC.calcAll();
            //munfTestC2.FillHistComp(eventtestCxC);
            monTestC2.FillHistogramswCuts(eventtestCxC);
        }
        if (simCxC_MC.has_value())                          {//CxC sim
            Event eventsimCxC_MC = simCxC_MC.value();
            eventsimCxC_MC.SetTargetType(1);
            eventsimCxC_MC.calcMCAll();
                if ( simCxC.has_value()) {
                    Event eventsimCxC = simCxC.value();
                    eventsimCxC.SetTargetType(1);
                    eventsimCxC.calcAll();
                    munfSimC2.FillHistCompwCuts(eventsimCxC, eventsimCxC_MC);
                    monSimC2.FillHistogramswCuts(eventsimCxC);

                }
        }
        if (testLD2.has_value()) {
            Event eventtestLD2 = testLD2.value();
            eventtestLD2.SetTargetType(0);
            eventtestLD2.calcAll();
            monTestLD2.FillHistogramswCuts(eventtestLD2);
        }
        if (simLD2_MC.has_value()) {
            Event eventsimLD2_MC = simLD2_MC.value();
            eventsimLD2_MC.SetTargetType(0);
            eventsimLD2_MC.calcMCAll();
            if (simLD2.has_value()) {
                Event eventsimLD2 = simLD2.value();
                eventsimLD2.SetTargetType(0);
                eventsimLD2.calcAll();
                munfSimLD2.FillHistCompwCuts(eventsimLD2, eventsimLD2_MC);
                monSimLD2.FillHistogramswCuts(eventsimLD2);
            }
        }
            //else{ counter_restCxC++;}
            files.displayProgress(i + 1, totalevts);
        
    
    
    }
    std::cout << "\nProcessing completed \n";

    monTestC2.SaveHistRoot("resetC2_test");
    monSimC2.SaveHistRoot("resetC2_sim");
    //monSimC1.SaveHistRoot("augC1_sim");
    //monSimC1.DrawHistograms("augC1_sim//");

    monTestLD2.SaveHistRoot("resetLD2_test");
    monSimLD2.SaveHistRoot("resetLD2_sim");

    simratC2.calcR();
    //ratC2.calcR();
    ratC2.DrawHistos(simratC2);
    simratC2.writeMatrixToFile("Rmatrix");
    ratC2.writeMatrixToFile("RmatrixRGD");
    simratC2.multiplotR();
    //munfSimC1.DrawCompRECMC("augcompC1_sim");
    munfSimC2.DrawCompRECMC("zwolfcompC2_sim");
    munfSimLD2.DrawCompRECMC("zwolfcompLD2_sim");
    //munfSimC2.DrawMomentainSim("augmomentumC2_sim");



    return 0;
}
