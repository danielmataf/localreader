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

    //files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/", filenamesLD2);
    //files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesCuSn);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/C_v0/018451/", filenamesCxC);
    ////files.ParDir2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/", simufilesLD2);   //do not uncomment this line
    
    //files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/Deut", simufilesLD2);  //uncomment for sim
    
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/C", simufilesCxC);  //uncomment for sim
    //files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/Cu", simufilesCu);  //uncomment for sim
    //files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/Sn", simufilesSn);  //uncomment for sim

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
   


    std::optional<Event> simCxC;
    std::optional<Event> simCxC_MC; //testing this if we can loop twice top read MC events//
    std::optional<Event> testCxC;


    CutSet simC1cuts;   //C1
    CutSet simC2cuts;   //C2
    CutSet testC2cuts;  //C2 RGD    


    simC1cuts.SetCutVz(Constants::RcutminVzC1,Constants::RcutminVzC2);     //vz cut for C1 target
    simC2cuts.SetCutVz(Constants::RcutminVzC2, Constants::RcutmaxVzC2);     //vz cut for C2 target
    testC2cuts.SetCutVz(Constants::RcutminVzC2, Constants::RcutmaxVzC2);    //vz cut for C2 target
    simC1cuts.SetCutGen4Rat();
    simC2cuts.SetCutGen4Rat();
    testC2cuts.SetCutGen4Rat();

    int sumevts = 0;
    Monunfold munfSimC1(simC1cuts, "C1_simf");
    Monunfold munfSimC2(simC2cuts, "C2_simf");

    Monitoring monSimC1(simC1cuts, "C1_sim");     //This needs to be figured out ASAP
    Monitoring monSimC2(simC2cuts, "C2_sim");     //This needs to be figured out ASAP
    Monitoring monTestC2(testC2cuts, "C2_true");  //This needs to be figured out ASAP

    int totalevts = 40000;

    for (int i=0; i<totalevts; i++){
        //simCxC  = Sim_CxC.ProcessEventsInFile();
        //simCxC_MC  = Sim_CxC.ProcessEventsInFileMC();
        testCxC = RGD_CxC.ProcessEventsInFile(); 
        if (testCxC.has_value()) {
            Event eventtestCxC = testCxC.value();
            eventtestCxC.SetTargetType(1);
            eventtestCxC.calcAll();
            //munfTestC2.FillHistComp(eventtestCxC);
            monTestC2.FillHistogramswCuts(eventtestCxC);
        }
        //if (simCxC_MC.has_value())                          {//CxC sim
        //    Event eventsimCxC_MC = simCxC_MC.value();
        //    eventsimCxC_MC.SetTargetType(1);
        //    eventsimCxC_MC.calcMCAll();
        //    //eventsimCxC_MC.ReadRunconfig();
        //    //std::cout << "---------------- " << std::endl;
        //    //std::cout << "Event count prt MC: " << i << std::endl;
        //    //std::cout << "Event count fct MC: " << eventsimCxC_MC.GetEvtnbr() << std::endl;
        //        if ( simCxC.has_value()) {
        //            Event eventsimCxC = simCxC.value();
        //            eventsimCxC.SetTargetType(1);
        //            eventsimCxC.calcAll();
        //            munfSimC1.FillHistCompwCuts(eventsimCxC, eventsimCxC_MC);
        //            munfSimC2.FillHistComp(eventsimCxC, eventsimCxC_MC);
        //            monSimC2.FillHistogramswCuts(eventsimCxC);
        //            //std::cout << "Event count prt SIM: " << i << std::endl;
        //            //std::cout << "Event count fct SIM: " << eventsimCxC.GetEvtnbr() << std::endl;
//
        //        }
        //}

            //else{ counter_restCxC++;}
            files.displayProgress(i + 1, totalevts);
            
    
    
    }
    std::cout << "\nProcessing completed \n";

    monTestC2.SaveHistRoot("augC2_test");
    //monSimC1.SaveHistRoot("augC1_sim");
    //monSimC1.DrawHistograms("augC1_sim//");

    monSimC2.SaveHistRoot("augC2_sim");

    //munfSimC1.DrawCompRECMC("augcompC1_sim");
    //munfSimC2.DrawCompRECMC("augcompC2_sim");
    //munfSimC2.DrawMomentainSim("augmomentumC2_sim");



    return 0;
}
