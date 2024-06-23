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



int main() {
    FilePathGenerator files;
    std::vector<std::string> filenamesLD2;
    std::vector<std::string> simufilesLD2;
    std::vector<std::string> simufilesSn;
    std::vector<std::string> filenamesCuSn;
    std::vector<std::string> filenamesCxC;
    std::vector<std::string> simufilesCxC;

    //files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/", filenamesLD2);
    //files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesCuSn);
    //files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesCxC);
    ////files.ParDir2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/", simufilesLD2);   //do not uncomment this line
    
    //files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/Deut", simufilesLD2);  //uncomment for sim
    
    //files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/C", simufilesCxC);  //uncomment for sim

    //Uncomment 4 test on ifarm, comment all above
    files.Files2Vector("/cache/hallb/scratch/rg-d/production/prod/v4ob_aideLD2/dst/recon/018528/", filenamesLD2);
    files.Files2Vector("/cache/hallb/scratch/rg-d/production/prod/v4ob_aideCuSn/dst/recon/018624/", filenamesCuSn);
    files.Files2Vector("/cache/hallb/scratch/rg-d/production/prod/v4ob_aideCxC/dst/recon/018451/", filenamesCxC);
    //uncommment also below for sim  on ifarm
    files.SnDir2Vector("/volatile/clas12/dmat/gen/C/", simufilesLD2);  //uncomment for sim
    files.SnDir2Vector("/volatile/clas12/dmat/gen/C/", simufilesCxC);  //uncomment for sim 


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
    std::optional<Event> simLD2_MC;
    //std::optional<Event> simCuSn;
    std::optional<Event> simCxC;
    std::optional<Event> simCxC_MC; //testing this if we can loop twice top read MC events//

    CutSet simC1cuts;   //C1
    CutSet simC2cuts;   //C2
    CutSet truec1cuts;  //C1
    CutSet truec2cuts;
    CutSet simLD2cuts;  //LD2
    CutSet trueLD2cuts;
    CutSet trueCucuts;
    CutSet trueSncuts;


    simC1cuts.SetCutVz(-3,-2);     //vz cut for C1 target
    simC2cuts.SetCutVz(-8,-7);     //vz cut for C2 target
    simLD2cuts.SetCutVz(-8,-2);    //vz cut for LD2 target
    //Fix target borders in constants. remember is not the definite vz cut, but a cut for target separation 
    trueLD2cuts.SetCutVz(-8,-2);
    trueSncuts.SetCutGen4Rat();
    trueSncuts.SetCutVz(-3.5,-1.5);
    trueCucuts.SetCutGen4Rat();
    trueCucuts.SetCutVz(-8.5,-6.5);
    truec1cuts.SetCutVz(-3.5,-1.5);
    truec2cuts.SetCutVz(-8.5,-6.5);


    Ratio rat( trueLD2cuts, trueSncuts, "Sn_data");
    Ratio rat2( trueLD2cuts, trueCucuts, "Cu_data");
    Ratio rat3( trueLD2cuts, truec1cuts, "C1_data");
    Ratio rat4( trueLD2cuts, truec2cuts, "C2_data");
    Ratio simrat3( simLD2cuts, simC1cuts, "C1_sim"); //calling the class with the corresponding cuts
    Ratio simrat4( simLD2cuts, simC2cuts, "C2_sim"); // RE calling the class with the corresponding cuts for study with Cu


    

    int sumevts = 0;
    Monunfold munfTrueC1(truec1cuts, "C1_true");
    Monunfold munfSimC1(simC1cuts, "C1_simf");
    Monunfold munfSimC2(simC2cuts, "C2_simf");
    Monunfold munfsimLD2(simLD2cuts, "LD2_simf");
    Monitoring monSimC1(simC1cuts, "C1_sim");     //This needs to be figured out ASAP
    Monitoring monSimC2(simC2cuts, "C2_sim");     //This needs to be figured out ASAP

    //Ratio simratC1(simC1cuts, "C1_sim");

    int totalevts = 2000000;

    for (int i=0; i<totalevts; i++){
            simCxC  = Sim_CxC.ProcessEventsInFile();
            simCxC_MC  = Sim_CxC.ProcessEventsInFileMC();
	    simLD2 = Sim_LD2.ProcessEventsInFile();
	    simLD2_MC = Sim_LD2.ProcessEventsInFileMC();
	    testLD2 = MC_LD2.ProcessEventsInFile();
	    testCuSn = MC_CuSn.ProcessEventsInFile();
	    testCxC = MC_CxC.ProcessEventsInFile(); 
	    if (testLD2.has_value()){
	        Event eventtestLD2 = testLD2.value();
	        eventtestLD2.SetTargetType(0);
	        eventtestLD2.calcAll();
		rat.FillHistograms(eventtestLD2);
		rat2.FillHistograms(eventtestLD2);
		rat3.FillHistograms(eventtestLD2);
	    }
	    if (testCuSn.has_value()) {
		Event eventtestCuSn = testCuSn.value();
                eventtestCuSn.SetTargetType(1);
                eventtestCuSn.calcAll();
                rat.FillHistograms(eventtestCuSn);
                rat2.FillHistograms(eventtestCuSn);
	    } 
            if (testCxC.has_value()) {
                Event eventtestCxC = testCxC.value();
                eventtestCxC.SetTargetType(1);              //is it targettype 1 or targettype 2? needs to be checked and modified accordingly
                eventtestCxC.calcAll();
                rat3.FillHistograms(eventtestCxC);
	    }
            if (simCxC_MC.has_value()) {//simulated data
                Event eventsimCxC_MC = simCxC_MC.value();
                eventsimCxC_MC.SetTargetType(1);
                eventsimCxC_MC.calcMCAll();
                munfSimC1.FillHistogramsNoCutsMC(eventsimCxC_MC);
                    // We can Potentially create a function in the same class that takes as  arguments both
                    // eventsimCxC and eventsimCxC_MC (filling histos with both banks)
                    // probably eventsimCxC can exist outside the condition ? 
                    // if not the condition simCxC.has_value can be contained in simCxC_MC.has_value 
                    // after all simCxC_MC.has_value is supposed to work always, and from there we can read REC bank
                    // ---- is working. But need to generalize for it to be used in true data too. And we need to check if is simulated by checking on MC rows
                    // ssuggestion! see notebook 
                    if ( simCxC.has_value()) {
                        Event eventsimCxC = simCxC.value();

                        eventsimCxC.SetTargetType(1);
                        eventsimCxC.calcAll();
                        munfSimC1.FillHistogramsNoCuts(eventsimCxC);
                        monSimC1.FillHistogramswCuts(eventsimCxC);
                        munfSimC1.FillHistComp(eventsimCxC, eventsimCxC_MC);
                        simrat3.FillHistograms(eventsimCxC);
                        simrat4.FillHistograms(eventsimCxC);
                    }

            }

            //else{ counter_restCxC++;}
            files.displayProgress(i + 1, totalevts);
            
    
    
    }
    std::cout << "\nProcessing completed \n";

   // std::cout << "Events processed for simCtotal: " << counter_trueCxC << std::endl;

   // std::cout << "Events processed for simCxC: " << counter_elCxC << std::endl;
   // std::cout << "Events processed for simCrest: " << counter_restCxC << std::endl;
   // std::cout << "is it compatible?  " << counter_restCxC +  counter_elCxC <<" =? " << counter_trueCxC<< std::endl;


    munfSimC1.DrawHistograms("newC1monSIM_noCuts");
    munfSimC1.DrawHistoMC("newC1monSIM_noCutsMC");  
    munfSimC1.DrawHistoRec("newC1monSIM_noCutsREC");
    munfTrueC1.DrawHistoRec("newC1monTRUE_noCutsREC");
    munfSimC1.DrawCompRECMC("test2D");
    munfsimLD2.DrawHistoMC("test2D_LD2");
    simrat3.calcRcarbon(simrat4);
    simrat3.multiplotRbis();
    rat.calcR();
    rat2.calcR();
    rat2.writeMatrixToFile("matrix2_test.txt");
    rat3.calcR();
    rat.multiplotR(rat2, rat3);
    return 0;
}
