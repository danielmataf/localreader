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

    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/", filenamesLD2);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesCuSn);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesCxC);
    ////files.ParDir2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/", simufilesLD2);   //do not uncomment this line
    
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/Deut", simufilesLD2);  //uncomment for sim
    
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/C", simufilesCxC);  //uncomment for sim

    //Uncomment 4 test on ifarm, comment all above
    //files.Files2Vector("/volatile/clas12/rg-d/production/prod/v3_ob_LD2/dst/recon/019033/", filenamesLD2);
    //files.Files2Vector("/volatile/clas12/rg-d/production/prod/v3_ob_CuSn/dst/recon/019130/", filenamesCuSn);
    //files.Files2Vector("/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018835/", filenamesCxC);
    //uncommment also below for sim  on ifarm
    //files.SnDir2Vector("/volatile/clas12/dmat/gen/Deut/", simufilesLD2);  //uncomment for sim
    //files.SnDir2Vector("/volatile/clas12/dmat/gen/C/", simufilesCxC);  //uncomment for sim 


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
    CutSet simLD2cuts;  //LD2

    Ratio simrat3( simLD2cuts, simC1cuts, "C1_sim"); //calling the class with the corresponding cuts
    Ratio simrat4( simLD2cuts, simC2cuts, "C2_sim"); // RE calling the class with the corresponding cuts for study with Cu


    simC1cuts.SetCutVz(-3,-2);     //vz cut for C1 target
    simC2cuts.SetCutVz(-8,-7);     //vz cut for C2 target
    simLD2cuts.SetCutVz(-8,-2);    //vz cut for LD2 target
    //Fix target borders in constants. remember is not the definite vz cut, but a cut for target separation 



    int sumevts = 0;
    Monunfold munfTrueC1(truec1cuts, "C1_true");
    Monunfold munfSimC1(simC1cuts, "C1_simf");
    Monunfold munfSimC2(simC2cuts, "C2_simf");
    Monunfold munfsimLD2(simLD2cuts, "LD2_simf");
    Monitoring monSimC1(simC1cuts, "C1_sim");     //This needs to be figured out ASAP
    Monitoring monSimC2(simC2cuts, "C2_sim");     //This needs to be figured out ASAP

    //Ratio simratC1(simC1cuts, "C1_sim");


    int counter_el= 0.0;
    int counter_elLD2 = 0;
    int counter_elSn= 0.0;
    int counter_elCxC= 0.0;
    int counter_restCxC= 0.0;

    int counter_trueCxC= 0.0;
    int counter_elsimLD2= 0.0;
    int totalevts = 100;
    for (int i=0; i<totalevts; i++){
            simCxC  = Sim_CxC.ProcessEventsInFile();
            simCxC_MC  = Sim_CxC.ProcessEventsInFileMC();
            testCxC = MC_CxC.ProcessEventsInFile();
            simLD2 = Sim_LD2.ProcessEventsInFile();
            simLD2_MC = MC_LD2.ProcessEventsInFile();
            counter_trueCxC++;
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
                        counter_elCxC++;
                        Event eventsimCxC = simCxC.value();

                        eventsimCxC.SetTargetType(1);
                        eventsimCxC.calcAll();
                        //eventsimCxC.calcMCAll();
                        //eventsimCxC.Print();
                        //eventsimCxC.PrintMC();
                        //monSimC1.FillHistogramsNoCuts(eventsimCxC);
                        //monSimC2.FillHistogramsNoCuts(eventsimCxC);
                        //munfSimC1.Fill(eventsimCxC);    //this supposed to fil both rec and MC histos
                        //munfSim   C1.FillHistogramsNoCutsMC(eventsimCxC);
                        munfSimC1.FillHistogramsNoCuts(eventsimCxC);
                        monSimC1.FillHistogramswCuts(eventsimCxC);
                        munfSimC1.FillHistComp(eventsimCxC, eventsimCxC_MC);
                        simrat3.FillHistograms(eventsimCxC);
                        simrat4.FillHistograms(eventsimCxC);
                    }

            }

            if ( simLD2_MC.has_value()) {
                counter_elsimLD2++;
                Event eventsimLD2_MC = simLD2.value();
                eventsimLD2_MC.SetTargetType(1);
                eventsimLD2_MC.calcMCAll();
                munfsimLD2.FillHistogramsNoCutsMC(eventsimLD2_MC);
                if ( simLD2.has_value()) {
                    counter_elLD2++;
                    Event eventLD2 = testLD2.value();
                    eventLD2.SetTargetType(0);
                    eventLD2.calcAll();
                    munfsimLD2.FillHistogramsNoCuts(eventLD2);
                    munfsimLD2.FillHistComp(eventLD2, eventsimLD2_MC);
                    //simrat3.Fill(eventLD2);
                    //simrat4.Fill(eventLD2);
                }

            }
            else{ counter_restCxC++;}
            files.displayProgress(i + 1, totalevts);
            
    
    
    }
    std::cout << "\nProcessing completed \n";

    std::cout << "Events processed for simCtotal: " << counter_trueCxC << std::endl;

    std::cout << "Events processed for simCxC: " << counter_elCxC << std::endl;
    std::cout << "Events processed for simCrest: " << counter_restCxC << std::endl;
    std::cout << "is it compatible?  " << counter_restCxC +  counter_elCxC <<" =? " << counter_trueCxC<< std::endl;


    //munfSimC1.DrawHistograms("newC1monSIM_noCuts");
    munfSimC1.DrawHistoMC("newC1monSIM_noCutsMC");  
    munfSimC1.DrawHistoRec("newC1monSIM_noCutsREC");
    munfTrueC1.DrawHistoRec("newC1monTRUE_noCutsREC");
    munfSimC1.DrawCompRECMC("test2D");
    munfsimLD2.DrawHistoMC("test2D_LD2");
    simrat3.calcRcarbon(simrat4);
    simrat3.multiplotRbis();


    return 0;
}
