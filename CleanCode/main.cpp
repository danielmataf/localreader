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
    std::vector<std::string> simufilesCu;
    std::vector<std::string> filenamesCuSn;
    std::vector<std::string> filenamesCxC;
    std::vector<std::string> simufilesCxC;

    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/018428/", filenamesLD2);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/CuSn_v0/018348/", filenamesCuSn);
    files.Files2Vector("/home/matamoros/Desktop/LumiScanDta/C_v0/018451/", filenamesCxC);
    ////files.ParDir2Vector("/home/matamoros/Desktop/LumiScanDta/LD2_v0/", simufilesLD2);   //do not uncomment this line
    
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/Deut", simufilesLD2);  //uncomment for sim
    
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/C", simufilesCxC);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/Cu", simufilesCu);  //uncomment for sim
    files.SnDir2Vector("/home/matamoros/Desktop/LumiScanDta/simtestfolder/Sn", simufilesSn);  //uncomment for sim

    //Uncomment 4 test on ifarm, comment all above
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/prod/v4ob_aideLD2/dst/recon/018528/", filenamesLD2);
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/prod/v4ob_aideCuSn/dst/recon/018624/", filenamesCuSn);
    //files.Files2Vector("/cache/hallb/scratch/rg-d/production/prod/v4ob_aideCxC/dst/recon/018451/", filenamesCxC);
    ////uncommment also below for sim  on ifarm
    //files.SnDir2Vector("/volatile/clas12/dmat/genfiles/LD2/", simufilesLD2);  //uncomment for sim
    //files.SnDir2Vector("/volatile/clas12/dmat/genfiles/C/", simufilesCxC);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/genfiles/Cu/", simufilesCu);  //uncomment for sim 
    //files.SnDir2Vector("/volatile/clas12/dmat/genfiles/Sn/", simufilesSn);  //uncomment for sim 


    std::cout<< "Hello world \n";
    //EventReader MC(filenames);
    EventReader MC_LD2(filenamesLD2);
    EventReader MC_CuSn(filenamesCuSn);
    EventReader MC_CxC(filenamesCxC);
    EventReader Sim_LD2(simufilesLD2);
    EventReader Sim_Cu(simufilesSn);
    EventReader Sim_Sn(simufilesCu);
    EventReader Sim_CxC(simufilesCxC);
   


    //std::optional<Event> test;
    std::optional<Event> testLD2;
    std::optional<Event> testCuSn;
    std::optional<Event> testCxC;
    std::optional<Event> simLD2;
    std::optional<Event> simLD2_MC;
    std::optional<Event> simCu;
    std::optional<Event> simCu_MC;
    std::optional<Event> simSn;
    std::optional<Event> simSn_MC;
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
    CutSet simCucuts;   //Cu
    CutSet simSncuts;   //Sn


    simC1cuts.SetCutVz(-3,-2);     //vz cut for C1 target
    simC2cuts.SetCutVz(-8,-7);     //vz cut for C2 target
    simLD2cuts.SetCutVz(-8,-2);    //vz cut for LD2 target
    simCucuts.SetCutVz(-8.5,-6.5); //vz cut for Cu target
    simSncuts.SetCutVz(-3.5,-1.5); //vz cut for Sn target
    simC1cuts.SetCutGen4Rat();
    simC2cuts.SetCutGen4Rat();
    simLD2cuts.SetCutGen4Rat();
    simCucuts.SetCutGen4Rat();
    simSncuts.SetCutGen4Rat();
    //Fix target borders in constants. remember is not the definite vz cut, but a cut for target separation 
    trueLD2cuts.SetCutVz(-8,-2);
    trueLD2cuts.SetCutGen4Rat();
    trueSncuts.SetCutGen4Rat();
    trueSncuts.SetCutVz(-3.5,-1.5);
    trueCucuts.SetCutGen4Rat();
    trueCucuts.SetCutVz(-8.5,-6.5);
    truec1cuts.SetCutVz(-3.5,-1.5);
    truec1cuts.SetCutGen4Rat();
    truec2cuts.SetCutVz(-8.5,-6.5);
    truec2cuts.SetCutGen4Rat();

    int sumevts = 0;
    Monunfold munfTrueC1(truec1cuts, "C1_truef");
    Monunfold munfTrueC2(truec2cuts, "C2_truef");
    Monunfold munfTrueLD2(trueLD2cuts, "LD2_truef");
    Monunfold munfTrueCu(trueCucuts, "Cu_truef");
    Monunfold munfTrueSn(trueSncuts, "Sn_truef");
    Monunfold munfSimC1(simC1cuts, "C1_simf");
    Monunfold munfSimC2(simC2cuts, "C2_simf");
    Monunfold munfsimLD2(simLD2cuts, "LD2_simf");
    //Monunfold munfSimCu(simCucuts, "Cu_simf");
    Monunfold munfSimSn(simSncuts, "Sn_simf");
    Monitoring monSimC1(simC1cuts, "C1_sim");     //This needs to be figured out ASAP
    Monitoring monSimC2(simC2cuts, "C2_sim");     //This needs to be figured out ASAP
    Monitoring monSimLD2(simLD2cuts, "LD2_sim");
    Monitoring monSimCu(simCucuts, "Cu_sim");
    Monitoring monSimSn(simSncuts, "Sn_sim");
    Monitoring monTrueC1(truec1cuts, "C1_true");
    Monitoring monTrueC2(truec2cuts, "C2_true");


    Ratio rat( trueLD2cuts, trueSncuts, "Sn_data");
    Ratio rat2( trueLD2cuts, trueCucuts, "Cu_data");
    Ratio rat3( trueLD2cuts, truec1cuts, "C1_data");
    Ratio rat4( trueLD2cuts, truec2cuts, "C2_data");
    Ratio simrat( simLD2cuts, simSncuts, "Cu_sim");
    Ratio simrat2( simLD2cuts, simCucuts, "Sn_sim");
    Ratio simrat3( simLD2cuts, simC1cuts, "C1_sim"); //calling the class with the corresponding cuts
    Ratio simrat4( simLD2cuts, simC2cuts, "C2_sim"); // RE calling the class with the corresponding cuts for study with Cu

    cratio simcratC1(simLD2cuts, simC1cuts, "C1_sim");
    cratio simcratC2(simLD2cuts, simC2cuts, "C2_sim");
    cratio simcratCu(simLD2cuts, simCucuts, "Cu_sim");
    cratio simcratSn(simLD2cuts, simSncuts, "Sn_sim");
    cratio truecratC1(trueLD2cuts, truec1cuts, "C1_true");
    cratio truecratC2(trueLD2cuts, truec2cuts, "C2_true");
    cratio truecratCu(trueLD2cuts, trueCucuts, "Cu_true");
    cratio truecratSn(trueLD2cuts, trueSncuts, "Sn_true");     

    deltaptsq broadsimC1(simLD2cuts, simC1cuts, "C1_sim");
    deltaptsq broadsimC2(simLD2cuts, simC2cuts, "C2_sim");
    deltaptsq broadtrueC1(trueLD2cuts, truec1cuts, "C1_true");
    deltaptsq broadtrueC2(trueLD2cuts, truec2cuts, "C2_true");
    deltaptsq broadtrueCu(trueLD2cuts, trueCucuts, "Cu_true");
    deltaptsq broadtrueSn(trueLD2cuts, trueSncuts, "Sn_true");
    deltaptsq broadsimCu(simLD2cuts, simCucuts, "Cu_sim");
    deltaptsq broadsimSn(simLD2cuts, simSncuts, "Sn_sim");     

    


    //Ratio simratC1(simC1cuts, "C1_sim");

    int totalevts = 5000;

    for (int i=0; i<totalevts; i++){
        simCxC  = Sim_CxC.ProcessEventsInFile();
        simCxC_MC  = Sim_CxC.ProcessEventsInFileMC();
	    simLD2 = Sim_LD2.ProcessEventsInFile();
	    simLD2_MC = Sim_LD2.ProcessEventsInFileMC();
        simCu = Sim_Cu.ProcessEventsInFile();
        simCu_MC = Sim_Cu.ProcessEventsInFileMC();
        simSn = Sim_Sn.ProcessEventsInFile();
        simSn_MC = Sim_Sn.ProcessEventsInFileMC();
	    testLD2 = MC_LD2.ProcessEventsInFile();
	    testCuSn = MC_CuSn.ProcessEventsInFile();
	    testCxC = MC_CxC.ProcessEventsInFile(); 
	    if (testLD2.has_value()){
	        Event eventtestLD2 = testLD2.value();
	        eventtestLD2.SetTargetType(0);
	        eventtestLD2.calcAll();
            munfTrueLD2.FillHistogramswCuts(eventtestLD2);
		    rat.FillHistograms(eventtestLD2);
		    rat2.FillHistograms(eventtestLD2);
		    rat3.FillHistograms(eventtestLD2);
            truecratSn.FillHistograms(eventtestLD2);
		    truecratCu.FillHistograms(eventtestLD2);
		    truecratC1.FillHistograms(eventtestLD2);
		    truecratC2.FillHistograms(eventtestLD2);
            broadtrueC1.FillHistograms(eventtestLD2);
            broadtrueC2.FillHistograms(eventtestLD2);
            broadtrueCu.FillHistograms(eventtestLD2);
            broadtrueSn.FillHistograms(eventtestLD2);
            broadtrueSn.FillOnlyptandz(eventtestLD2);


	    }
	    if (testCuSn.has_value()) {
		    Event eventtestCuSn = testCuSn.value();
            eventtestCuSn.SetTargetType(1);
            eventtestCuSn.calcAll();
            munfTrueCu.FillHistogramswCuts(eventtestCuSn);
            munfTrueSn.FillHistogramswCuts(eventtestCuSn);
            rat.FillHistograms(eventtestCuSn);
            rat2.FillHistograms(eventtestCuSn);
	        truecratSn.FillHistograms(eventtestCuSn);
            truecratCu.FillHistograms(eventtestCuSn);
            broadtrueCu.FillHistograms(eventtestCuSn);
            broadtrueSn.FillHistograms(eventtestCuSn);
            broadtrueSn.FillOnlyptandz(eventtestCuSn);

        } 
        if (testCxC.has_value()) {
            Event eventtestCxC = testCxC.value();
            eventtestCxC.SetTargetType(1);              //is it targettype 1 or targettype 2? needs to be checked and modified accordingly
            eventtestCxC.calcAll();
            munfTrueC1.FillHistogramswCuts(eventtestCxC);
            munfTrueC2.FillHistogramswCuts(eventtestCxC);
            monTrueC1.FillHistogramswCuts(eventtestCxC);
            monTrueC2.FillHistogramswCuts(eventtestCxC);
            rat3.FillHistograms(eventtestCxC);
            rat4.FillHistograms(eventtestCxC);
            truecratC1.FillHistograms(eventtestCxC);
            truecratC2.FillHistograms(eventtestCxC);
            broadtrueC1.FillHistograms(eventtestCxC);
            broadtrueC2.FillHistograms(eventtestCxC);


	    }
        if (simLD2_MC.has_value()){
            Event eventsimLD2_MC = simLD2_MC.value();
            eventsimLD2_MC.SetTargetType(0);
            eventsimLD2_MC.calcMCAll();
            if (simLD2.has_value()){
                Event eventsimLD2 = simLD2.value();
                eventsimLD2.SetTargetType(0);
                eventsimLD2.calcAll();
                munfsimLD2.FillHistogramswCuts(eventsimLD2);
                monSimLD2.FillHistogramsNoCuts(eventsimLD2);
                simrat.FillHistograms(eventsimLD2);
                simrat2.FillHistograms(eventsimLD2);
                simrat3.FillHistograms(eventsimLD2);
                simrat4.FillHistograms(eventsimLD2);
            }
        }
        if (simCu_MC.has_value()){
            Event eventsimCu_MC = simCu_MC.value();
            eventsimCu_MC.SetTargetType(1);
            eventsimCu_MC.calcMCAll();
            if (simCu.has_value()){
                Event eventsimCu = simCu.value();
                eventsimCu.SetTargetType(1);
                eventsimCu.calcAll();
                munfSimSn.FillHistogramswCuts(eventsimCu);
                monSimCu.FillHistogramsNoCuts(eventsimCu);
                simrat2.FillHistograms(eventsimCu);
                simcratCu.FillHistograms(eventsimCu);
                broadsimCu.FillHistograms(eventsimCu);

            }
        }
        if (simSn_MC.has_value()){
            Event eventsimSn_MC = simSn_MC.value();
            eventsimSn_MC.SetTargetType(1);
            eventsimSn_MC.calcMCAll();
            if (simSn.has_value()){
                Event eventsimSn = simSn.value();
                eventsimSn.SetTargetType(1);
                eventsimSn.calcAll();
                //eventsimSn.Print();
                //munfSimSn.FillHistogramswCuts(eventsimSn);
                monSimSn.FillHistogramsNoCuts(eventsimSn);
                simrat.FillHistograms(eventsimSn);
                simcratSn.FillHistograms(eventsimSn);
                broadsimSn.FillHistograms(eventsimSn);
            }
        }
        if (simCxC_MC.has_value()) {//simulated data
            Event eventsimCxC_MC = simCxC_MC.value();
            eventsimCxC_MC.SetTargetType(1);
            eventsimCxC_MC.calcMCAll();
                if ( simCxC.has_value()) {
                    Event eventsimCxC = simCxC.value();
                    eventsimCxC.SetTargetType(1);
                    eventsimCxC.calcAll();
                    munfSimC1.FillHistogramswCuts(eventsimCxC);
                    monSimC1.FillHistogramsNoCuts(eventsimCxC);
                    monSimC2.FillHistogramsNoCuts(eventsimCxC);
                    simrat3.FillHistograms(eventsimCxC);
                    simrat4.FillHistograms(eventsimCxC);
                    simcratC1.FillHistograms(eventsimCxC);
                    simcratC2.FillHistograms(eventsimCxC);
                    broadsimC1.FillHistograms(eventsimCxC);
                    broadsimC2.FillHistograms(eventsimCxC);


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


    monSimC1.DrawHistograms("MONC1SIM");
    monSimC2.DrawHistograms("MONC2SIM");
    monSimLD2.DrawHistograms("MONLD2SIM");
    monSimCu.DrawHistograms("MONCuSIM");
    monSimSn.DrawHistograms("MONSnSIM");
    simrat3.calcRcarbon(simrat4);
    rat3.calcRcarbon(rat4);
    //simrat3.multiplotRbis();
    //rat3.multiplotRbis();
    rat.calcR();
    rat2.calcR();
    //rat2.writeMatrixToFile("matrix2_test.txt");
    rat3.calcR();
    //rat.multiplotR(rat2, rat3);
    truecratC1.calcCratio();
    truecratC2.calcCratio();
    truecratCu.calcCratio();
    truecratSn.calcCratio();
    simcratC1.calcCratio();
    simcratC2.calcCratio();
    //truecratSn.multiplotCratio(truecratCu, truecratC1);
    broadsimC1.calcDpt();
    broadsimC2.calcDpt();
    broadtrueC1.calcDpt();
    broadtrueC2.calcDpt();
    broadtrueCu.calcDpt();
    broadtrueSn.calcDpt();
    broadtrueSn.writeMatrixToFile("matrixSn.txt");
    broadtrueSn.multiplotDpt(broadtrueCu, broadtrueC1);
    broadtrueSn.DrawOnlyptandz("Snptz");

    return 0;
}
