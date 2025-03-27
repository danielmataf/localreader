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
#include "TCanvas.h"


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
    EventReader RGD_CxC(filenamesCxC);
    EventReader RGD_LD2(filenamesLD2);   
    EventReader RGD_CuSn(filenamesCuSn);   
    EventReader Sim_CxC(simufilesCxC);



    std::optional<Event> simCxC;
    std::optional<Event> simCxC_MC; //testing this if we can loop twice top read MC events//
    std::optional<Event> testCxC;
    std::optional<Event> testLD2;
    std::optional<Event> testCu;
    std::optional<Event> testSn;



    CutSet testC1cuts;  //C1 RGD
    CutSet testC2cuts;  //C2 RGD  
    CutSet simC2cuts;  
    CutSet testLD2cuts;
    CutSet testCucuts;
    CutSet testSncuts;



    //Defining cuts here. Vz cuts for SIM and RGD may differ -> check Mon Vz
    //simC1cuts.SetCutVz(Constants::RcutminVzC1sim,Constants::RcutminVzC1sim);     //vz cut for C1 target
    //testC1cuts.SetCutVz(Constants::RcutminVzC1data,Constants::RcutmaxVzC1data);     //previous cut, b4jan2025. We change for uncalib data (?) peaks shifted of -1cm
    testC1cuts.SetCutVz(Constants::RcutminVzC1data,Constants::RcutmaxVzC1data);     //Changing values in constants, also widening the cut
    testC1cuts.SetCutGen4Rat();

    testC2cuts.SetCutVz(Constants::RcutminVzC2data, Constants::RcutmaxVzC2data);    //vz cut for C2 target
    testC2cuts.SetCutGen4Rat();
    simC2cuts.SetCutVz(Constants::RcutminVzC2sim    , Constants::RcutmaxVzC2sim);     //vz cut for C2 target
    simC2cuts.SetCutGen4Rat();

    testLD2cuts.SetCutVz(Constants::RcutminVzLD2,Constants::RcutmaxVzLD2);     //vz cut for LD2 target
    testLD2cuts.SetCutGen4Rat();

    testCucuts.SetCutVz(Constants::RcutminVzCu,Constants::RcutmaxVzCu);     //vz cut for Cu target
    testCucuts.SetCutGen4Rat();
    
    testSncuts.SetCutVz(Constants::RcutminVzSndata,Constants::RcutmaxVzSndata);     //vz cut for Sn target
    testSncuts.SetCutGen4Rat();

    int sumevts = 0;

    Monitoring monSimC2(simC2cuts, "C2_sim");     //This needs to be figured out ASAP
    Monitoring monTestC1(testC1cuts, "C1_RGD");  //This needs to be figured out ASAP
    Monitoring monTestC2(testC2cuts, "C2_RGD");  //This needs to be figured out ASAP
    Monitoring monTestLD2(testLD2cuts, "LD2_RGD");
    Monitoring monTestCu(testCucuts, "Cu_RGD");
    Monitoring monTestSn(testSncuts, "Sn_RGD");

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


    const int totalevts = 19000;

    //6 vectors for phi one for each sector and 6 vectors for theta, same, then with another target
    std::vector<float> theta_sn[6], phi_sn[6];
    std::vector<float> theta_cxc[6], phi_cxc[6];

    for (int i = 0; i < totalevts; ++i) {
        auto testCxC = RGD_CxC.ProcessEventsInFile();
        auto testSn = RGD_CuSn.ProcessEventsInFile();

        if (testCxC.has_value()) {
            Event evtCxC = testCxC.value();
            evtCxC.SetTargetType(1);
            evtCxC.calcAll();

            monTestC2.FillHistogramswCuts(evtCxC);  //usual procedure

            int sec = evtCxC.electron.GetCalSector();   //sec is for sectors 
            if (sec >= 1 && sec <= 6) {
                int s = sec - 1;
                theta_cxc[s].push_back(evtCxC.electron.GetMomentum().Theta() * 180 / Constants::PI);
                phi_cxc[s].push_back(evtCxC.electron.GetMomentum().Phi() * 180 / Constants::PI + 180);
            }
        }

        if (testSn.has_value()) {
            Event evtSn = testSn.value();
            evtSn.SetTargetType(1);
            evtSn.calcAll();

            monTestSn.FillHistogramswCuts(evtSn);

            int sec = evtSn.electron.GetCalSector();
            if (sec >= 1 && sec <= 6) {
                int s = sec - 1;
                theta_sn[s].push_back(evtSn.electron.GetMomentum().Theta() * 180 / Constants::PI);
                phi_sn[s].push_back(evtSn.electron.GetMomentum().Phi() * 180 / Constants::PI + 180);
            }
        }

        files.displayProgress(i + 1, totalevts);
    }

    //rescalc per sec
    TH1F* h_phi_res[6];
    TH1F* h_theta_res[6];
    for (int s = 0; s < 6; ++s) {
        h_phi_res[s] = new TH1F(Form("phi_res_s%d", s+1), Form("Phi Resolution Sector %d", s+1), 100, -5, 5);
        h_theta_res[s] = new TH1F(Form("theta_res_s%d", s+1), Form("Theta Resolution Sector %d", s+1), 100, -5, 5);

        size_t N = std::min(phi_cxc[s].size(), phi_sn[s].size());
        for (size_t i = 0; i < N; ++i) {
            float dphi = phi_cxc[s][i] - phi_sn[s][i];
            float dtheta = theta_cxc[s][i] - theta_sn[s][i];
            h_phi_res[s]->Fill(dphi);
            h_theta_res[s]->Fill(dtheta);
        }
    }

    //res histos
    TCanvas c("res_canvas", "Angular Resolutions", 1200, 800);
    c.Divide(3, 2);
    for (int s = 0; s < 6; ++s) {
        c.cd(s + 1);
        h_phi_res[s]->GetXaxis()->SetRangeUser(-1, 1);
        h_phi_res[s]->SetLineColor(kBlue);
        h_phi_res[s]->Draw("hist");
    }
    c.SaveAs("phi_resolution_per_sector.pdf");
    monTestC2.SaveHistRoot("janC2_test");
    monTestC2.DrawHistograms("monC2_test");
    monTestSn.SaveHistRoot("janSn_test");
    monTestSn.DrawHistograms("SnUrgent");

    return 0;
}
