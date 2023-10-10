#ifndef MONITORING_H
#define MONITORING_H

#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include "Event.h" // Include your Event class header here
#include "Cuts.h"

class Monitoring {
public:
    Monitoring(CutSet);

    void FillHistograms(const Event& event);


private:
    CutSet cut1;
    TH1F h_Q2;
    TH1F h_xb ;
    TH1F h_y ;
    TH1F h_nu ;
    TH1F h_W2 ;
    TH1F h_z ;
    TH1F h_pt2 ;
    TH1F h_phih ;
    //add more histograms for other variables here

    TFile* outputFile;
    //creating a root file 
};

#endif 