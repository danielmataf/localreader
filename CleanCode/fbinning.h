#ifndef FBINNING_H
#define FBINNING_H

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include "Event.h" // 
#include "CutSet.h"
#include "constants.h"

//This class is made to facilitate corrections. Goal is a priori to defnie an uniform binning for important variables in TMD
//These somewhat uniform binning will be used later for corrections such as unfolding. It is important to store the edges of this binning porpoerly
//variables we are interested in are Q2, xb, z and pt2. We will also be interested in the binning of phi_h. these var will be stored in a matrix after kincuts

class fbinning {
public:
    fbinning(CutSet, const std::string& targetName); //adding name to constructor to add unique names to histograms
    //fbinning(CutSet);
    ~fbinning();          


    void FillwCuts(const Event& );
    void FillswCutsMC(const Event& );

private:
    CutSet cut1;
    std::string targetName;     //adding target name to add to hist names
    
    std::vector<std::vector<std::vector<double>>> Matrix;    //three vectors for 3D matrix
    std::vector<std::vector<std::vector<double>>> errorMatrix;

    std::vector<double> eQ2, exb, ez, ept2, ephih; //defining a matrix for the the binning of each variable we want


    TFile outputFile;
    //creating a root file 


};

#endif 
