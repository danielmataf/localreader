#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLatex.h>  

#include <TPad.h>
#include <fstream>
#include <math.h>
#include <TPDF.h>
#include "Event.h" 
#include "fbinning.h"
#include "CutSet.h"
#include <TMultiGraph.h>
#include <TLegend.h>


fbinning::fbinning(CutSet a, const std::string& targetName): 
    cut1(a)
    targetName(targetName),
    SratioMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    errorSratioMatrix(Cratiobin, std::vector<std::vector<double>>(Cratiobin, std::vector<double>(Cratiobin, 0.0))),
    std::vector<double> eQ2, exb, ez, ept2, ephi_h; 



    {
    //cuta = cutsA;
    }



    void fbinning::FillwCuts(const Event& event ){
        if (cut1.PassCutsDetectors(event)) {
            if (cut1.PassCutsElectrons(event)==true) {
                for (const Particle& hadron : event.GetHadrons()) {
                    if (cut1.PassCutsHadrons(hadron)) {
                        if (hadron.GetPID() == Constants::PION_PLUS_PID) {
                            eQ2.push_back(event.GetElectron.GetQ2);
                            exb.push_back(event.GetElectron.GetQ2);
                            ez.push_back(hadron.Getz);
                            ept2.push_back(hadron.Getpt2);
                            ephih.push_back(hadron.Getphih);
                        }
                    }
                }
            }
        }
    }

