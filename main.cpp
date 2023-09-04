
#include "reader.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TCanvas.h>
#include <TH2.h>
#include <TF1.h>
#include <TRatioPlot.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TLine.h>
#include <TArrow.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define PI 3.14159265
//#include "dcFiducialCuts.hpp"
#include "linspace.hpp"
#include "variablecalc.hpp"
#include "histograms.hpp"
#include "vectors.hpp"
#include "variables.hpp"

VarHistograms myHists;
MCHistograms MCHists;
hist_electron elHists;
hist_pip pipHists;
hist_hadX hadHists;
MultipRatio R_hists;
h_Dpt hist_Dpt;
h_cratio hist_cratio;
h_sratio hist_sratio;
h_c2ratio hist_c2ratio;
WeightHistograms weights;
VarianceHistograms variance;
PhysicsVectors physicsVecs;
MyVariables variables;


using namespace std;
int main() {
    
    //cart coord for Pion+
    
    //
    
    int nubin = 100;  //--------------------------------------------------------------------binning -->USE 10 FOR RATIO PLOTS
    //===== Histogram creation =====//
    TH1F* h_scapiE = new TH1F("nrg(pi)", "nrg(pi)", nubin, 0.0, 10);
    // REC histograms  --->  see histograms.hpp
    //== counters ==//
    int counter_onlye =0;   //only electron counter, no coincidence
    int counter_el = 0;   //electron counter
    int counter_pi = 0;   //pion (+) counter
    int counter_pim = 0;  //pion (-) counter
    int counter_pro= 0;     //proton counter
    int counter_neu = 0;   //neutron counter
    int j ;   //photon saver
    int k;
    double e_sup_pho = 0;
    int coinci_De = 0;
    //ID	(chi2, eventID for each particle)
    double IDel, IDpip;

    //== bools ==/
    bool val_e, val_pip, val_pim, val_neu, val_pro,val_pit; //resets condition for electrons, and allow to count photons

    //== values  ==//
    double m_e = 0.000511, Mn = 0.939565, Mpro = 0.938272, M_He = 3.727379, m_pi = 0.134976, M_D2 = 1.875;
	//ELECTRON MASS ADJUST /!\ (???) //it was at 0.000000511 (not correct i guess)
    
    int typeoffile; 
    cout<<"which type of file ? 1=D , 2 = Sn, 3=Cu"<<endl;
    cin>> typeoffile;
    const int files = 2;    //100 in ifarm
    int RECoMC; 
    cout<<"which type of analysis ? 1=REC , 2 = MC"<<endl;
    cin>> RECoMC;
    //RECoMC = 1;

    std::string basePath;
    std::string EndPath;
    if (typeoffile==1)
    {
        //basePath = "/volatile/clas12/dmat/gen/Deut/";       //path of files in ifarm
        //EndPath = "/sidis_mc-master/r_out_rgd10k_d.hipo";
        basePath = "../../files2read/r_eD-0";       //path of files in ifarm
        EndPath = ".hipo";

    }
    else
    {
        //basePath = "/volatile/clas12/dmat/gen/Sn/";
        //EndPath = "/sidis_mc-master/r_out_rgd10k_Sn.hipo";
        basePath = "../../files2read/r_eSn-0";       //path of files in ifarm
        EndPath = ".hipo";
        
    }

    std::vector<std::string> filelist;
    for (int filenb = 1; filenb <= files; filenb++) {
        std::string filePath = basePath + std::to_string(filenb) + EndPath;
        filelist.push_back(filePath);
    }
    
    
    
    for (int filenbr = 0; filenbr<files ; filenbr++){
        //filelist[filenbr]; is a const char....most probably...
        hipo::reader reader;
	    cout<<(filenbr*100)/files<<"% complete"<<endl;
        //const char *filename = "../../path/to/file/here";
        reader.open(filelist[filenbr].c_str());
        hipo::dictionary factory;
        reader.readDictionary(factory);
        hipo::bank RECgen;
        if (RECoMC==1)
        {
            hipo::bank RECgen(factory.getSchema("REC::Particle"));  //call REC Bank    
        }
        else
        {
            hipo::bank RECgen(factory.getSchema("MC::Particle"));  //call MC Bank
        }
        hipo::bank REC_Particle_bank(factory.getSchema("REC::Particle"));  //call REC Bank
        hipo::bank MC_Particle_bank(factory.getSchema("MC::Particle"));       //call MC Bank
        hipo::bank RUNevt(factory.getSchema("RUN::config"));    //call RUN bank 4 evt_nbr
        //***************//
        hipo::bank RECTraj(factory.getSchema("REC::Traj"));
	    hipo::bank RECCal(factory.getSchema("REC::Calorimeter"));
        //**************//
        hipo::event event;
        while (reader.next()) {
            e_sup_pho = 0; //reset the photon max energy at the beginning of every new event.
            j = 0;   //reset the photon saver
            k = 0;
            val_e=false;
	        val_pip = false;
	        val_pim = false;
	        val_pit = false;
	        val_neu = false;
	        val_pro = false;
            reader.read(event);
            //event.getStructure(RECgen);  //general particle bank call for MC and REC
            event.getStructure(REC_Particle_bank);
            event.getStructure(MC_Particle_bank);
            event.getStructure(RUNevt);
	        event.getStructure(RECTraj);
	        event.getStructure(RECCal);
            int numrows    = RECgen.getRows();
            int numRECrows = REC_Particle_bank.getRows();
            int numMCRows  = MC_Particle_bank.getRows();
            int numRUNRows  = RUNevt.getRows();
            int numRECcal  = RECCal.getRows();
	        //cout<<numRECcal<<endl;
	        
	        int numTRAJrows= RECTraj.getRows();
	        if (numTRAJrows==-1 || 0){continue;}
            //-----Q2//
            //==   Vectors ==//
            //Vectors are defined HERE in order to be reset after every iteration
            physicsVecs.v_incil->SetPxPyPzE(0.0, 0.0, 11.0, 11.0); //has no MC version, since it is a predefined vector, REC and MC versions are the same
            TVector3 p_incil = physicsVecs.v_incil->Vect();  //3D
            physicsVecs.v_nucleontarg->SetPxPyPzE(0.0, 0.0, 0.0, Mn);
	        int evt_nbr = RUNevt.getInt("run",0); //goes inside the event loop
            //======GET REC EVENTS=======//*
	        for (int i = 0; i < numRECrows; ++i) {
		        if (REC_Particle_bank.getInt("pid", i)==211 ){
                    //Get the REC PION      //detect PION
		            val_pip = true;
                    counter_pi += 1;
		            IDpip = REC_Particle_bank.getFloat("chi2pid", i);
		            physicsVecs.v_scapip->SetPx(REC_Particle_bank.getFloat("px", i));
                    physicsVecs.v_scapip->SetPy(REC_Particle_bank.getFloat("py", i));
                    physicsVecs.v_scapip->SetPz(REC_Particle_bank.getFloat("pz", i));
                    TVector3 p_scapip = physicsVecs.v_scapip->Vect();
                    physicsVecs.v_scapip->SetVectM(p_scapip,m_pi);
                }
		        if (REC_Particle_bank.getInt("pid", i)==-211 ){
                    //Get the REC PION      //detect PION MINUS
                    val_pim = true;
		            counter_pim += 1;
		            physicsVecs.v_scapim->SetPx(REC_Particle_bank.getFloat("px", i));
                    physicsVecs.v_scapim->SetPy(REC_Particle_bank.getFloat("py", i));
                    physicsVecs.v_scapim->SetPz(REC_Particle_bank.getFloat("pz", i));
                    TVector3 p_scapim = physicsVecs.v_scapim->Vect();
                    physicsVecs.v_scapim->SetVectM(p_scapim,m_pi);
                }
		        if (REC_Particle_bank.getInt("pid", i)==2212 ){
                    //Get the REC PROTON
                    val_pro = true;
		            counter_pro += 1;
		            physicsVecs.v_scapro->SetPx(REC_Particle_bank.getFloat("px", i));
                    physicsVecs.v_scapro->SetPy(REC_Particle_bank.getFloat("py", i));
                    physicsVecs.v_scapro->SetPz(REC_Particle_bank.getFloat("pz", i));
                    TVector3 p_scapro = physicsVecs.v_scapro->Vect();
                    physicsVecs.v_scapro->SetVectM(p_scapro,Mpro);
                }
		        if (REC_Particle_bank.getInt("pid", i)==2112 ){
                    //Get the REC NEUTRON      //
                    val_neu = true;
		            counter_neu += 1;
		            physicsVecs.v_scaneu->SetPx(REC_Particle_bank.getFloat("px", i));
                    physicsVecs.v_scaneu->SetPy(REC_Particle_bank.getFloat("py", i));
                    physicsVecs.v_scaneu->SetPz(REC_Particle_bank.getFloat("pz", i));
                    TVector3 p_scaneu = physicsVecs.v_scaneu->Vect();
                    physicsVecs.v_scaneu->SetVectM(p_scaneu,Mn);
                }
                if (REC_Particle_bank.getInt("pid", i)==11 && val_e==false ){
                    //Get the REC electron      //detect electron
                    val_e = true;
                    counter_el += 1;
        		    IDel = REC_Particle_bank.getFloat("chi2pid", i);
                    physicsVecs.v_scal->SetPx(REC_Particle_bank.getFloat("px", i));
                    physicsVecs.v_scal->SetPy(REC_Particle_bank.getFloat("py", i));
                    physicsVecs.v_scal->SetPz(REC_Particle_bank.getFloat("pz", i));
                    TVector3 p_scal = physicsVecs.v_scal->Vect();
                    physicsVecs.v_scal->SetVectM(p_scal,m_e);

                }
            }
            				/////////////////////
                // ---BEGIN CALC---//
				/////////////////////
                //get the vqlues on qny possible vector contained in the event
                //suppose: ward against empty vectors due to the condition on val_e =true for ex
                //theta_el = v_scal->Theta();	//coord_el
            
			
                //======= START ANALYSIS 4 ANY TYPE OF FILE =======//
		        //==== CALCS and FILLINGS for MC (Deuterium) ===//
                    //add an option in cin at the beginning of the code in order to have an option for MC analysis.... TBD !
                    //DELETED, if needed check previous versions//

	            //==== CALCS and FILLINGS for REC ===//		---[De]
	        if (val_e==true){	
                counter_onlye +=1;
                variables.theta_el = cordtheta(physicsVecs.v_scal);   
                variables.phi_el = physicsVecs.v_scal->Phi();		//coord_el
                *physicsVecs.v_vipho = *physicsVecs.v_incil - *physicsVecs.v_scal;          //-----keep
                TVector3 p_scal = physicsVecs.v_scal->Vect();
		        variables.p_el = p_scal.Mag();
                variables.theta_e= acos( p_incil.Dot(p_scal) / ( p_incil.Mag() * p_scal.Mag() ) ); //useful

                cout<<"theta_e"<<variables.theta_e<<endl;

		        variables.Q2 = calcQ(physicsVecs.v_incil,physicsVecs.v_scal);
                variables.gamnu= calcgamnu(physicsVecs.v_incil, physicsVecs.v_scal);
                variables.y = calcy(physicsVecs.v_incil,variables.gamnu);
                variables.x_b = calcx(variables.Q2, Mn, variables.gamnu);
		        variables.W = calcW(variables.Q2, Mn, variables.gamnu);
                TVector3 p_vipho= physicsVecs.v_vipho->Vect();
		            
			    TVector3 n1 = p_vipho.Cross(p_incil);
			        
                variables.px_scal=physicsVecs.v_scal->Px();
 		        variables.py_scal=physicsVecs.v_scal->Py();
		        variables.pz_scal=physicsVecs.v_scal->Pz();
		            
                TVector3 p_scaph= physicsVecs.v_scaph->Vect();
		        variables.p_pho = p_scaph.Mag();
                variables.Epho = physicsVecs.v_scaph->E();
                    
		        //variable.cone_ang =acos(( p_misspho.Dot(p_scaph) ) / ( p_misspho.Mag()*p_scaph.Mag() )) *180/PI;//?
		        //bool dccut = DC_cuts_electron(RECTraj, 0);				//NEEDED (TBD)
	            variables.scalX = physicsVecs.v_scal->X();
		        variables.scalY = physicsVecs.v_scal->Y();			
	            if (variables.Q2>1.5){
		            R_hists.hist_Q_e->Fill(variables.Q2);
		            if (variables.gamnu<8.5 && variables.gamnu>2.5){
                        R_hists.hist_v_e->Fill(variables.gamnu);
				        if (variables.W>2) {		//!!!!!!!!//
				        }
			        }
		        }
                if (val_pip==true  ){	//this condition involves the condition in val_e
                    coinci_De += 1;
                    variables.theta_pip = physicsVecs.v_scapip->Theta()  ;	//bof
                    variables.phi_pip = physicsVecs.v_scapip->Phi() ;	//bof
                    variables.z = calcz(physicsVecs.v_scapip,variables.gamnu);
                    variables.t = physicsVecs.v_diff->M2()  ;
                        
                    TVector3 p_scapip = physicsVecs.v_scapip->Vect();
                    TVector3 n2 = p_vipho.Cross(p_scapip);
                    TVector3 n3 = n1.Cross(n2);
                    *physicsVecs.v_diffmass= *physicsVecs.v_incil + *physicsVecs.v_nucleontarg -*physicsVecs.v_scal - *physicsVecs.v_scapip;	//this vector corresponds to the 	X-HADRON
		            TVector3 p_diffmass= physicsVecs.v_diffmass->Vect();
                    variables.mmass = physicsVecs.v_diffmass->M2(); //SQUARED		//bof
                    variables.px_scapi=physicsVecs.v_scapip->Px();
 		            variables.py_scapi=physicsVecs.v_scapip->Py();
		            variables.pz_scapi=physicsVecs.v_scapip->Pz();
                    variables.P_t =(n2.Mag() )/ (p_vipho.Mag() );
		            TVector3 crossX = p_vipho.Cross(p_diffmass);
                    variables.angle_norm=( n1.Dot(n2) ) / ( n1.Mag()*n2.Mag() );
                    variables.sign = n3.Dot(p_vipho);
                    variables.P_trX = calcP_trX(physicsVecs.v_vipho,physicsVecs.v_diffmass );
		            variables.pipX = physicsVecs.v_scapip->X();
                    variables.pipY = physicsVecs.v_scapip->Y();
                    
                    if (variables.sign<0 ){
	                    variables.phih =acos(variables.angle_norm) *180.0/PI+180;
                        //if (filenbr%2 == 0) {h_Phphi->Fill(phih);}	//histogral filling for eventual BSA TBD
                        //if (filenbr%2 != 0) {h_Nhphi->Fill(phih);}	//keeping them, useful
                        //phih = phih +180;
                        variables.rawphih =acos(variables.angle_norm)*180.0/PI;
                        variables.rawphih2 = variables.angle_norm;
	                }
                    if (variables.sign>0 ){					
	                    variables.phih =acos(variables.angle_norm) *-180.0/PI+180;
                        //if (filenbr%2 == 0) {h_Phphi->Fill(phih);}	//TBD
                        //if (filenbr%2 != 0) {h_Nhphi->Fill(phih);}	//
                        //phih = phih +180;
                        variables.rawphih =-acos(variables.angle_norm)*180.0/PI;
                        variables.rawphih2 = variables.angle_norm;
                    } 
                    
                    elHists.h_theta_el->Fill(variables.theta_el*180/PI);
		            elHists.h_phi_el->Fill(variables.phi_el*180/PI);
		            elHists.coor2D_el->Fill(variables.phi_el*180/PI,variables.theta_el*180/PI);
		            pipHists.h_theta_pip->Fill(variables.theta_pip*180/PI);
		            pipHists.h_phi_pip->Fill(variables.phi_pip*180/PI);
		            pipHists.coor2D_pip->Fill(variables.phi_pip*180/PI,variables.theta_pip*180/PI);
		            elHists.hID_el->Fill(IDel);
                    if (abs(IDel)<3){		//first condition
		                pipHists.hID_pip->Fill(IDpip);
		                if (abs (IDpip)<2){
		                    //hist_Q->Fill(Q2);			//-------------------Q---//
		                    if (variables.Q2>1.5){
                                //myHists.hist1->Fill(y);
                                myHists.hist_y->Fill(variables.y);
		                        if (variables.y<0.85 && variables.y>0.25){
	    	    	                //hist_z->Fill(z);		//-------------------z---//
			                        if (variables.z>0.3 && variables.z<0.7){
			    	                    //if (dccut == true){	// /!\ implementing DC CUT!!!!		(compatibility TBD)
				                            if (variables.theta_el*180/PI>6){	// /!\ implementing CUT on theta coordinate for electrons!!!!
            		        		        //hist_v->Fill(gamnu);	//-------------------v---//
	    	    	   	                        myHists.hist_x->Fill(variables.x_b);
            	    	    	                myHists.hist_W->Fill(variables.W);
					                            if (variables.W>2) {
	    	    	    	    	                myHists.hist_phih->Fill(variables.phih);
	    	    	    	    	                //hist_P_t->Fill(P_t);	//-------------------pt--//
	    	    	    	    	                myHists.hist_MX->Fill(variables.mmass);
					                                //h_cphi->Fill(cos(phih));
					                                if (variables.mmass>1.2){				//To remove the proton. Don't forget to draw the ct in the plot
				        	                            elHists.hpx_el->Fill(variables.px_scal);
				        	                            elHists.hpy_el->Fill(variables.py_scal);
				        	                            elHists.hpz_el->Fill(variables.pz_scal);
						                                myHists.hist_t->Fill(variables.t);
				        	                            pipHists.hpx_pi->Fill(variables.px_scapi);
				        	                            pipHists.hpy_pi->Fill(variables.py_scapi);
				        	                            pipHists.hpz_pi->Fill(variables.pz_scapi);
				        	                            hadHists.hpx_X->Fill(physicsVecs.v_diffmass->Px());
				        	                            hadHists.hpy_X->Fill(physicsVecs.v_diffmass->Py());
				        	                            hadHists.hpz_X->Fill(physicsVecs.v_diffmass->Pz());
  				        	                            hadHists.h_P_trX->Fill(variables.P_trX);
						                                elHists.h_scalE->Fill(physicsVecs.v_scal->E());
                        						        h_scapiE->Fill(physicsVecs.v_scapip->E());      //replace this one with a class histogram in hpp TBD
					                                    if (physicsVecs.v_scapip->E()<3 && physicsVecs.v_scapip->E()>1){            //+++++++++++++++++++++++++++//
				        	                                //if (lv>90 && lw>90) {				//calorimeter cut on 9cm (?)
						                                    //h_cphi->Fill(cos(phih));			//Cos & Sin analysis
						                                    myHists.hist_Q->Fill(variables.Q2);		//-----------------------Q--///
						                                    myHists.hist_v->Fill(variables.gamnu);		//-----------------------v--///
						                                    myHists.hist_z->Fill(variables.z);		//-----------------------z--///
						                                    myHists.hist_P_t->Fill(variables.P_t);		//----------------------pt--//
                                                        
                                                            //all P_t used below have to be squared even in the weights
                                                            weights.hist_wQp->Fill(variables.Q2,variables.P_t*variables.P_t);
                                                            weights.hist_wvp->Fill(variables.gamnu,variables.P_t*variables.P_t);
                                                            weights.hist_wzp->Fill(variables.z,variables.P_t*variables.P_t);
                                                            weights.hist_wQcos->Fill(variables.Q2 , cos(variables.phih));
                                                            weights.hist_wvcos->Fill(variables.gamnu , cos(variables.phih));
                                                            weights.hist_wP_tcos->Fill( variables.P_t*variables.P_t, cos(variables.phih));
                                                            weights.hist_wzcos->Fill( variables.z, cos(variables.phih));
                                                            weights.hist_wQsin->Fill( variables.Q2, sin(variables.phih));
                                                            weights.hist_wvsin->Fill( variables.gamnu , sin(variables.phih));
                                                            weights.hist_wP_tsin->Fill(variables.P_t*variables.P_t, sin(variables.phih) );
                                                            weights.hist_wzsin->Fill(variables.z, sin(variables.phih));
                                                            weights.hist_wQc2->Fill(variables.Q2, cos(2*variables.phih)  );
                                                            weights.hist_wvc2->Fill( variables.gamnu, cos(2*variables.phih));
                                                            weights.hist_wP_tc2->Fill(variables.P_t*variables.P_t, cos(2*variables.phih));
                                                            weights.hist_wzc2->Fill( variables.z, cos(2*variables.phih));
                                                            variance.hist_varQp->Fill(variables.Q2,pow(variables.P_t,4));
                                                            variance.hist_varvp->Fill(variables.gamnu,pow(variables.P_t,4));
                                                            variance.hist_varzp->Fill(variables.z,pow(variables.P_t,4));
                                                            variance.hist_varQcos->Fill(variables.Q2 , pow(cos(variables.phih),2));
                                                            variance.hist_varvcos->Fill(variables.gamnu , pow(cos(variables.phih),2));
                                                            variance.hist_varP_tcos->Fill( variables.P_t*variables.P_t, pow(cos(variables.phih),2));
                                                            variance.hist_varzcos->Fill( variables.z, pow(cos(variables.phih),2));
                                                            variance.hist_varQsin->Fill( variables.Q2, pow(sin(variables.phih),2));
                                                            variance.hist_varvsin->Fill( variables.gamnu ,pow(sin(variables.phih),2));
                                                            variance.hist_varP_tsin->Fill(variables.P_t*variables.P_t, pow(cos(variables.phih),2) );
                                                            variance.hist_varzsin->Fill(variables.z, pow(sin(variables.phih),2));
                                                            variance.hist_varQc2->Fill(variables.Q2, pow(cos(2*variables.phih),2) );
                                                            variance.hist_varvc2->Fill( variables.gamnu, pow(cos(2*variables.phih),2));
                                                            variance.hist_varP_tc2->Fill(variables.P_t*variables.P_t, pow(cos(2*variables.phih),2));
                                                            variance.hist_varzc2->Fill( variables.z, pow(cos(2*variables.phih),2));



                                                            hist_Dpt.h_pQ->Fill(variables.P_t*variables.P_t,variables.Q2);
                                                            hist_Dpt.h_pv->Fill(variables.P_t*variables.P_t,variables.gamnu);
                                                            hist_Dpt.h_pz->Fill(variables.P_t*variables.P_t,variables.z);
                                                            hist_cratio.h_cQ->Fill(cos(variables.phih),variables.Q2);
                                                            hist_cratio.h_cv->Fill(cos(variables.phih),variables.gamnu);
                                                            hist_cratio.h_cz->Fill(cos(variables.phih),variables.z);
                                                            hist_cratio.h_cpt->Fill(cos(variables.phih),variables.P_t); 
                                                            hist_sratio.h_sQ->Fill(sin(variables.phih),variables.Q2);
                                                            hist_sratio.h_sv->Fill(sin(variables.phih),variables.gamnu);
                                                            hist_sratio.h_sz->Fill(sin(variables.phih),variables.z);
                                                            hist_sratio.h_spt->Fill(sin(variables.phih),variables.P_t);
                                                            hist_c2ratio.h_c2Q->Fill(cos(2*variables.phih),variables.Q2);
                                                            hist_c2ratio.h_c2v->Fill(cos(2*variables.phih),variables.gamnu);
                                                            hist_c2ratio.h_c2z->Fill(cos(2*variables.phih),variables.z);
                                                            hist_c2ratio.h_c2pt->Fill(cos(2*variables.phih),variables.P_t*variables.P_t );
					            	                        //q condition was here for the MC comparisons DELETED
							                            }
                       			       	                // }					//end if on MC vs REC
 			 	            	                    }
					                            }
				                            }
			   	                        //}			//DC cut (compatibility TBD)
			                        }
		                        }
		                    }
		                }
		            }
                }
		    }       //condition for true e 
            
		    //===============================================   Sn loop erased   =============================================//
        }
    }

	cout<<"D : "<<coinci_De<<endl;
    //the couts for Sn and D are useless now

    std::map<int, std::string> fileNames = {
        {1, "../output_D.root"},
        {2, "../output_Sn.root"},
        {3, "../output_Cu.root"}
    };

    if (fileNames.find(typeoffile) == fileNames.end()) {
        cout << "Invalid option" << endl;
        // Handle invalid option error or provide appropriate action
    }
    else {
        std::string filename = fileNames[typeoffile];
        TFile* file = new TFile(filename.c_str(), "RECREATE");

        myHists.hist_Q->Write();
        myHists.hist_v->Write();
        myHists.hist_y->Write();
        myHists.hist_x->Write();
        myHists.hist_W->Write();
        myHists.hist_phih->Write();
        myHists.hist_t->Write();
        myHists.hist_P_t->Write();
        myHists.hist_z->Write();
        myHists.hist_MX->Write();
        R_hists.hist_Q_e->Write();
        R_hists.hist_v_e->Write();
        hist_Dpt.h_pQ->Write();
        hist_Dpt.h_pv->Write();
        hist_Dpt.h_pz->Write();
        hist_cratio.h_cQ->Write();
        hist_cratio.h_cv->Write();
        hist_cratio.h_cz->Write();
        hist_cratio.h_cpt->Write(); 
        hist_sratio.h_sQ->Write();
        hist_sratio.h_sv->Write();
        hist_sratio.h_sz->Write();
        hist_sratio.h_spt->Write();
        hist_c2ratio.h_c2Q->Write();
        hist_c2ratio.h_c2v->Write();
        hist_c2ratio.h_c2z->Write();
        hist_c2ratio.h_c2pt->Write();
        weights.hist_wQp->Write();
        weights.hist_wvp->Write();
        weights.hist_wzp->Write();
        weights.hist_wQcos->Write();
        weights.hist_wvcos->Write();
        weights.hist_wP_tcos->Write();
        weights.hist_wzcos->Write();
        weights.hist_wQsin->Write();
        weights.hist_wvsin->Write();
        weights.hist_wP_tsin->Write();
        weights.hist_wzsin->Write();
        weights.hist_wQc2->Write();
        weights.hist_wvc2->Write();
        weights.hist_wP_tc2->Write();
        weights.hist_wzc2->Write();
        variance.hist_varQp->Write();
        variance.hist_varvp->Write();
        variance.hist_varzp->Write();
        variance.hist_varQcos->Write();
        variance.hist_varvcos->Write();
        variance.hist_varP_tcos->Write();
        variance.hist_varzcos->Write();
        variance.hist_varQsin->Write();
        variance.hist_varvsin->Write();
        variance.hist_varP_tsin->Write();
        variance.hist_varzsin->Write();
        variance.hist_varQc2->Write();
        variance.hist_varvc2->Write();
        variance.hist_varP_tc2->Write();
        variance.hist_varzc2->Write();
        TTree tree("treecounter","treecounter");
        tree.Branch("counter_onlye", &counter_onlye, "counter_onlye/I");
        tree.Fill();
        file->Write();

        file->Close();
    }
    // we have one only loop that could be used to calc all variables for any nuclear target 

    
}
