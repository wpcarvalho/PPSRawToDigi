/*************************************************
  
   .L InternalMisalGenerator.C  
   InternalMisalGenerator(10,/home/users/berretti/SL/WorkingArea/Analysis/TestRandomMisal)

   Output: N_scenarios misalignment configuration files in the format "DetId: "DX:" "DY:" "DPhi:"
           Name will be : RandomInternalMisal_ScenarioI.dat where I=0..N_Scenario. 

   /afs/cern.ch/exp/totem/scratch/berretti/tmp/testSplitMerge/CMSSW_3_1_1/src/SimTotem/T2Digitizer/data/InternalMisalignSimulation
   Authors: M.Berretti
  
   Feb 2011, TOTEM Experiment

process.T2TrkBasedIntAl.IdreferencedetMille=cms.vuint32(0,6,12,18,23,28,32,37)
process.T2TrkBasedIntAl.XreferencedetMille=cms.vdouble(-0.094,0.105, 0.027,-0.057, -0.102,0.0297, -0.1542,0.0036)
process.T2TrkBasedIntAl.YreferencedetMille=cms.vdouble(-0.0013,-0.2048,-0.0588,0.0558,-0.0315,0.0557,-0.108,-0.2723)
*************************************************/

#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <cstdio> 
#include <sstream>
#include <TSystem.h> 
#include <TRandom.h> 
//#include <boost>


//////////////////////////////////////////////////////
double XYTollerance=1.0; //Mechanical tollerance in mm
//////////////////////////////////////////////////////


void InternalMisalGenerator(int N_scenarios, std::string outputFolder);




void InternalMisalGenerator(int N_scenarios, std::string outputFolder){



  for(unsigned int i=0;i<=N_scenarios;i++){
    
    ofstream outputFileDispl;
    stringstream st_I;//create a stringstream
    st_I << (i);//add number to the stream
    
    std::string Absnomefile= outputFolder+"/RandomInternalMisal_Scenario"+st_I.str()+".dat";
    

    outputFileDispl.open(Absnomefile.c_str());

    
    TRandom *r0 = new TRandom();
    r0->SetSeed(i);

    for(unsigned int det=0;det<40;det++){

     
      double x = r0->Uniform(-(1.0)*XYTollerance,XYTollerance);
      double y = r0->Uniform(-(1.0)*XYTollerance,XYTollerance);
      double phi = 0.;
      
      if(det==0){
	x =0.;
	y =0.;
      }
      if(det==6){
	x =0.;
	y =0.;
      }
      if(det==12){
	x =0.;
	y =0.;
      }
      if(det==18){
	x =0.;
	y =0.;
      }
      if(det==23){
	x =0.;
	y =0.;
      }
      if(det==28){
	x =0.; 
	y =0.;
      }
      if(det==32){
	x =0.;
	y =0.;
      }
      if(det==37){
	x =0.;
	y =0.;
      }


      outputFileDispl<<"DetId: "<<det<<"  DX: "<<x<<"  DY: "<<y<<"  DPhi: "<<phi<<"\n";

    }
    

    outputFileDispl.close();
    std::cout<<Absnomefile.c_str()<<" Saved "<<std::endl;
      
    
  }
	    
}
