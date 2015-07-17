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
double XYTollerance=2.0; //Mechanical tollerance in mm
double TILTTollerance=0.006; //Mechanical tollerance in mm
//////////////////////////////////////////////////////

void InternalMisalGeneratorGLOB(int N_scenarios, std::string outputFolder);




void InternalMisalGeneratorGLOB(int N_scenarios, std::string outputFolder){



  for(unsigned int i=0;i<=N_scenarios;i++){
    
    ofstream outputFileDispl;
    stringstream st_I;//create a stringstream
    st_I << (i);//add number to the stream
    
    std::string Absnomefile= outputFolder+"/RandomGlobalMisal_Scenario"+st_I.str()+".dat";
    

    outputFileDispl.open(Absnomefile.c_str());

    
    TRandom *r0 = new TRandom();
    r0->SetSeed(i);

    for(unsigned int det=0;det<4;det++){

     
       
      double x = r0->Uniform(-(1.0)*XYTollerance,XYTollerance);
      double y = r0->Uniform(-(1.0)*XYTollerance,XYTollerance);
      double tx = r0->Uniform(-(1.0)*TILTTollerance,TILTTollerance);
      double ty = r0->Uniform(-(1.0)*TILTTollerance,TILTTollerance);
      
     

      outputFileDispl<<"QuarterId: "<<det<<"  TiltX: "<<tx<<"  TiltY: "<<ty<<"  ShiftX0: "<<x<<"  ShiftY0:     "<<y<<"\n";

    }
    

    outputFileDispl.close();
    std::cout<<Absnomefile.c_str()<<" Saved "<<std::endl;
      
    
  }
	    
}
/*
QuarterId:  0  TiltX:  0.0033  TiltY:  -0.002  ShiftX0: -2.0   ShiftY0:  -2.0
QuarterId:  1  TiltX:  0.001  TiltY:  -0.0067 ShiftX0:  2.0   ShiftY0:  0.0
QuarterId:  2  TiltX:  -0.0052  TiltY:  -0.0005 ShiftX0:  0.0   ShiftY0:  1.0
QuarterId:  3  TiltX:  -0.0015  TiltY:  -0.0005  ShiftX0:  -1.0   ShiftY0:  -1.0
*/
