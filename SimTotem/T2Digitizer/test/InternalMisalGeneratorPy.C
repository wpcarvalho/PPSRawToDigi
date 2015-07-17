

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

void InternalMisalGeneratorPy(int N_scenarios, std::string outputFolder);




void InternalMisalGeneratorPy(int N_scenarios, std::string outputFolder){



  for(unsigned int i=0;i<=N_scenarios;i++){
    
    ofstream outputFileDispl;
    stringstream st_I;//create a stringstream
    st_I << (i);//add number to the stream PythiaGlobalMisal_Scenario
    
    std::string Absnomefile= outputFolder+"/PythiaGlobalMisal_ScenarioRec"+st_I.str()+".dat";
    TRandom *r0 = new TRandom();
    r0->SetSeed(i);
    double intalres=0.00;
    outputFileDispl.open(Absnomefile.c_str());
    double pert=0.;double pert2=0.;
    pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);    
    outputFileDispl<<"DetId: 0   DX: "<<pert+ -0.094   <<"    DY: "<<-0.00130582   +pert2 << "    EDX: 0.0177629   EDY: 0.0146356  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 1   DX: "<<pert+ 0.0544365   <<"    DY: "<<0.118436   +pert2 << "    EDX: 0.0183491   EDY: 0.0149959  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 2   DX: "<<pert+ 0.180176   <<"    DY: "<<-0.0729026   +pert2 << "    EDX: 0.0157517   EDY: 0.0125472  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 3   DX: "<<pert+ -0.176028   <<"    DY: "<<0.0768611   +pert2 << "    EDX: 0.0158455   EDY: 0.012052  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 4   DX: "<<pert+ 0.0720577   <<"    DY: "<<-0.101375   +pert2 << "    EDX: 0.0150138   EDY: 0.0123267  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 5   DX: "<<pert+ 0.0169654   <<"    DY: "<<0.0302686   +pert2 << "    EDX: 0.01893   EDY: 0.0128677  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 6   DX: "<<pert+ 0.105   <<"    DY: "<<-0.214095   +pert2 << "    EDX: 0.0159336   EDY: 0.0134207  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 7   DX: "<<pert+ -0.137541   <<"    DY: "<<0.0407741   +pert2 << "    EDX: 0.0167705   EDY: 0.0133451  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 8   DX: "<<pert+ 0.120441   <<"    DY: "<<-0.0306582   +pert2 << "    EDX: 0.0175432   EDY: 0.0151516  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 9   DX: "<<pert+ -0.111904   <<"    DY: "<<0.412776   +pert2 << "    EDX: 0.0181478   EDY: 0.0212197  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 10   DX: "<<pert+ -0.176777   <<"    DY: "<<0.0158931   +pert2 << "    EDX: 0.0202218   EDY: 0.0235347  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 11   DX: "<<pert+ -0.0544126   <<"    DY: "<<0.273946   +pert2 << "    EDX: 0.0208537   EDY: 0.0231178  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 12   DX: "<<pert+ 0.027   <<"    DY: "<<-0.0641708   +pert2 << "    EDX: 0.0183085   EDY: 0.0193093  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 13   DX: "<<pert+ 0.0340103   <<"    DY: "<<0.193815   +pert2 << "    EDX: 0.0261153   EDY: 0.0268773  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 14   DX: "<<pert+ 0.0628499   <<"    DY: "<<-0.263689   +pert2 << "    EDX: 0.0184397   EDY: 0.0221075  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 15   DX: "<<pert+ 0.013539   <<"    DY: "<<-0.0488511   +pert2 << "    EDX: 0.0197318   EDY: 0.0192994 "<<"\n";  pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 16   DX: "<<pert+ -0.0157598   <<"    DY: "<<-0.227755   +pert2 << "    EDX: 0.0170816   EDY: 0.0212668  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 17   DX: "<<pert+ 0.117013   <<"    DY: "<<0.0555559   +pert2 << "    EDX: 0.019321   EDY: 0.021407  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 18   DX: "<<pert+ -0.057   <<"    DY: "<<0.0800558   +pert2 << "    EDX: 0.0199645   EDY: 0.0226092  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 19   DX: "<<pert+ -0.157459   <<"    DY: "<<0.289785   +pert2 << "    EDX: 0.0207546   EDY: 0.0256126 "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);   
    outputFileDispl<<"DetId: 20   DX: "<<pert+ 0.242043   <<"    DY: "<<0.118509   +pert2 << "    EDX: 0.0509099   EDY: 0.0238636  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 21   DX: "<<pert+ 0.271205   <<"    DY: "<<0.237057   +pert2 << "    EDX: 0.0398959   EDY: 0.0249258  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 22   DX: "<<pert+ -0.0421631   <<"    DY: "<<-0.213763   +pert2 << "    EDX: 0.0295261   EDY: 0.0211119  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 23   DX: "<<pert+ -0.102   <<"    DY: "<<-0.0552067   +pert2 << "    EDX: 0.0326146   EDY: 0.0202966  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 24   DX: "<<pert+ -0.238918   <<"    DY: "<<0.0860672   +pert2 << "    EDX: 0.0317614   EDY: 0.0253245 "<<"\n";  pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 25   DX: "<<pert+ -0.125347   <<"    DY: "<<-0.112258   +pert2 << "    EDX: 0.0315969   EDY: 0.0203659 "<<"\n";  pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 26   DX: "<<pert+ 0.113563   <<"    DY: "<<-0.180858   +pert2 << "    EDX: 0.0267729   EDY: 0.0232102 "<<"\n";  pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 27   DX: "<<pert+ -0.0704475   <<"    DY: "<<0.148889   +pert2 << "    EDX: 0.0266202   EDY: 0.0230308 "<<"\n";  pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 28   DX: "<<pert+ 0.0297   <<"    DY: "<<0.0730924   +pert2 << "    EDX: 0.0297419   EDY: 0.0257488 "<<"\n";  pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 29   DX: "<<pert+ 0.26683   <<"    DY: "<<0.174256   +pert2 << "    EDX: 0.0335964   EDY: 0.0367055  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 30   DX: "<<pert+ 0.266428   <<"    DY: "<<0.3171   +pert2 << "    EDX: 0.028405   EDY: 0.0672571 "<<"\n";  pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 31   DX: "<<pert+ -0.00534646   <<"    DY: "<<-0.119534   +pert2 << "    EDX: 0.0250299   EDY: 0.0474638  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 32   DX: "<<pert+ -0.1542   <<"    DY: "<<-0.101934   +pert2 << "    EDX: 0.0237587   EDY: 0.0456516 "<<"\n";  pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 33   DX: "<<pert+ -0.0267682   <<"    DY: "<<0.022554   +pert2 << "    EDX: 0.0231409   EDY: 0.0443086 "<<"\n";  pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 34   DX: "<<pert+ -0.00618162   <<"    DY: "<<0.120274   +pert2 << "    EDX: 0.0233663   EDY: 0.0453296 "<<"\n";  pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 35   DX: "<<pert+ 0   <<"    DY: "<<0   +pert2 << "    EDX: 0.   EDY: 0.  "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 36   DX: "<<pert+ -0.0204563   <<"    DY: "<<-0.129652   +pert2 << "    EDX: 0.0238055   EDY: 0.0542502 "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);   
    outputFileDispl<<"DetId: 37   DX: "<<pert+ 0.0036   <<"    DY: "<<-0.246896   +pert2 << "    EDX: 0.0290415   EDY: 0.0567559 "<<"\n";  pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 38   DX: "<<pert+ 0.0561515   <<"    DY: "<<0.348367   +pert2 << "    EDX: 0.0266533   EDY: 0.0859148 "<<"\n";  pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  
    outputFileDispl<<"DetId: 39   DX: "<<pert+ 0.0316245   <<"    DY: "<<0.108396   +pert2 << "    EDX: 0.0261034   EDY: 0.0517604 "<<"\n"; pert=r0->Uniform(-(1.0)*intalres,intalres);    pert2=r0->Uniform(-(1.0)*intalres,intalres);  

  

    for(unsigned int det=0;det<4;det++){

     
       
      double x = r0->Uniform(-(1.0)*XYTollerance,XYTollerance);
      double y = r0->Uniform(-(1.0)*XYTollerance,XYTollerance);
      double tx = r0->Uniform(-(1.0)*TILTTollerance,TILTTollerance);
      double ty = r0->Uniform(-(1.0)*TILTTollerance,TILTTollerance);
      
     

      //      outputFileDispl<<"QuarterId: "<<det<<"  TiltX: "<<tx<<"  TiltY: "<<ty<<"  ShiftX0: "<<x<<"  ShiftY0:     "<<y<<"\n";

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
