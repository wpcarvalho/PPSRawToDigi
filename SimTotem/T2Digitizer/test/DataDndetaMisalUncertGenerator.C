

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
double TILTTollerance=0.0004; //Mechanical tollerance in mm
//////////////////////////////////////////////////////

void DataDndetaMisalUncertGenerator(std::string outputFolder);




void DataDndetaMisalUncertGenerator(std::string outputFolder){


  double startingaddtilt=-0.0004;
  double increm=0.0001;
  double tx=startingaddtilt; double ty=startingaddtilt;
  for(unsigned int i=0;i<=8;i++){
    tx=startingaddtilt+i*increm;
  for(unsigned int j=0;j<=8;j++){

    ty=startingaddtilt+j*increm;
    ofstream outputFileDispl;
    stringstream st_I;//create a stringstream
    st_I << (i*9+j);//add number to the stream PythiaGlobalMisal_Scenario
    
    std::string Absnomefile= outputFolder+"/Data2012GlobalMisal_ScenarioRec"+st_I.str()+".dat";
    TRandom *r0 = new TRandom();
    r0->SetSeed(i);
    double intalres=0.00;
    outputFileDispl.open(Absnomefile.c_str());
    double pert=0.;double pert2=0.;
      
    outputFileDispl<<"DetId: 0   DX: 0   DY: 0   EDX: 0.00628654   EDY: 0.00519129"<<"\n";
    outputFileDispl<<"DetId: 1   DX: 0.109883   DY: 0.139714   EDX: 0.00688938   EDY: 0.00565935"<<"\n";
    outputFileDispl<<"DetId: 2   DX: 0.249521   DY: -0.060385   EDX: 0.00530346   EDY: 0.00440499"<<"\n";
    outputFileDispl<<"DetId: 3   DX: -0.181116   DY: 0.0879289   EDX: 0.00528719   EDY: 0.00429201"<<"\n";
    outputFileDispl<<"DetId: 4   DX: 0.0913549   DY: -0.112072   EDX: 0.00525001   EDY: 0.00426748"<<"\n";
    outputFileDispl<<"DetId: 5   DX: 0.0538705   DY: 0.00251201   EDX: 0.00546854   EDY: 0.00428894"<<"\n";
    outputFileDispl<<"DetId: 6   DX: 0.11467   DY: -0.26168   EDX: 0.00534813   EDY: 0.00438869"<<"\n";
    outputFileDispl<<"DetId: 7   DX: -0.161602   DY: -0.00109845   EDX: 0.00543115   EDY: 0.00443948"<<"\n";
    outputFileDispl<<"DetId: 8   DX: 0.105276   DY: -0.0888045   EDX: 0.00610674   EDY: 0.00491889"<<"\n";
    outputFileDispl<<"DetId: 9   DX: -0.13468   DY: 0.42308   EDX: 0.00630319   EDY: 0.00524005"<<"\n";
    outputFileDispl<<"DetId: 10   DX: -0.173021   DY: 0.0215625   EDX: 0.00625044   EDY: 0.00512344"<<"\n";
    outputFileDispl<<"DetId: 11   DX: -0.0268242   DY: 0.288521   EDX: 0.00613805   EDY: 0.00497029"<<"\n";
    outputFileDispl<<"DetId: 12   DX: 0.0510472   DY: -0.0947231   EDX: 0.00567855   EDY: 0.0045147"<<"\n";
    outputFileDispl<<"DetId: 13   DX: 0.0368867   DY: 0.231111   EDX: 0.00728036   EDY: 0.00528245"<<"\n";
    outputFileDispl<<"DetId: 14   DX: 0.0439104   DY: -0.332416   EDX: 0.00596278   EDY: 0.00462728"<<"\n";
    outputFileDispl<<"DetId: 15   DX: 0.0278239   DY: -0.0676228   EDX: 0.0056851   EDY: 0.00442306"<<"\n";
    outputFileDispl<<"DetId: 16   DX: 0.00472105   DY: -0.243941   EDX: 0.00557882   EDY: 0.00444284"<<"\n";
    outputFileDispl<<"DetId: 17   DX: 0.122297   DY: 0.0549978   EDX: 0.00566709   EDY: 0.00447238"<<"\n";
    outputFileDispl<<"DetId: 18   DX: -0.0455546   DY: 0.0840613   EDX: 0.006113   EDY: 0.00490237"<<"\n";
    outputFileDispl<<"DetId: 19   DX: -0.158152   DY: 0.309336   EDX: 0.00646846   EDY: 0.00521199"<<"\n";
    outputFileDispl<<"DetId: 20   DX: 0.395829   DY: 0.112834   EDX: 0.0193911   EDY: 0.00849627"<<"\n";
    outputFileDispl<<"DetId: 21   DX: 0.439858   DY: 0.238498   EDX: 0.0101601   EDY: 0.00812125"<<"\n";
    outputFileDispl<<"DetId: 22   DX: 0.00651766   DY: -0.236361   EDX: 0.00866745   EDY: 0.00722846"<<"\n";
    outputFileDispl<<"DetId: 23   DX: -0.178073   DY: -0.0194215   EDX: 0.00852637   EDY: 0.00695195"<<"\n";
    outputFileDispl<<"DetId: 24   DX: -0.279815   DY: 0.115125   EDX: 0.00820661   EDY: 0.00669243"<<"\n";
    outputFileDispl<<"DetId: 25   DX: -0.200654   DY: -0.147422   EDX: 0.00829252   EDY: 0.00691366"<<"\n";
    outputFileDispl<<"DetId: 26   DX: 0.154162   DY: -0.238579   EDX: 0.00810406   EDY: 0.00688069"<<"\n";
    outputFileDispl<<"DetId: 27   DX: -0.0639781   DY: 0.15756   EDX: 0.00826867   EDY: 0.00703872"<<"\n";
    outputFileDispl<<"DetId: 28   DX: 0.0338314   DY: 0.0436531   EDX: 0.00937458   EDY: 0.00812184"<<"\n";
    outputFileDispl<<"DetId: 29   DX: 0.332508   DY: 0.184904   EDX: 0.00974124   EDY: 0.010562"<<"\n";
    outputFileDispl<<"DetId: 30   DX: 0.23356   DY: 0.197356   EDX: 0.0118458   EDY: 0.0158145"<<"\n";
    outputFileDispl<<"DetId: 31   DX: -0.0454532   DY: -0.133409   EDX: 0.0108498   EDY: 0.0161747"<<"\n";
    outputFileDispl<<"DetId: 32   DX: -0.161353   DY: -0.0969199   EDX: 0.0102276   EDY: 0.0126811"<<"\n";
    outputFileDispl<<"DetId: 33   DX: -0.030238   DY: 0.0194859   EDX: 0.010612   EDY: 0.0123552"<<"\n";
    outputFileDispl<<"DetId: 34   DX: -0.0115858   DY: 0.175983   EDX: 0.00980165   EDY: 0.0121778"<<"\n";
    outputFileDispl<<"DetId: 35   DX: 0   DY: 0   EDX: 0   EDY: 0"<<"\n";
    outputFileDispl<<"DetId: 36   DX: -0.0180024   DY: -0.210479   EDX: 0.0101133   EDY: 0.0131674"<<"\n";
    outputFileDispl<<"DetId: 37   DX: 0.0189391   DY: -0.398533   EDX: 0.0106669   EDY: 0.0132965"<<"\n";
    outputFileDispl<<"DetId: 38   DX: 0.0444098   DY: 0.492166   EDX: 0.0117324   EDY: 0.0144723"<<"\n";
    outputFileDispl<<"DetId: 39   DX: -0.0293981   DY: 0.0449121   EDX: 0.0112256   EDY: 0.0143539"<<"\n";


    outputFileDispl<<"QuarterId:  0  TiltX:  "<<0.0035+tx<<"  TiltY:  "<<-0.0013+ty<<"    ShiftX0: -1.52   ShiftY0:  -0.62"<<"\n";
    outputFileDispl<<"QuarterId:  1  TiltX:  "<<0.0012+tx<<"  TiltY:  "<<-0.0065+ty<<"   ShiftX0:  -2.34   ShiftY0:  -0.80"<<"\n";
    outputFileDispl<<"QuarterId:  2  TiltX:  "<<-0.0082+tx<<"  TiltY:  "<<-0.0007+ty<<"    ShiftX0:  0.13   ShiftY0: -1.30"<<"\n";
    outputFileDispl<<"QuarterId:  3  TiltX:  "<<-0.0030+tx<<"   TiltY:  "<<0.0026+ty<<"  ShiftX0:  1.40   ShiftY0:  -1.80"<<"\n";


    
    

    outputFileDispl.close();
    std::cout<<Absnomefile.c_str()<<" Saved "<<std::endl;
      
    
  }
}
	    
}
/*
QuarterId:  0  TiltX:  0.0033  TiltY:  -0.002  ShiftX0: -2.0   ShiftY0:  -2.0
QuarterId:  1  TiltX:  0.001  TiltY:  -0.0067 ShiftX0:  2.0   ShiftY0:  0.0
QuarterId:  2  TiltX:  -0.0052  TiltY:  -0.0005 ShiftX0:  0.0   ShiftY0:  1.0
QuarterId:  3  TiltX:  -0.0015  TiltY:  -0.0005  ShiftX0:  -1.0   ShiftY0:  -1.0
*/
