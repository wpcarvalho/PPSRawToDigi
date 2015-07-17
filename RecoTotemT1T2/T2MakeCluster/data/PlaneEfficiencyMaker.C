/*
#include "TProfile.h"
#include <stdio.h>
#include <TFile.h>
#include <string.h>
#include <stddef.h>
#include <iostream>

#include <fstream>
#include "TH2.h"
#include "TPaletteAxis.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"

*/

#include "Getline.h"

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
#include "TH2.h"
#include "TPaletteAxis.h"
#include "TCanvas.h"
#include "TStyle.h" 
#include "TAxis.h"

#include "TProfile.h"
/*
bool StripChannelExcluded(int ch, int iid){
  bool isExcluded=false;

  if(iid==0)
    if(((ch>=52)&&(ch<=60))||(ch<=3)||((ch>=106)&&(ch<=110))||(ch>=120))
      isExcluded=true;

  if(iid==1)
    if((ch<=3)||((ch>=56)&&(ch<=60)))
      isExcluded=true;

  if(iid==15)
    if((ch>=66)&&(ch<=72))
       isExcluded=true;
      
  if(iid==16)
    if((ch>=16)&&(ch<=20))
      isExcluded=true;
  
  return isExcluded
}
*/

void PlaneEfficiencyMaker(std::string inputfilename, int quarterSel)
{
  std::cout<<"------------------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<"------------------------------------------------------------------------------------------------"<<std::endl;

  std::cout<<"Strategy For Strip DEAD CH Finder: Lower than 20% of the Y average in the Cumulative. Average Multiplicity bigger than 10 otherwise flag the Vfat as problematic. Low occupancy above 10 should be reproduced by normal inefficiency"<<std::endl;
  std::cout<<"------------------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<"------------------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<"Strategy For Pad DEAD CH  Finder: Lower than 20% of the Y average in the Cumulative. Average Multiplicity bigger than 2.5 otherwise flag the Vfat as problematic"<<std::endl;
  std::cout<<"------------------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<"------------------------------------------------------------------------------------------------"<<std::endl; 


  std::cout<<"Strategy For Strip Noisy CH Finder: Lower than 80% of the Y average in the Cumulative. Average Multiplicity lower than 16 otherwise flag the Vfat as problematic. Low occupancy above 10 should be reproduced by normal inefficiency"<<std::endl;
  std::cout<<"------------------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<"------------------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<"Strategy For Pad Noisy CH Finder: Lower than 80% of the Y average in the Cumulative. Average Multiplicity lower than 8 otherwise flag the Vfat as problematic"<<std::endl;
  std::cout<<"------------------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<"------------------------------------------------------------------------------------------------"<<std::endl; 


  int maxPadMultToBeNOTCompletelyNoisy=8;
  int minPadMultToBeNOTCompletelyDEAD=2.1;

  std::vector<std::pair<int,int> > NoisyVfatForcedToCalculate_PlIId;
  /*
  NoisyVfatForcedToCalculate_PlIId.push_back(std::pair<int,int>(19,10));
  NoisyVfatForcedToCalculate_PlIId.push_back(std::pair<int,int>(19,9));
  NoisyVfatForcedToCalculate_PlIId.push_back(std::pair<int,int>(19,11));
  NoisyVfatForcedToCalculate_PlIId.push_back(std::pair<int,int>(19,12));
  NoisyVfatForcedToCalculate_PlIId.push_back(std::pair<int,int>(19,13));

  NoisyVfatForcedToCalculate_PlIId.push_back(std::pair<int,int>(11,9));
  NoisyVfatForcedToCalculate_PlIId.push_back(std::pair<int,int>(11,10));
  NoisyVfatForcedToCalculate_PlIId.push_back(std::pair<int,int>(11,11));
  NoisyVfatForcedToCalculate_PlIId.push_back(std::pair<int,int>(11,12));
  NoisyVfatForcedToCalculate_PlIId.push_back(std::pair<int,int>(11,13));
  */

  std::vector<std::pair<int,int> > IgnoredVfat_PlIId;
  
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,0));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,1));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,2));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,3));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,4));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,5));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,6));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,7));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,8));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,9));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,10));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,11));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,12));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,13));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,14));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(20,15));  
  
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,0));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,1));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,2));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,3));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,4));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,5));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,6));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,7));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,8));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,9));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,10));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,11));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,12));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,13));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,14));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(21,15));  
  
  
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,0));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,1));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,2));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,3));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,4));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,5));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,6));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,7));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,8));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,9));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,10));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,11));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,12));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,13));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,14));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(22,15));  
  
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,0));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,1));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,2));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,3));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,4));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,5));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,6));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,7));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,8));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,9));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,10));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,11));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,12));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,13));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,14));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(23,15));  
  

 
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,0));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,1));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,2));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,3));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,4));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,5));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,6));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,7));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,8));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,9));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,10));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,11));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,12));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,13));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,14));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(24,15));  
  
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,0));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,1));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,2));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,3));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,4));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,5));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,6));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,7));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,8));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,9));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,10));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,11));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,12));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,13));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,14));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(25,15));  


 
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,0));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,1));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,2));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,3));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,4));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,5));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,6));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,7));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,8));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,9));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,10));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,11));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,12));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,13));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,14));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(26,15));  
  
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,0));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,1));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,2));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,3));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,4));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,5));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,6));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,7));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,8));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,9));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,10));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,11));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,12));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,13));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,14));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(27,15));  

  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,0));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,1));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,2));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,3));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,4));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,5));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,6));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,7));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,8));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,9));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,10));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,11));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,12));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,13));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,14));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(28,15));  
  
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,0));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,1));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,2));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,3));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,4));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,5));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,6));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,7));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,8));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,9));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,10));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,11));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,12));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,13));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,14));
  IgnoredVfat_PlIId.push_back(std::pair<int,int>(29,15));  




  //May_PHYS_RUNS_Tracks_NoAlign_H0123_5525-5535.root

  //H0_4384.root H1_4384.root H2_4384.root H3_4384.root
  // GeometricalEffi.root
  //inputfilename="H0EG.root";

  string outputfilename="OutEFFG.root";
  outputfilename.erase (outputfilename.end()-5, outputfilename.end());
  outputfilename=outputfilename+".dat";
  
  std::cout<<"Base outputFileName:"<<outputfilename.c_str()<<std::endl;

  TH1F vFatCumulative[40][17]; 
  TH1F vFatMultiplicty[40][17];
    //TProfile VFATOccup[40][17];
  //string namefile="~/SL/WorkingArea/Analysis/TrkDevelop/testDigi/SimuRootFiles/TriggerDigi/H0All.root";

  //string namefile="~/SL/WorkingArea/Analysis/DNDeta/Unfolding/root_file/Data/EfficiencyFiles/H1_3706.root";
  string namefile=inputfilename;

  //H0MuonTuneG.root
   //"~/SL/WorkingArea/Analysis/TrkDevelop/testDigi/SimuRootFiles/TriggerDigi/H0SingleMuon_50GeV_reco4380.root";

  TFile *f1 = new TFile(namefile.c_str());  
  int SelectedHalf=0;

  string str2="H0";
  size_t found=namefile.find(str2);
  if (found!=string::npos)
    SelectedHalf=0;
    
  str2="H1";
  size_t found1=namefile.find(str2);
  if (found1!=string::npos)
    SelectedHalf=1;
  
  str2="H2";
  size_t found2=namefile.find(str2);
  if (found2!=string::npos)
    SelectedHalf=2;
  
  str2="H3";
  size_t found3=namefile.find(str2);
  if (found3!=string::npos)
    SelectedHalf=3;
  
  SelectedHalf=quarterSel;
   

  TCanvas *a3 = new TCanvas("shapescaleh","shapescaleh",800,600);
  a3->Divide(1,1);

 
  TH2D* HistoMulti;
  HistoMulti=(new TH2D("HistoMulti","Vfat Average Multiplicity Map",40,-0.5,39.5,17,-0.5,16.5));
  HistoMulti->GetXaxis()->SetTitle("Plane") ;
  HistoMulti->GetYaxis()->SetTitle("VFAT position") ;

  
  TH2D* HistoOccup;
  HistoOccup=(new TH2D("HistoOccup","Vfat Occupancy Map",40,-0.5,39.5,17,-0.5,16.5));
  HistoOccup->GetXaxis()->SetTitle("Plane") ;
  HistoOccup->GetYaxis()->SetTitle("VFAT position") ;

  TH2D* NoisyHits;
  NoisyHits=(new TH2D("NoisyHits","Number of Noisy channel per VFAT ",40,-0.5,39.5,17,-0.5,16.5));
  NoisyHits->GetXaxis()->SetTitle("Plane") ;
  NoisyHits->GetYaxis()->SetTitle("VFAT position") ;

  TH2D* DeadChannels;
  DeadChannels=(new TH2D("DeadChannels","Number of dead Channels per VFAT",40,-0.5,39.5,17,-0.5,16.5));
  DeadChannels->GetXaxis()->SetTitle("Plane") ;
  DeadChannels->GetYaxis()->SetTitle("VFAT position") ;

  TH2D* ProblematicDeadsVFAT;
  ProblematicDeadsVFAT=(new TH2D("ProblematicDeadsVFAT","Problematic Low Occupancy VFATs",40,-0.5,39.5,17,-0.5,16.5));
  ProblematicDeadsVFAT->GetXaxis()->SetTitle("Plane") ;
  ProblematicDeadsVFAT->GetYaxis()->SetTitle("VFAT position") ;


  
  TH2D* NoisyChannels;
  NoisyChannels=(new TH2D("NoisyChannels","Number of Noisy Channels per VFAT",40,-0.5,39.5,17,-0.5,16.5));
  NoisyChannels->GetXaxis()->SetTitle("Plane") ;
  NoisyChannels->GetYaxis()->SetTitle("VFAT position") ;

  
  TH2D* ProblematicNoisyVFAT;
  ProblematicNoisyVFAT=(new TH2D("ProblematicNoisyVFAT","Problematic Noisy VFATs",40,-0.5,39.5,17,-0.5,16.5));
  ProblematicNoisyVFAT->GetXaxis()->SetTitle("Plane") ;
  ProblematicNoisyVFAT->GetYaxis()->SetTitle("VFAT position") ;


  //Good vfats plane 2,H0
   
  TH2D* RelativeEfficiencyFromOccup;
  RelativeEfficiencyFromOccup=(new TH2D("RelativeEfficiencyFromOccup","Relative Efficiency From Occupancy (respect to pl 2) - Abs diff",40,-0.5,39.5,17,-0.5,16.5));
  RelativeEfficiencyFromOccup->GetXaxis()->SetTitle("Plane") ;
  RelativeEfficiencyFromOccup->GetYaxis()->SetTitle("VFAT position") ;


  TH2D* RelativeEfficiencyFromOccup2;
  RelativeEfficiencyFromOccup2=(new TH2D("RelativeEfficiencyFromOccup2","Relative Efficiency From Occupancy (respect to pl 2) ",40,-0.5,39.5,17,-0.5,16.5));
  RelativeEfficiencyFromOccup2->GetXaxis()->SetTitle("Plane") ;
  RelativeEfficiencyFromOccup2->GetYaxis()->SetTitle("VFAT position") ;

 f1->cd("EffiGeometry");

 char sZname2[1024];
 char sZnamehist[1024];
std::cout<<"HereA"<<std::endl;
 

//std::ofstream myfile2;
//std::string outputfilenamemulti="Multipl_"+outputfilename;
////myfile2.open (outputfilenamemulti.c_str());

 TH2D EffiStripN[40]; 
 TH2D EffiPadN[40];
 TH2D EffiStripD[40]; 
 TH2D EffiPadD[40];

 TH2D EffiStrip[40]; 
 TH2D EffiPad[40];

 for(unsigned int m=0;m<40; m++)
   // for(unsigned int k=0;k<17; k++)
   {
     
     sprintf(sZname2, "HGeometry_StripEfficiency_Num %d", m); 
     sprintf(sZnamehist, "HGeometry_StripEfficiency_Num %d", m);
     TH2D  *h0 = (TH2D*)gDirectory->Get(sZname2);
     EffiStripN[m]=(*h0);
     //HGeometry_StripEfficiency_Num[m]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,numRsectEffi,-0.5,numRsectEffi-0.5,numPhisectEffi,-0.5,numPhisectEffi-0.5));

     
     sprintf(sZname2, "HGeometry_StripEfficiency_Den %d", m); 
     sprintf(sZnamehist, "HGeometry_StripEfficiency_Den %d", m);
     TH2D  *h1 = (TH2D*)gDirectory->Get(sZname2);
     EffiStripD[m]=(*h1);
     //HGeometry_StripEfficiency_Den[m]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,numRsectEffi,-0.5,numRsectEffi-0.5,numPhisectEffi,-0.5,numPhisectEffi-0.5)); 

     
     sprintf(sZname2, "HGeometry_PadEfficiency_Num %d", m); 
     sprintf(sZnamehist, "HGeometry_PadEfficiency_Num %d", m);
     TH2D  *h2 = (TH2D*)gDirectory->Get(sZname2);
     EffiPadN[m]=(*h2);
     // HGeometry_PadEfficiency_Num[m]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,numRsectEffi,-0.5,numRsectEffi-0.5,numPhisectEffi,-0.5,numPhisectEffi-0.5));

     
     sprintf(sZname2, "HGeometry_PadEfficiency_Den %d", m); 
     sprintf(sZnamehist, "HGeometry_PadEfficiency_Den %d", m);
     TH2D  *h3 = (TH2D*)gDirectory->Get(sZname2);     
     EffiPadD[m]=(*h3);
     //HGeometry_PadEfficiency_Den[m]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,numRsectEffi,-0.5,numRsectEffi-0.5,numPhisectEffi,-0.5,numPhisectEffi-0.5)); 
     
   }
  int numxbin=EffiPadD[0].GetXaxis()->GetNbins();
  int numybin=EffiPadD[0].GetYaxis()->GetNbins();
  int numRsectEffi=numxbin;
  int numPhisectEffi=numybin;

 for(unsigned int m=0;m<40; m++){
   
   sprintf(sZname2, "HGeometry_PadEfficiency %d", m); 
   sprintf(sZnamehist, "Pad Efficiency Plane %d", m);
   TH2D  *EffiPad1=new TH2D(sZname2,sZnamehist,numRsectEffi,-0.5,numRsectEffi-0.5,numPhisectEffi,-0.5,numPhisectEffi-0.5);
   EffiPad1->GetXaxis()->SetTitle("#eta Sector");
   EffiPad1->GetYaxis()->SetTitle("#phi Sector");
 
   EffiPad[m]=(*EffiPad1);
 
   sprintf(sZname2, "HGeometry_StripEfficiency %d", m); 
   sprintf(sZnamehist, "Strip Efficiency Plane %d", m);
   TH2D  *EffiStrip1=new TH2D(sZname2,sZnamehist,numRsectEffi,-0.5,numRsectEffi-0.5,numPhisectEffi,-0.5,numPhisectEffi-0.5);
   EffiStrip1->GetXaxis()->SetTitle("#eta Sector");
   EffiStrip1->GetYaxis()->SetTitle("#phi Sector");
   EffiStrip[m]=(*EffiStrip1);
}


 std::cout<<"HereB"<<std::endl;

 double avgStripMult=0;
 double avgPadMult=0; 
 unsigned int countgdpad=0;
 unsigned int countgdstrip=0;

 

 
 for(unsigned int m=0;m<40; m++)
   {
     int planeID=m;
     //if((m==plane1)||(m==plane2))
     for(unsigned int xx=0;xx<numxbin;xx++)
       for(unsigned int yy=0;yy<numybin;yy++)
	 {
	   double strN=EffiStripN[m].GetBinContent(xx+1,yy+1); 
	   double padN=EffiPadN[m].GetBinContent(xx+1,yy+1);
	   double strD=EffiStripD[m].GetBinContent(xx+1,yy+1); 
	   double padD=EffiPadD[m].GetBinContent(xx+1,yy+1);
	   
	   double effi=0.;
	   
	   if(strD>0)
	     effi= strN/strD;

	   effi=effi*1000;//In order to improve precision.

	   EffiStrip[m].Fill(xx,yy,effi);
	   
	   effi=0.;
	   if(padD>0)
	     effi= padN/padD;

	   effi=effi*1000;//In order to improve precision.	   
	   EffiPad[m].Fill(xx,yy,effi);

	 }
     
   }


 
 // myfile2.close();
 TFile *MyFile; 
 MyFile= new TFile("EffiG.root","UPDATE");
 if(MyFile->IsOpen()) 
   printf("File opened successfully\n");
 else
   MyFile= new TFile("EffiG.root","RECREATE");
				
for(unsigned int i=0;i<40;i++)
   if(i<SelectedHalf*10+10)
     if(i>=SelectedHalf*10){
        std::cout<<"HereC"<<i<<std::endl;
       EffiPad[i].SetDirectory(0);
       EffiPad[i].Write();
       EffiStrip[i].SetDirectory(0);
       EffiStrip[i].Write();
     }
 MyFile->Close();
 std::cout<<"HereC"<<std::endl;

 /*
 for(unsigned int i=0;i<40;i++)
   if(i<SelectedHalf*10+10)
     if(i>=SelectedHalf*10)
     {
       TCanvas *plotter;
       if(i%4 == 0){
	 
	 sprintf(sZname2," canvas PAD i %d", i); 
	 
	 plotter= new TCanvas(sZname2,sZname2,800,600);
	 plotter->Divide(2,2);
	 
       }
      
       int posincanv=(i%10)%4 + 1;
       std::cout<<"HereD "<<i<<" "<<posincanv<<std::endl;
       
       plotter->cd(posincanv);
       std::cout<<"HereD "<<i<<std::endl;
       
       TH2D * histo = new TH2D(EffiPad[i]);
       std::cout<<"HereD "<<i<<std::endl;
       histo->Scale(1.0/(1000));
       histo->SetDirectory(0);
       histo->Draw(("COLZ"));
       
     }






for(unsigned int i=0;i<40;i++)
   if(i<SelectedHalf*10+10)
     if(i>=SelectedHalf*10)
     {
       TCanvas *plotter;
       if(i%4 == 0){
	 
	 sprintf(sZname2," canvas STRIP i %d", i); 
	 
	 plotter= new TCanvas(sZname2,sZname2,800,600);
	 plotter->Divide(2,2);
	 
       }
       
       int posincanv=i%4 + 1;
       plotter->cd(posincanv);
       TH2D * histo = new TH2D(EffiStrip[i]);
       histo->Scale(1.0/(1000));
       histo->SetDirectory(0);
       histo->Draw(("COLZ"));
       //vFatCumulative[vfatToPlot.at(i).first][vfatToPlot.at(i).second].Draw("A");
       //std::cout<<vfatToPlot.at(i).first<<" "<<vfatToPlot.at(i).second<<" "<<vFatCumulative[vfatToPlot.at(i).first][vfatToPlot.at(i).second].GetEntries()<<std::endl;
     }


 */

 // HistoOccup->Fill(row,col,weight);
std::cout<<"Here!!"<<std::endl;
/*

 
 
 HistoOccup->Draw("COLZ");
 TCanvas *a4 = new TCanvas("sh4apescaleh","sh4apescaleh",800,600);
 a4->Divide(1,1);
 HistoMulti->Draw("COLZ");

 TCanvas *aDeadChannels = new TCanvas("aDeadChannels","aDeadChannels",800,600);
 aDeadChannels->Divide(1,1);
 DeadChannels->Draw("COLZ");

 TCanvas *aProblematicDeadsVFAT = new TCanvas("aProblematicDeadsVFAT","aProblematicDeadsVFAT",800,600);
 aProblematicDeadsVFAT->Divide(1,1);
 ProblematicDeadsVFAT->Draw("COLZ");


 TCanvas *aNoisyChannels = new TCanvas("aNoisyChannels","aNoisyChannels",800,600);
 aNoisyChannels->Divide(1,1);
 NoisyChannels->Draw("COLZ");

 TCanvas *aProblematicNoisysVFAT = new TCanvas("aProblematicNoisysVFAT","aProblematicNoisysVFAT",800,600);
 aProblematicNoisysVFAT->Divide(1,1);
 ProblematicNoisyVFAT->Draw("COLZ");

 TCanvas *aRelativeEfficiencyFromOccup= new TCanvas("RelativeEfficiencyFromOccup","RelativeEfficiencyFromOccup",800,600);
 aRelativeEfficiencyFromOccup->Divide(1,1);
 RelativeEfficiencyFromOccup->Draw("COLZ");

 //TCanvas *aRelativeEfficiencyFromOccup2= new TCanvas("RelativeEfficiencyFromOccup2","RelativeEfficiencyFromOccup2",800,600);
 //aRelativeEfficiencyFromOccup2->Divide(1,1);
 //RelativeEfficiencyFromOccup2->Draw("COLZ");


 TCanvas *plotBest= new TCanvas("plotBest","plotBest",800,600);
 plotBest->Divide(2,2);
 plotBest->cd(1);
 StripBestProfile0->Draw();
 plotBest->cd(2);
 StripBestProfile1->Draw();
 plotBest->cd(3);
 StripBestProfile15->Draw();
 plotBest->cd(4);
 StripBestProfile16->Draw();

 TCanvas *plotBestPad= new TCanvas("plotBestPad","plotBestPad",800,600);
 plotBestPad->Divide(2,2);
 plotBestPad->cd(1);
 PadBestProfile2->Draw();
 plotBestPad->cd(2);
 PadBestProfile14->Draw();
 plotBestPad->cd(3);
 PadBestProfileCentrals->Draw();
 

 //Vfat to plots
 std::vector<std::pair<int,int> > vfatToPlot;
 vfatToPlot.push_back(std::pair<int,int>(3,14));
 vfatToPlot.push_back(std::pair<int,int>(19,14));
 vfatToPlot.push_back(std::pair<int,int>(29,14));
 vfatToPlot.push_back(std::pair<int,int>(39,14));

 vfatToPlot.push_back(std::pair<int,int>(7,16));
 vfatToPlot.push_back(std::pair<int,int>(17,16));
 vfatToPlot.push_back(std::pair<int,int>(34,16));

 vfatToPlot.push_back(std::pair<int,int>(7,1));
 vfatToPlot.push_back(std::pair<int,int>(28,1));

 vfatToPlot.push_back(std::pair<int,int>(34,14));
 vfatToPlot.push_back(std::pair<int,int>(38,14));
 vfatToPlot.push_back(std::pair<int,int>(24,14));

 vfatToPlot.push_back(std::pair<int,int>(36,0));
 vfatToPlot.push_back(std::pair<int,int>(36,1));
 vfatToPlot.push_back(std::pair<int,int>(36,15));
 vfatToPlot.push_back(std::pair<int,int>(36,16));


 for(unsigned int i=0;NoisyVFATByEyes.size();i++){
   std::cout<<"CONFIRMED Noisy Plane: "<<NoisyVFATByEyes.at(i).first<<"  VFAT: "<<NoisyVFATByEyes.at(i).second<<std::endl;
 }
 
 for(unsigned int i=0;i<vfatToPlot.size();i++){
   TCanvas *plotter;
   if(i%4 == 0){
      sprintf(sZnamehiste," canvas i %d", i);
      
      plotter= new TCanvas(sZnamehiste,sZnamehiste,800,600);
      plotter->Divide(2,2);
      
   }

   int posincanv=i%4 + 1;
   plotter->cd(posincanv);
   TH1F * histo = new TH1F(vFatCumulative[vfatToPlot.at(i).first][vfatToPlot.at(i).second]);
   histo->SetDirectory(0);
   histo->Draw();
   //vFatCumulative[vfatToPlot.at(i).first][vfatToPlot.at(i).second].Draw("A");
   std::cout<<vfatToPlot.at(i).first<<" "<<vfatToPlot.at(i).second<<" "<<vFatCumulative[vfatToPlot.at(i).first][vfatToPlot.at(i).second].GetEntries()<<std::endl;
 }
 */
 
}
//VFATOccupNormalized
	 //file_op<<"DetId: "<<planeID<<"  VFat_iid: "<<iid<<"   Eff: "<< VFATOccup[m]->GetBinContent(iid+1)<<"        pm: "<< VFATOccup[m]->GetBinError(iid+1)<<" ( stat: "<<VFATOccup[m]->GetBinEntries(iid+1)<<" ) \n";
