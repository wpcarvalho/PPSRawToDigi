//Author Mirko Berretti (mirko.berretti@gmail.com)
//University of Siena and Pisa INFN

//3_1_1_version
#include "TotemAnalysis/T2HalfQuarterTrkEfficiency/interface/T2HalfQuarterTrkEfficiency.h"
#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include "CLHEP/Random/RandGaussQ.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/T2Digi/interface/T2StripDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2StripDigi.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2PadDigi.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/T2Cluster/interface/cluster_entry.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "FWCore/Framework/interface/MakerMacros.h"



T2HalfQuarterTrkEfficiency::T2HalfQuarterTrkEfficiency(const edm::ParameterSet& iConfig){
  
  t2PadDigiCollectionLabel = iConfig.getParameter<edm::InputTag>("T2PadDigiCollectionLabel");
  t2StripDigiCollectionLabel = iConfig.getParameter<edm::InputTag>("T2StripDigiCollectionLabel");
  t2VfatInformationLabel = iConfig.getParameter<edm::InputTag>("T2VfatInformationLabel");
  verbosity=iConfig.getParameter<bool>("verbosity");
  RefTrkMultiplicity=iConfig.getParameter<unsigned int>("RefTrkMultiplicity");
  MaxPad=iConfig.getParameter<unsigned int>("MaxPad");
  MaxStrip=iConfig.getParameter<unsigned int>("MaxStrip");
  MaxDphi=iConfig.getParameter<double>("MaxDphi");
  maxdrhit=iConfig.getParameter<double>("maxdrhit");
  maxdphihit=iConfig.getParameter<double>("maxdphihit");
  MaxTrkInProcess=iConfig.getParameter<unsigned int>("MaxTrkInProcess");

  Refstring= iConfig.getParameter<std::string>("Refstring");
  outputFileName = iConfig.getUntrackedParameter<std::string>("OutputFile");
  CluLabel = iConfig.getParameter<std::string>("CluLabel");
  HitLabel = iConfig.getParameter<std::string>("HitLabel");
  RoadLabel = iConfig.getParameter<std::string>("RoadLabel");
  TrackLabel= iConfig.getParameter<std::string>("TrackLabel");
  
  HitLabelTestedQ= iConfig.getParameter<std::string>("HitLabelTestedQ");
  RefCluLabel= iConfig.getParameter<std::string>("RefCluLabel");
  
  RoadLabelH1 = iConfig.getParameter<std::string>("RoadLabelH1");
  RoadLabelH0 = iConfig.getParameter<std::string>("RoadLabelH0");
  

  TrackLabelFirstPlanes= iConfig.getParameter<std::string>("TrackLabelFirstPlanes");
  TrackLabelLastPlanes= iConfig.getParameter<std::string>("TrackLabelLastPlanes");
  PadCluLabelFirstPlanes= iConfig.getParameter<std::string>("PadCluLabelFirstPlanes"); 
  PadCluLabelLastPlanes= iConfig.getParameter<std::string>("PadCluLabelLastPlanes");


  RoadLabelH2 = iConfig.getParameter<std::string>("RoadLabelH2");
  TrackLabelH2= iConfig.getParameter<std::string>("TrackLabelH2");
  
  RoadLabelH3 = iConfig.getParameter<std::string>("RoadLabelH3");
  TrackLabelH3= iConfig.getParameter<std::string>("TrackLabelH3");
   
  Trk_PhiSeparation= iConfig.getParameter<double>("Trk_PhiSeparation");
  Trk_RSeparation= iConfig.getParameter<double>("Trk_RSeparation");
  Trk_EtaSeparation= iConfig.getParameter<double>("Trk_EtaSeparation");
  MaxPadCluInOverlapUporDown= iConfig.getParameter<unsigned int>("MaxPadCluInOverlapUporDown");
  Chi2ProbRefTrk= iConfig.getParameter<double>("Chi2ProbRefTrk");
  UseRestrictedOverlap= iConfig.getParameter<bool>("UseRestrictedOverlap"); 
  // MaxNumbCl1HitInOverlap= iConfig.getParameter<unsigned int>("MaxPadCluInOverlapUporDown");


  LookToRawEvent= iConfig.getParameter<bool>("LookToRawEvent");
  UseUncorrupetdEventMap= iConfig.getParameter<bool>("UseUncorrupetdEventMap");
  xmlfilenameUsed_NotDead= iConfig.getParameter<std::string>("xmlfilenameUsed_NotDead");
  xmlfilenameFull= iConfig.getParameter<std::string>("xmlfilenameFull");

  ReferenceQuarters=iConfig.getParameter<std::vector<unsigned int> >("ReferenceQuarters");  

  TrkEtamin=iConfig.getParameter<double>("TrkEtamin");
  TrkEtaMAX=iConfig.getParameter<double>("TrkEtaMAX");
  AllowedDRTrackDistance=iConfig.getParameter<double>("AllowedDRTrackDistance");   

  ShadowTan=iConfig.getParameter<double>("ShadowTan");
  OnlyShadowAnalysis= iConfig.getParameter<bool>("OnlyShadowAnalysis");//Don't caluclate efficiency, only shadow histo

  UnassHitThr=iConfig.getParameter<double>("UnassHitThr");      
  RefTrkHitMult=iConfig.getParameter<unsigned int>("RefTrkHitMult");
  AnType=iConfig.getParameter<std::string>("AnType");

}


void T2HalfQuarterTrkEfficiency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace HepMC;

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  /* LOADING OF ALL THE RECORDS FROM THE EVENT */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  if((numevent%200)==0)
  std::cout<<"---> numevent: "<<numevent<<std::endl;

    /* :::::::::::::TakeDigi::::::::::::*/
  edm::Handle<T2PadDigiCollection> t2paddigicoll;
  iEvent.getByLabel(t2PadDigiCollectionLabel, t2paddigicoll);
  t2paddigicoll.product();
  edm::Handle<T2StripDigiCollection> t2stripdigicoll;
  iEvent.getByLabel(t2StripDigiCollectionLabel, t2stripdigicoll);
  t2stripdigicoll.product(); 
  DigiContainerIterator<T2DetId, T2PadDigi> itp;
  DigiContainerIterator<T2DetId, T2StripDigi> its;


  edm::Handle<T2VfatInformation> t2vfatinfo;
  
  /* :::::::::::::Take The Clusters::::::::::::*/
  Handle<T2StripClusterCollection> t2strclcoll;
  iEvent.getByLabel(CluLabel,"T2StripClusters",t2strclcoll);
  Handle<T2PadClusterCollection> t2padclcoll;
  iEvent.getByLabel(CluLabel,"T2PadClusters",t2padclcoll);

   
  /*::::::Take  T2  Hits::::::*/
  Handle<T2HitCollection> t2hitcoll;
  iEvent.getByLabel(HitLabel,"T2Hits",t2hitcoll);
  
   //std::cout<<"Analyze-A"<<std::endl;
  Handle<T2PadClusterCollection> Reft2padclcoll;
  iEvent.getByLabel(RefCluLabel,"T2PadClusters",Reft2padclcoll);

  
  if((UseUncorrupetdEventMap))   //Added 1 august 2010 in order to work also with digi
    { 
  
      iEvent.getByLabel(t2VfatInformationLabel, t2vfatinfo);
  
      t2vfatinfo.product();
    
      MaxplaneCorruptInquarter=5;
 
   }

  //std::cout<<"Analyze-B"<<std::endl;


  
 /*:::::: Take Single Arm reconstructed T2 tracks ::::::*/

  //Handle<T2RoadCollection> t2roadcollH0;
  //Handle<T2RoadCollection> t2roadcollH1;
   // iEvent.getByLabel(RoadLabelH1,"T2RoadColl",t2roadcollH1); 
      // iEvent.getByLabel(RoadLabelH0,"T2RoadColl",t2roadcollH0);//T2RoadPadFinderA,"T2Roads",
  
 
  
  //std::cout<<"Analyze-c"<<std::endl;

  //std::cout<<"@@analyze-1"<<std::endl;
  Handle<T1T2TrackCollection> trackCollectionFirsts; 
  Handle<T1T2TrackCollection> trackCollectionLasts;
  
  Handle<T2PadClusterCollection> PadCollectionFirsts; 
  Handle<T2PadClusterCollection> PadCollectionLasts;

    {
     
      iEvent.getByLabel(TrackLabelFirstPlanes,"T2TrackColl",trackCollectionFirsts);
      iEvent.getByLabel(TrackLabelLastPlanes,"T2TrackColl",trackCollectionLasts);

      iEvent.getByLabel(PadCluLabelFirstPlanes,"T2PadClusters",PadCollectionFirsts); 
      iEvent.getByLabel(PadCluLabelLastPlanes,"T2PadClusters",PadCollectionLasts); 
    }

   trackCollectionH0FirstsValid=false;
   trackCollectionH1FirstsValid=false;
   trackCollectionH2FirstsValid=false;
   trackCollectionH3FirstsValid=false;
   trackCollectionH0LastsValid=false;
   trackCollectionH1LastsValid=false;
   trackCollectionH2LastsValid=false;
   trackCollectionH3LastsValid=false;



 //------------------------------------------------------------------------------------
 //                            Track Collection assignment -BEGIN-
 //------------------------------------------------------------------------------------
   
  


 //Track Check
  if(trackCollectionFirsts.isValid())
    {
       if(verbosity)
	 std::cout<<"HFirsts Trk Is assgined "<<std::endl;
      trackCollectionH3FirstsValid=true; trackCollectionH2FirstsValid=true;  trackCollectionH1FirstsValid=true; trackCollectionH0FirstsValid=true; 
      
    }
  else
    {
      if(verbosity)
	std::cout<<"HFirst Trk Is NOT assgined "<<endl;
    }

 



  if(trackCollectionLasts.isValid())
    {
       if(verbosity)
	 std::cout<<"HLasts Trk Is assgined "<<std::endl;
        trackCollectionH3LastsValid=true; trackCollectionH2LastsValid=true;  trackCollectionH1LastsValid=true; trackCollectionH0LastsValid=true;
       
    }
  else
    {
      if(verbosity)
	std::cout<<"HLasts Trk Is NOT assgined "<<endl;
    }




    
  
  Handle<T1T2TrackCollection> trackCollectionRef;
  Handle<T1T2TrackCollection> trackCollectionToMeasure;

  T2SelectionCutUtils T2CutsUtil;
  
  Handle<T2PadClusterCollection> PadCollectionRef; 
  Handle<T2PadClusterCollection> PadCollectionToMeasure;

   
  if((trackCollectionH0FirstsValid)&&(trackCollectionH0LastsValid))
    {
      
      if(Refstring=="A"){
	trackCollectionRef=trackCollectionFirsts;
	trackCollectionToMeasure=trackCollectionLasts; 
	PadCollectionRef=PadCollectionFirsts;
	PadCollectionToMeasure=PadCollectionLasts;  
      }
      else{
	trackCollectionRef=trackCollectionLasts;
	trackCollectionToMeasure=trackCollectionFirsts;
	PadCollectionRef=PadCollectionLasts;
	PadCollectionToMeasure=PadCollectionFirsts;  	
      }
    }
  
  



  
 
 //------------------------------------------------------------------------------------
 //                            Track Collection assignment -END-
 //------------------------------------------------------------------------------------
   
   //std::cout<<"Analyze-f"<<std::endl;
 

 

 //------------------------------------------------------------------------------------
 //                             Check CleanEVENT -BEGIN-
 //------------------------------------------------------------------------------------
 
    double  RsepMin=0.;//was 20
    double PhisepMin=0.;
    unsigned int maxtrkcount=1000;//was 25//was 15
    
unsigned int totalCL1HitInArm0Down=0;
unsigned int totalCL1HitInArm0Up=0;

 std::vector<T1T2Track> CleanRefHalf0Trk;
 std::vector<T1T2Track> CleanRefHalf1Trk; 
 std::vector<T1T2Track> QuarterTestHalf0Trk;
 std::vector<T1T2Track> QuarterTestHalf1Trk;
 
 std::vector<T1T2Track> CleanRefHalf2Trk;
 std::vector<T1T2Track> CleanRefHalf3Trk; 
 std::vector<T1T2Track> QuarterTestHalf2Trk;
 std::vector<T1T2Track> QuarterTestHalf3Trk;

 T1T2TrackCollection::const_iterator TrkCit;
 T2SelectionCutUtils T2CutsUtilForAl;
 //I'm looking if the event is too schifoso (new controls)
  unsigned int count_t2trackVectorALL_H0=0; 
  unsigned int count_t2trackVectorALL_H1=0; 
  unsigned int count_t2trackVectorALL_H2=0; 
  unsigned int count_t2trackVectorALL_H3=0;
  unsigned int trkcounterallEVERYQUARTER=0;
  std::vector<int> dummy_0; dummy_0.push_back(0);	  
  std::vector<int> dummy_1; dummy_1.push_back(1);
  std::vector<int> dummy_2; dummy_2.push_back(2);
  std::vector<int> dummy_3; dummy_3.push_back(3);
  T1T2TrackCollection::const_iterator TrkCit2;
  bool separateTrkInH0=true;  bool separateTrkInH1=true;  bool separateTrkInH2=true;  bool separateTrkInH3=true;

  //trackCollectionRefH0 contain all the ref tracks, 0,1,2,3
  
  //Check how many hit not associated to the trk in the REF. Take only event with Less than 30% unassociated  
  unsigned int sizeref=0;

  unsigned int totHitTrackH0=0;
  unsigned int totHitTrackH1=0;
  unsigned int totHitTrackH2=0;
  unsigned int totHitTrackH3=0;
  std::vector<T2Hit> hits2; 
  for(TrkCit=trackCollectionRef->begin(); TrkCit!=trackCollectionRef->end(); TrkCit++){

    //   sizeref=((*TrkCit).GetT2TrackHits()).size();
    hits2.clear(); 
    T1T2Track test=(*TrkCit);
    hits2=test.GetT2TrackHits();
    sizeref=hits2.size();
    
    
    if((*TrkCit).Eta()>0){
      if(T2CutsUtilForAl.TrkInQuarter((*TrkCit),0)){
	count_t2trackVectorALL_H0++;

      }

      if(T2CutsUtilForAl.TrkInQuarter((*TrkCit),1)){
	count_t2trackVectorALL_H1++;
      }
    }
    else{
      if(T2CutsUtilForAl.TrkInQuarter((*TrkCit),2)){
	count_t2trackVectorALL_H2++;
      }
      if(T2CutsUtilForAl.TrkInQuarter((*TrkCit),3)){
	count_t2trackVectorALL_H3++;
      }
    }

    for(unsigned int kk=0;kk<sizeref;kk++){

      if((RawtoSymb((*TrkCit).GetHitT2(kk).GetHitDetRawId())/10)==0)
	{
	  totHitTrackH0++;
	}
      if((RawtoSymb((*TrkCit).GetHitT2(kk).GetHitDetRawId())/10)==1)
	{
	  totHitTrackH1++;
	}
      if((RawtoSymb((*TrkCit).GetHitT2(kk).GetHitDetRawId())/10)==2)
	{
	  totHitTrackH2++;
	}
      if((RawtoSymb((*TrkCit).GetHitT2(kk).GetHitDetRawId())/10)==3)
	{
	  totHitTrackH3++;
	}        
    }

    
  }

  unsigned int totHitAllH0=0;
  unsigned int totHitAllH1=0;
  unsigned int totHitAllH2=0;
  unsigned int totHitAllH3=0;
  
  int ActivePlH0[10]= {0,0,0,0,0,0,0,0,0,0};
  int ActivePlH1[10]= {0,0,0,0,0,0,0,0,0,0};
  int ActivePlH2[10]= {0,0,0,0,0,0,0,0,0,0}; 
  int ActivePlH3[10]= {0,0,0,0,0,0,0,0,0,0};

  unsigned int intplane=0;

 

  for(T2PadClusterCollection::const_iterator itpad= PadCollectionRef->begin(); itpad != PadCollectionRef->end(); itpad++){//Reft2padclcoll
    vector<T2Cluster> padClv = itpad->second;
    T2DetId *detID =new T2DetId(itpad->first);
    
    uint32_t cmsswdId= detID->calculateRawId(detID->arm(),detID->halfTelescope(),detID->plane(),detID->planeSide());
    unsigned int symb=RawtoSymb(cmsswdId);
    //unsigned int quarter=symbol/10;
    //unsigned int plane=symbol%10;

    if((symb/10)==0){
	totHitAllH0=totHitAllH0+padClv.size();
	intplane=(symb%10);
	ActivePlH0[intplane]=ActivePlH0[intplane]+padClv.size();
      }
      if((symb/10)==1){	
	totHitAllH1=totHitAllH1+padClv.size();
	intplane=(symb%10);
	ActivePlH1[intplane]=ActivePlH1[intplane]+padClv.size();

      }
      if((symb/10)==2){
	totHitAllH2=totHitAllH2+padClv.size();
	intplane=(symb%10);
	ActivePlH2[intplane]=ActivePlH2[intplane]+padClv.size();
      }
      if((symb/10)==3){
	totHitAllH3=totHitAllH3+padClv.size();     
	intplane=(symb%10);
	ActivePlH3[intplane]=ActivePlH3[intplane]+padClv.size();
      }
 }





 /*
  for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++)
    {            
      symb=RawtoSymb(ithit->GetHitDetRawId());

      if((symb/10)==0){
	totHitAllH0++;
	intplane=(symb%10);
	ActivePlH0[intplane]++;
      }
      if((symb/10)==1){	
	totHitAllH1++;
	intplane=(symb%10);
	ActivePlH1[intplane]++;

      }
      if((symb/10)==2){
	totHitAllH2++;
	intplane=(symb%10);
	ActivePlH2[intplane]++;
      }
      if((symb/10)==3){
	totHitAllH3++;     
	intplane=(symb%10);
	ActivePlH3[intplane]++;
      }
    }
 */

  double ActiveH0=0.;
  for(unsigned int i=0;i<10;i++)
    if(ActivePlH0[i]>0){
      ActiveH0=ActiveH0+1.0;
      
    }
  double AvgHitMultH0=0;
  H0RefHitActivePlane->Fill(ActiveH0);


  if(ActiveH0>0)
    AvgHitMultH0=totHitAllH0/ActiveH0;

  double ActiveH1=0.;
  for(unsigned int i=0;i<10;i++)
    if(ActivePlH1[i]>0)
      ActiveH1=ActiveH0+1.0;
  double AvgHitMultH1=0;
  if(ActiveH1>0)
    AvgHitMultH1=totHitAllH1/ActiveH1;
  double ActiveH2=0.;
  for(unsigned int i=0;i<10;i++)
    if(ActivePlH2[i]>0)
      ActiveH2=ActiveH2+1.0;
  double AvgHitMultH2=0;
  if(ActiveH2>0)
    AvgHitMultH2=totHitAllH2/ActiveH2;

  double ActiveH3=0.;
  for(unsigned int i=0;i<10;i++)
    if(ActivePlH3[i]>0)
      ActiveH3=ActiveH3+1.0;
  double AvgHitMultH3=0;
  if(ActiveH3>0)
    AvgHitMultH3=totHitAllH3/ActiveH3;
  
 
 
  //Put here a cut on the clean-event

  AvgRefMultH0->Fill(AvgHitMultH0);
  AvgRefMultH1->Fill(AvgHitMultH1);
  AvgRefMultH2->Fill(AvgHitMultH2);
  AvgRefMultH3->Fill(AvgHitMultH3);
 

  //std::cout<<"  Clean:"<<((totHitAllH0-totHitTrackH0)/((double)totHitTrackH0))<<" "<<separateTrkInH0<<" Allhit:"<<totHitAllH0<<" trkHit:"<<totHitTrackH0<<std::endl;


//def=INCLUSIVE, DOUBLE SINGLE


  if(AnType=="DOUBLE")
    {

      bool doublerespect=false;
      if((count_t2trackVectorALL_H0+count_t2trackVectorALL_H1)>0) 
	if((count_t2trackVectorALL_H2+count_t2trackVectorALL_H3)>0) {
	  doublerespect=true;
	}
      if(doublerespect==false){// A WAY TO DISCARD THE EVENT
	  separateTrkInH0=false;separateTrkInH1=false;separateTrkInH2=false;separateTrkInH3=false; 
      }
    }

  if(AnType=="SINGLE")
    {

      bool singlerespect=false;
      if((count_t2trackVectorALL_H0+count_t2trackVectorALL_H1)>0) 
	if((count_t2trackVectorALL_H2+count_t2trackVectorALL_H3)==0) {
	  singlerespect=true;
	}
      if((count_t2trackVectorALL_H0+count_t2trackVectorALL_H1)==0) 
	if((count_t2trackVectorALL_H2+count_t2trackVectorALL_H3)>0) {
	  singlerespect=true;
	}


      if(singlerespect==false){// A WAY TO DISCARD THE EVENT
	  separateTrkInH0=false;separateTrkInH1=false;separateTrkInH2=false;separateTrkInH3=false; 
      }
    }


    if(count_t2trackVectorALL_H0>maxtrkcount) 
      separateTrkInH0=false;
    if(count_t2trackVectorALL_H1>maxtrkcount)
      separateTrkInH1=false;
    if(count_t2trackVectorALL_H2>maxtrkcount)
      separateTrkInH2=false;
    if(count_t2trackVectorALL_H3>maxtrkcount)
      separateTrkInH3=false;


  for(TrkCit=trackCollectionRef->begin(); TrkCit!=trackCollectionRef->end(); TrkCit++){   

 
   Reference_RawTrkNumHit_NoCut->Fill(TrkCit->GetHitEntries());
  
    
  

  
    
    int quartid=-1;
    if((RawtoSymb((*TrkCit).GetHitT2(0).GetHitDetRawId())/10)==0)
      quartid=0;
    if((RawtoSymb((*TrkCit).GetHitT2(0).GetHitDetRawId())/10)==1)
      quartid=1;
    if((RawtoSymb((*TrkCit).GetHitT2(0).GetHitDetRawId())/10)==2)
      quartid=2;
    if((RawtoSymb((*TrkCit).GetHitT2(0).GetHitDetRawId())/10)==3)
      quartid=3;
    

 

    if((separateTrkInH0)&&(quartid==0))//Continue the check of trk separation in H0
      for(TrkCit2=trackCollectionRef->begin(); TrkCit2!=trackCollectionRef->end(); TrkCit2++){

	if((RawtoSymb((*TrkCit2).GetHitT2(0).GetHitDetRawId())/10)==0)
	  if(TrkCit2!=TrkCit)
	    if(fabs((*TrkCit).GetHitT2(0).GetHitR()-(*TrkCit2).GetHitT2(0).GetHitR())<RsepMin)//Was 50  
	      if(fabs((*TrkCit).GetHitT2(0).GetHitPhi()-(*TrkCit2).GetHitT2(0).GetHitPhi())<PhisepMin)
		{
		  separateTrkInH0=false;		  
		}
      }
    
    if((separateTrkInH1)&&(quartid==1))//Continue the check of trk separation in H1
      for(TrkCit2=trackCollectionRef->begin(); TrkCit2!=trackCollectionRef->end(); TrkCit2++){
	
	if((RawtoSymb((*TrkCit2).GetHitT2(0).GetHitDetRawId())/10)==1)
	  if(TrkCit2!=TrkCit)
	    if(fabs((*TrkCit).GetHitT2(0).GetHitR()-(*TrkCit2).GetHitT2(0).GetHitR())<RsepMin)//Was 50  
	      if(fabs((*TrkCit).GetHitT2(0).GetHitPhi()-(*TrkCit2).GetHitT2(0).GetHitPhi())<PhisepMin)
		{
		  separateTrkInH1=false;
		}
      }

      if((separateTrkInH2)&&(quartid==2))//Continue the check of trk separation in H2
	for(TrkCit2=trackCollectionRef->begin(); TrkCit2!=trackCollectionRef->end(); TrkCit2++){
	
	  if((RawtoSymb((*TrkCit2).GetHitT2(0).GetHitDetRawId())/10)==2)
	    if(TrkCit2!=TrkCit)
	      if(fabs((*TrkCit).GetHitT2(0).GetHitR()-(*TrkCit2).GetHitT2(0).GetHitR())<RsepMin)//Was 50  
		if(fabs((*TrkCit).GetHitT2(0).GetHitPhi()-(*TrkCit2).GetHitT2(0).GetHitPhi())<PhisepMin)
		  {
		    separateTrkInH2=false;
		  }
	}

      if((separateTrkInH3)&&(quartid==3))//Continue the check of trk separation in H3
	for(TrkCit2=trackCollectionRef->begin(); TrkCit2!=trackCollectionRef->end(); TrkCit2++){
	  
	  if((RawtoSymb((*TrkCit2).GetHitT2(0).GetHitDetRawId())/10)==3)
	    if(TrkCit2!=TrkCit)
	      if(fabs((*TrkCit).GetHitT2(0).GetHitR()-(*TrkCit2).GetHitT2(0).GetHitR())<RsepMin)//Was 50  
		if(fabs((*TrkCit).GetHitT2(0).GetHitPhi()-(*TrkCit2).GetHitT2(0).GetHitPhi())<PhisepMin)
		  {
		    separateTrkInH3=false;
		  }
	}
      




    trkcounterallEVERYQUARTER++;

//    double R0xy=sqrt(((*TrkCit).bx_)*((*TrkCit).bx_)+((*TrkCit).by_)*((*TrkCit).by_));      
//    double Z0xy=(*TrkCit).Z_at_Rmin();

    if(quartid==0)
    if((separateTrkInH0==true)){
      unsigned int cl1hitMult=0;
      for(unsigned int i=0;i<TrkCit->GetHitEntries();i++)
	{
	  if(TrkCit->GetHitT2(i).GetHitClass()==1)
	    cl1hitMult++;
	}
	     
      //if((fabs(Z0xy)<8000)&&(R0xy<260))
      if(T2CutsUtil.ChiCutCond((*TrkCit), true,0.01, 0.01))//0.05 in the results sent by email
	if(cl1hitMult>=RefTrkHitMult)				 
	  {
	    //   std::cout<<"Arm0, ov. trks saved !! "<<std::endl;
	    CleanRefHalf0Trk.push_back((*TrkCit));
	  }
    }


    if(quartid==1)
    if((separateTrkInH1==true)){
      unsigned int cl1hitMult=0;
      for(unsigned int i=0;i<TrkCit->GetHitEntries();i++)
	{
	  if(TrkCit->GetHitT2(i).GetHitClass()==1)
	    cl1hitMult++;
	}
	     
      	     
      //if((fabs(Z0xy)<8000)&&(R0xy<260))
      if(T2CutsUtil.ChiCutCond((*TrkCit), true,0.01, 0.01))//0.05 in the results sent by email
	if(cl1hitMult>=RefTrkHitMult)				 
	  {
	    //std::cout<<"Arm0, ov. trks saved !! "<<std::endl;
	    CleanRefHalf1Trk.push_back((*TrkCit));
	  }
    }
    
    if(quartid==2)
    if((separateTrkInH2==true)){
      unsigned int cl1hitMult=0;
      for(unsigned int i=0;i<TrkCit->GetHitEntries();i++)
	{
	  if(TrkCit->GetHitT2(i).GetHitClass()==1)
	    cl1hitMult++;
	}
	     
      	     
      //if((fabs(Z0xy)<8000)&&(R0xy<260))
      if(T2CutsUtil.ChiCutCond((*TrkCit), true,0.01, 0.01))//0.05 in the results sent by email
	if(cl1hitMult>=RefTrkHitMult)				 
	  {
	    //std::cout<<"Arm0, ov. trks saved !! "<<std::endl;
	    CleanRefHalf2Trk.push_back((*TrkCit));
	  }
    }

   if(quartid==3)
    if((separateTrkInH3==true)){
      unsigned int cl1hitMult=0;
      for(unsigned int i=0;i<TrkCit->GetHitEntries();i++)
	{
	  if(TrkCit->GetHitT2(i).GetHitClass()==1)
	    cl1hitMult++;
	}
	     
      	     
      //if((fabs(Z0xy)<8000)&&(R0xy<260))
      if(T2CutsUtil.ChiCutCond((*TrkCit), true,0.01, 0.01))//0.05 in the results sent by email
	if(cl1hitMult>=RefTrkHitMult)				 
	  {
	    //std::cout<<"Arm0, ov. trks saved !! "<<std::endl;
	    CleanRefHalf3Trk.push_back((*TrkCit));
	  }
    }


    
  }//Ref Trk loop End
  
  // std::cout<<"Analyze-f;"<<std::endl;
  NumRefTrk_vs_RefAvgHitMultH0->Fill(count_t2trackVectorALL_H0,AvgHitMultH0);
  NumRefTrk_vs_RefAvgHitMultH1->Fill(count_t2trackVectorALL_H0,AvgHitMultH1);
  NumRefTrk_vs_RefAvgHitMultH2->Fill(count_t2trackVectorALL_H0,AvgHitMultH2);
  NumRefTrk_vs_RefAvgHitMultH3->Fill(count_t2trackVectorALL_H0,AvgHitMultH3);
  
  Num_RefHalfTrkH0->Fill(count_t2trackVectorALL_H0);
  Num_RefHalfTrkH1->Fill(count_t2trackVectorALL_H1);
  Num_RefHalfTrkH2->Fill(count_t2trackVectorALL_H2);
  Num_RefHalfTrkH3->Fill(count_t2trackVectorALL_H3);
  
   
  
 //------------------------------------------------------------------------------------
 //                             Check CleanEVENT -END-
 //------------------------------------------------------------------------------------
 

 //-----------------------------------------------------------------------------------------------
 //                            Create a vector of good reference track and testing track  -BEGIN-
 //------------------------------------------------------------------------------------------------

   
 
   NumCl1HitInOvRegionArm0_Up->Fill(totalCL1HitInArm0Up);
   NumCl1HitInOvRegionArm0_Down->Fill(totalCL1HitInArm0Down);
   
  


   
       //-----------------------------------------------------------------------------------------------
       //                            End of Create a vector of good reference track and testing track  -END-
       //------------------------------------------------------------------------------------------------
 
       
       //std::cout<<"Reference Track Loaded"<<std::endl;

   //Filling the Tests samples
 
       
	 for(T1T2TrackCollection::const_iterator TrkCit=trackCollectionToMeasure->begin(); TrkCit!=trackCollectionToMeasure->end(); TrkCit++)
	   {   

	      if((*TrkCit).Eta()>0){
		if(T2CutsUtilForAl.TrkInQuarter((*TrkCit),0)){
		  QuarterTestHalf0Trk.push_back((*TrkCit));
		}

		if(T2CutsUtilForAl.TrkInQuarter((*TrkCit),1)){
		  QuarterTestHalf1Trk.push_back((*TrkCit));
		}
	      }
	      else{
		if(T2CutsUtilForAl.TrkInQuarter((*TrkCit),2)){
		  QuarterTestHalf2Trk.push_back((*TrkCit));
		}
		if(T2CutsUtilForAl.TrkInQuarter((*TrkCit),3)){
		  QuarterTestHalf3Trk.push_back((*TrkCit));
		}
	      }
	          
    
	   }
   
   

  
   //------------------------------------------------------------------------------------
   //                             EFFICIENCY CALCULATION -BEGIN-
   //------------------------------------------------------------------------------------
  
  //Check conditions for include this event in the sample used for efficiency calculation
 
 //Test the H1 efficiency: condition only in H0 Occupancy
   // std::cout<<"Arm 0 clean reference track size: "<<CleanRefHalf0Trk.size()<<std::endl;



 if(CleanRefHalf0Trk.size()>0)
   {               
     double investigatethisevent=0.;
     // (*investigatethisevent)=0.;    
     // std::cout<<"HERE-1"<<std::endl;
     CalculateEfficiency(&CleanRefHalf0Trk ,&QuarterTestHalf0Trk,0,0,"--",t2hitcoll,count_t2trackVectorALL_H0,AvgHitMultH0,investigatethisevent); 
     // std::cout<<"HERE-2"<<std::endl;std::cout<<"HERE-2"<<std::endl;
     
     if(investigatethisevent==1){        
       //  std::cout<<"HERE-3"<<std::endl;std::cout<<"HERE-3"<<std::endl;
       for(T2PadClusterCollection::const_iterator itpad= PadCollectionToMeasure->begin(); itpad !=PadCollectionToMeasure->end(); itpad++){
	 vector<T2Cluster> padClv = itpad->second;
	 T2DetId *detID =new T2DetId(itpad->first);
	 
	 uint32_t cmsswdId= detID->calculateRawId(detID->arm(),detID->halfTelescope(),detID->plane(),detID->planeSide());
	 unsigned int symb=RawtoSymb(cmsswdId);
	 //unsigned int quarter=symbol/10;
	 //unsigned int plane=symbol%10;
	 //std::cout<<"HERE-4"<<std::endl;
	 if((symb/10)==0){ 
	   unsigned int intplane=(symb%10);
	   std::cout<<"Plane: "<<intplane<<std::endl;
	   for(unsigned ff=0;ff<padClv.size();ff++){

	     double rr=padClv.at(ff).GetClusterR(); double phii=padClv.at(ff).GetClusterPhi();
	     double xx=rr*cos(phii*3.14159/180.); double yy=rr*sin(phii*3.14159/180.);
	     	   
	     
	     std::cout<<" | X-Y: "<<xx<<" "<<yy<<"      "<<std::endl;
	   }
	 }
	 
       }
     }

   }
	
 if(CleanRefHalf1Trk.size()>0)
   {     
     double investigatethisevent=0.;    
     CalculateEfficiency(&CleanRefHalf1Trk ,&QuarterTestHalf1Trk,1,1,"--",t2hitcoll,count_t2trackVectorALL_H1,AvgHitMultH1,investigatethisevent);     
     }

 if(CleanRefHalf2Trk.size()>0)
   {               
     double investigatethisevent=0.;
     CalculateEfficiency(&CleanRefHalf2Trk ,&QuarterTestHalf2Trk,2,2,"--",t2hitcoll,count_t2trackVectorALL_H2,AvgHitMultH2,investigatethisevent);     
   }

 if(CleanRefHalf3Trk.size()>0)
   {               
     double investigatethisevent=0.;
     CalculateEfficiency(&CleanRefHalf3Trk ,&QuarterTestHalf3Trk,3,3,"--",t2hitcoll,count_t2trackVectorALL_H3,AvgHitMultH3,investigatethisevent);     
   }

  
 
  
 //------------------------------------------------------------------------------------
 //                             EFFICIENCY CALCULATION -END-
 //------------------------------------------------------------------------------------
  


  
 if(verbosity) std::cout<<"numevent"<<numevent<<std::endl;

 
numevent++;


}












unsigned int T2HalfQuarterTrkEfficiency::Howmanyplanes(std::vector<unsigned int> *total_cl1Hit_CloseInArm0Up,unsigned int testedQuarter)
{
  unsigned int start=testedQuarter*10;
  unsigned int totcount=0;
  for(unsigned int i=start;i<start+10;i++)
    {
      if(total_cl1Hit_CloseInArm0Up->at(i)>0)
	totcount++;
    }
  return totcount;			
}






std::vector<double> T2HalfQuarterTrkEfficiency::MyLinearfitY(std::vector<T2Hit> hitvec,unsigned int UseJointProb)
{  

double Sy=0.;
double Syz=0.;
double Szz_y=0.;
double Sz_y=0.; 
double S0_y=0.; 



unsigned int  sizeHitv=hitvec.size();
std::vector<double> y;
std::vector<double> z;

std::vector<double> ey;
 std::vector<double> ez;

 double a_yz, b_yz;

  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      //   std::cout<<hitvec[jj].GetHitR()<<" "<<hitvec[jj].GetHitPhi()<<" "<<hitvec[jj].GetHitZ()<<std::endl;
      y.push_back(hitvec[jj].GetHitY());
      

      z.push_back(hitvec[jj].GetHitZ());    
 
      double phirad=hitvec[jj].GetHitPhi()*3.14159/180.0;
      double sigmay=sin(phirad)*sin(phirad)*0.12*0.12+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
      
      //     if(UseJointProb==1)
      //sigmay=(0.12*hitvec[jj].GetHitR())*(0.12/**hitvec[jj].GetHitR()*/);
	
      ey.push_back(sqrt(sigmay));
      
      // er.push_back(0.1);
      
      //ez.push_back(hitvec[jj].GetHitDZ());
            
      Syz += y[jj]*z[jj]/ey[jj]/ey[jj];
      Szz_y += z[jj]*z[jj]/ey[jj]/ey[jj];
      Sz_y += z[jj]/ey[jj]/ey[jj];
      Sy += y[jj]/ey[jj]/ey[jj];
      S0_y += 1.0/ey[jj]/ey[jj];

  }


  a_yz = (Syz*S0_y - Sz_y*Sy) / (Szz_y*S0_y - Sz_y*Sz_y);   // angulay coefficient
  b_yz = (Sy*Szz_y - Sz_y*Syz) / (Szz_y*S0_y - Sz_y*Sz_y);  // intercept   R=(a_rz)Z + b_rz
 double e_a_yz = sqrt( S0_y / (S0_y*Szz_y - Sz_y*Sz_y) );
 double e_b_yz = sqrt( Szz_y / (S0_y*Szz_y - Sz_y*Sz_y) );

  

 std::vector<double> vect;
 for (unsigned int i=0;i<5;i++)
   vect.push_back(0.0);

 vect[0]=a_yz;
 vect[1]=b_yz;
 vect[2]=e_a_yz;
 vect[3]=e_b_yz;

 double correl=(-1.0)*Sz_y*(1.0/(S0_y*Szz_y-(Sz_y*Sz_y)));
 vect[4]=correl;


 return vect;


}
std::vector<double> T2HalfQuarterTrkEfficiency::MyLinearfitX(std::vector<T2Hit> hitvec,unsigned int UseJointProb)
{  

double Sx=0.;
double Sxz=0.;
double Szz_x=0.;
double Sz_x=0.; 
double S0_x=0.; 



unsigned int  sizeHitv=hitvec.size();
std::vector<double> x;
std::vector<double> z;

std::vector<double> ex;
 std::vector<double> ez;

 double a_xz, b_xz;

  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      //   std::cout<<hitvec[jj].GetHitR()<<" "<<hitvec[jj].GetHitPhi()<<" "<<hitvec[jj].GetHitZ()<<std::endl;
      x.push_back(hitvec[jj].GetHitX());
      

      z.push_back(hitvec[jj].GetHitZ());    
 

      //ex.push_back(0.09*fabs(hitvec[jj].GetHitX()));
      double phirad=hitvec[jj].GetHitPhi()*3.14159/180.0;
      double sigmax=cos(phirad)*cos(phirad)*0.12*0.12+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
      
      //if(UseJointProb==1)
      //sigmax=(0.12*hitvec[jj].GetHitR())*(0.12/**hitvec[jj].GetHitR()*/);
      
      ex.push_back(sqrt(sigmax));
 
            
      Sxz += x[jj]*z[jj]/ex[jj]/ex[jj];
      Szz_x += z[jj]*z[jj]/ex[jj]/ex[jj];
      Sz_x += z[jj]/ex[jj]/ex[jj];
      Sx += x[jj]/ex[jj]/ex[jj];
      S0_x += 1.0/ex[jj]/ex[jj];

  }


  a_xz = (Sxz*S0_x - Sz_x*Sx) / (Szz_x*S0_x - Sz_x*Sz_x);   // angular coefficient
  b_xz = (Sx*Szz_x - Sz_x*Sxz) / (Szz_x*S0_x - Sz_x*Sz_x);  // intercept   R=(a_rz)Z + b_rz

 double e_a_xz = sqrt( S0_x / (S0_x*Szz_x - Sz_x*Sz_x) );
 double e_b_xz = sqrt( Szz_x / (S0_x*Szz_x - Sz_x*Sz_x) );

 std::vector<double> vect;
 for (unsigned int i=0;i<5;i++)
   vect.push_back(0.0);

 vect[0]=a_xz;
 vect[1]=b_xz;
 vect[2]=e_a_xz;
 vect[3]=e_b_xz;

 double correl=(-1.0)*Sz_x*(1.0/(S0_x*Szz_x-(Sz_x*Sz_x)));
 vect[4]=correl;

 return vect;


}


std::vector<double> T2HalfQuarterTrkEfficiency::TrkPredictedPoint(T1T2Track &thetrk, double planeZ)
{

  std::vector<T2Hit> hitvec;
  unsigned int sizeHitv=thetrk.GetHitEntries();
  for(unsigned int i=0;i<sizeHitv;i++)
    {
      hitvec.push_back(thetrk.GetHitT2(i));
    }

  std::vector<double> vparx=MyLinearfitX(hitvec,false); 
  std::vector<double> vpary=MyLinearfitY(hitvec,false);
  std::vector<double> PointPrediction;
  
 

  double e_a_xz =vparx.at(2) ;
 double e_b_xz = vparx.at(3) ;

 double e_a_yz = vpary.at(2);
 double e_b_yz = vpary.at(3);

double correlx=vparx.at(4);
double correly=vpary.at(4);

 


 double xpred=vparx.at(0)*(planeZ)+vparx.at(1);
 double ypred=vpary.at(0)*(planeZ)+vpary.at(1);;
 double Expred=e_b_xz*e_b_xz + planeZ*planeZ*e_a_xz*e_a_xz + 2*planeZ*correlx;
 double Eypred=e_b_yz*e_b_yz + planeZ*planeZ*e_a_yz*e_a_yz + 2*planeZ*correly;

 

 if((Expred*Eypred)<0)
   std::cout<<"T2RoadPadFinder::TubePredictedPoint error: error on point prediction <0 "<<std::endl;
 
 Expred=sqrt(Expred);
 Eypred=sqrt(Eypred);
 

 /*
if(((Expred)>40)||((Eypred)>40)){
   std::cout<<"T2RoadPadFinder::TubePredictedPoint Warning: Expred-Eypred:"<<Expred<<"-"<<Eypred<<" NumGoodHit:"<<countergoodHit<<" EAx-EAy: "<<e_a_xz<<" "<<e_a_yz<<" Eb_x Eb_y: "<<e_b_xz<<" "<<e_b_yz<<" Distance:"<<(planeZ-zgrav)<<std::endl;
   for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      if(atube.ClustV[jj].GetNoOfEntries()<MinCluSize_considered_asBlobs){
	std::cout<<"X:"<<atube.VectX[jj]<<" Ex:"<<atube.VectEX[jj]<<" Y:"<<atube.VectY[jj]<<" Ey:"<<atube.VectEY[jj]<<std::endl;
      }
    }
 }
 */


   // x, y, ex, ey
 PointPrediction.push_back(xpred); 
 PointPrediction.push_back(ypred);
 PointPrediction.push_back(Expred); 
 PointPrediction.push_back(Eypred);
    
 //std::cout<<"Prediction X-Y | Ex-Ey "<<xpred<<"-"<<ypred<<" | "<<Expred<<"-"<<Eypred<<"  CorrX:"<<correlx<<"CorrY:"<<correly<<"  E_ax="<<e_a_xz<<"  E_bx="<<e_b_xz<<" Zjump:"<<(planeZ-zgrav)<<std::endl;
 return PointPrediction;
}





void T2HalfQuarterTrkEfficiency::CalculateEfficiency(std::vector<T1T2Track>* CleanRefHalfXTrk ,std::vector<T1T2Track>* TestArmXTrk, unsigned int refQuarter, unsigned int testedQuarter,std::string Up_or_DownY,edm::Handle<T2HitCollection> t2hitcoll_,int numreftrks_raw, double avgMult,double &investigatethisevent)
{
  /*
  std::cout<<"CalculateEfficiency for quarter H"<<testedQuarter<<" using reference H"<<refQuarter<<" "<<std::endl;
  std::cout<<"Reference: "<<CleanRefHalfXTrk->size()<<" Testing "<<TestArmXTrk->size()<<std::endl;
  std::cout<<"-----------------------------------------------------------------------------"<<std::endl;
  std::cout<<"numevent: "<<numevent<<std::endl;
  std::cout<<"-----------------------------------------------------------------------------"<<std::endl;
*/
  T2SelectionCutUtils T2CutUtils;

  unsigned int numreftrks= 0;  
  unsigned int trackMatch= 0;
  
  unsigned int numreftrksEtaCut= 0;
  unsigned int numreftrksZImpactCut= 0;
  unsigned int trackMatchEtaCut= 0;
  unsigned int trackMatchZImpactCut= 0;
  
   
  double trkR1=0.;
    for(unsigned int i=0;i<CleanRefHalfXTrk->size();i++)
      {

	//	unsigned int numreftrksZImpactCut_forEta= 0;
	//unsigned int numreftrksZImpactCut_forEta_56Division= 0;
	//	unsigned int internalcounterZImpactCut_forEta= 0;

	double trkphi1=T2CutUtils.PhiFitAverage(CleanRefHalfXTrk->at(i).GetT2TrackHits());
	double trketa1=T2CutUtils.EtaFromAveragesVtxCorr(CleanRefHalfXTrk->at(i),0.,0.,0.);
	trkR1=0.;
	unsigned int sizeref=CleanRefHalfXTrk->at(i).GetT2TrackHits().size();
	double trkX1=0.;
	double trkY1=0.;


	TrkRinRefQuarter->Fill(CleanRefHalfXTrk->at(i).GetHitT2(0).GetHitR());
	TrkPhiinRefQuarter->Fill(trkphi1);
	ReferenceTrackXY->Fill(CleanRefHalfXTrk->at(i).GetHitT2(0).GetHitX(),CleanRefHalfXTrk->at(i).GetHitT2(0).GetHitY());
	TrkXinRefQuarter->Fill(CleanRefHalfXTrk->at(i).GetHitT2(0).GetHitX());
	TrkYinRefQuarter->Fill(CleanRefHalfXTrk->at(i).GetHitT2(0).GetHitY());

	for(unsigned int jj=0;jj<sizeref;jj++)
	  {
	    trkR1+=CleanRefHalfXTrk->at(i).GetHitT2(jj).GetHitR();
	    trkX1+=CleanRefHalfXTrk->at(i).GetHitT2(jj).GetHitX();
	    trkY1+=CleanRefHalfXTrk->at(i).GetHitT2(jj).GetHitY();  
	  }
	
	trkR1=trkR1/((double)sizeref);  
	trkY1=trkY1/((double)sizeref);  
	trkX1=trkX1/((double)sizeref);  

	/*
	  std::cout<<"New Reference Track <eta>-phi-R: "<<trketa1<<" "<<trkphi1<<" "<<trkR1<<", looking for a match"<<std::endl;
	*/
	double trkphi2=0.;
	double trkR2=0.;

	double TrkDX=0.;
	double TrkDY=0.;
	unsigned int internalcounter=0;
	unsigned int internalcounterEtaCut=0;
	unsigned int internalcounterZImpactCut=0;
	bool cutRefEta=false;

	numreftrks++;
	if(fabs(CleanRefHalfXTrk->at(i).Eta())>4.5){
	  numreftrksEtaCut++;
	  cutRefEta=true;
	}

	


	bool cutRefZImpact=false; double Zimp=0.;
	cutRefZImpact=RecoTrk_UnfoldingCut(CleanRefHalfXTrk->at(i),Zimp);

	if(refQuarter==0)	  
	  ZIMPdistrRefH0->Fill(Zimp);
	if(refQuarter==1)
	  ZIMPdistrRefH1->Fill(Zimp);
	if(refQuarter==2)	  
	  ZIMPdistrRefH2->Fill(Zimp);
	if(refQuarter==3)
	  ZIMPdistrRefH3->Fill(Zimp);
	
	if(cutRefZImpact){
	  numreftrksZImpactCut++;		  

	  //  numreftrksZImpactCut_forEta++;
	  if(refQuarter==0)	  
	    ReferenceZImpVsEta2->Fill(Zimp,trketa1);
	}

	//bool cutRefZImpactDivision=false;
	//	cutRefZImpactDivision=RecoTrk_UnfoldingCut_56Division(CleanRefHalfXTrk->at(i));
	//	if(cutRefZImpactDivision==true){
	// numreftrksZImpactCut_forEta_56Division++;

	
	  //std::cout<<"so it pass"<<std::endl;
	//}

	for(unsigned int j=0;j<TestArmXTrk->size();j++)//TestArmXTrk
	  {
	    //double C0=((the_t2ev->TrkEntryX.at(j)*the_t2ev->TrkEntryX.at(j)+the_t2ev->TrkEntryY.at(j)*the_t2ev->TrkEntryY.at(j))/(the_t2ev->TrkAx.at(j)*the_t2ev->TrkEntryX.at(j)+the_t2ev->TrkAy.at(j)*the_t2ev->TrkEntryY.at(j)));
	    //double Z0impact=the_t2ev->TrkEntryZ.at(j)-C0;
	    
	    //double C0=(atrk.GetHitT2(0).GetHitX()*atrk.GetHitT2(0).GetHitX()+atrk.GetHitT2(0).GetHitY()*atrk.GetHitT2(0).GetHitY())/(atrk.GetTx()*atrk.GetHitT2(0).GetHitX()+atrk.GetTy()*atrk.GetHitT2(0).GetHitY());
	    //double Z0impact=atrk.GetHitT2(0).GetHitZ()-C0;

	    //
	    unsigned int sizetest=TestArmXTrk->at(j).GetT2TrackHits().size();
	    bool TrkMatch=false;
	    unsigned int numhitmatch=0;

	    for(unsigned kk=0;kk<TestArmXTrk->at(j).GetT2TrackHits().size();kk++){

	      double z=TestArmXTrk->at(j).GetHitT2(kk).GetHitZ();	    	    
	      double xp=z*CleanRefHalfXTrk->at(i).GetTx()+CleanRefHalfXTrk->at(i).X0();
	      double yp=z*CleanRefHalfXTrk->at(i).GetTy()+CleanRefHalfXTrk->at(i).Y0();
	      trkR1=sqrt(xp*xp+yp*yp);
	  

	      trkphi2=T2CutUtils.PhiFitAverage(TestArmXTrk->at(j).GetT2TrackHits());
	      trkR2=TestArmXTrk->at(j).GetHitT2(kk).GetHitR();//TestArmXTrk->at(j).GetHitT2(0).GetHitR();


	      //Index to change Accord to the studied Quarter
	      //  trkR2=CleanRefHalfXTrk->at(i).GetHitT2(0).GetHitR();
	      TrkDY=yp-TestArmXTrk->at(j).GetHitT2(kk).GetHitY();//TestArmXTrk->at(j).GetHitT2(0).GetHitY();//CleanRefHalfXTrk->at(i).GetHitT2(0).GetHitY();
	      TrkDX=xp-TestArmXTrk->at(j).GetHitT2(kk).GetHitX();//TestArmXTrk->at(j).GetHitT2(0).GetHitX();
	      
	      DPhiTrk_HRefHTest->Fill((trkphi1-trkphi2));
	      DRTrk_HRefHTest->Fill((trkR2-trkR1)/TestArmXTrk->at(j).GetHitT2(kk).GetHitDR());
	      
	      
	      //if((fabs(trkR2-trkR1)<3*TestArmXTrk->at(j).GetHitT2(kk).GetHitDR())&&(fabs(TrkDY)<3*TestArmXTrk->at(j).GetHitT2(kk).GetHitDY())&&(fabs(TrkDX)<3*TestArmXTrk->at(j).GetHitT2(kk).GetHitDX())){
	      //	TrkMatch=true;
	      //}
		
		std::vector<double> predPointAndERR=TrkPredictedPoint(CleanRefHalfXTrk->at(i),z);
		//PointPrediction.push_back(xpred); 
		//PointPrediction.push_back(ypred);
		//PointPrediction.push_back(Expred); 
		//PointPrediction.push_back(Eypred);
		

 
		TrkDY=predPointAndERR.at(1)-TestArmXTrk->at(j).GetHitT2(kk).GetHitY();
		TrkDX=predPointAndERR.at(0)-TestArmXTrk->at(j).GetHitT2(kk).GetHitX();
		double erry=predPointAndERR.at(3);
		double errx=predPointAndERR.at(2);
		
		errx=errx*errx+TestArmXTrk->at(j).GetHitT2(kk).GetHitDX()*TestArmXTrk->at(j).GetHitT2(kk).GetHitDX();
		erry=erry*erry+TestArmXTrk->at(j).GetHitT2(kk).GetHitDY()*TestArmXTrk->at(j).GetHitT2(kk).GetHitDY();
		 
		errx=sqrt(errx);
		erry=sqrt(erry);
		
		if((fabs(TrkDY)<3*erry)&&(fabs(TrkDX)<3*errx))
		  numhitmatch++;
	      
	    }//HitInTrack for end

	    //WARNINGGGGGGGGGGGGGGG
	    if(numhitmatch>2)
	      TrkMatch=true;
	 
	    if(TrkMatch)
	      {
		
		DPhiDRTrk_HRefHTest->Fill(fabs(trkphi1-trkphi2),fabs(trkR1-trkR2));
		DPhiTrk_HRefHTest2CloseToRef->Fill((trkphi1-trkphi2));
		TrkRinTestedQuarterCloseToRef->Fill(trkR2);
		TrkPhiinTestedQuarterCloseToRef->Fill(TestArmXTrk->at(j).Phi());
		TrkXinTestedQuarterCloseToRef->Fill(TestArmXTrk->at(j).GetHitT2(0).GetHitX());
		TrkYinTestedQuarterCloseToRef->Fill(TestArmXTrk->at(j).GetHitT2(0).GetHitY());

		QuarterQuarterDYExit->Fill(CleanRefHalfXTrk->at(i).GetHitT2(sizeref-1).GetHitY()-TestArmXTrk->at(j).GetHitT2(sizetest-1).GetHitY());
		QuarterQuarterDYEntry->Fill(CleanRefHalfXTrk->at(i).GetHitT2(0).GetHitY()-TestArmXTrk->at(j).GetHitT2(0).GetHitY());
		QuarterQuarterDXExit->Fill(CleanRefHalfXTrk->at(i).GetHitT2(sizeref-1).GetHitX()-TestArmXTrk->at(j).GetHitT2(sizetest-1).GetHitX());
		QuarterQuarterDXEntry->Fill(CleanRefHalfXTrk->at(i).GetHitT2(0).GetHitX()-TestArmXTrk->at(j).GetHitT2(0).GetHitX());
		
		for(unsigned int jj=0;jj<sizetest;jj++)
		  {
		    DPhiTrkHit_HRefHTestCloseVsZTest->Fill(fabs(TestArmXTrk->at(j).GetHitT2(jj).GetHitZ()),trkphi1-TestArmXTrk->at(j).GetHitT2(jj).GetHitPhi());
		    DRTrkHit_HRefHTestCloseVsZTest->Fill(fabs(TestArmXTrk->at(j).GetHitT2(jj).GetHitZ()),trkR1-TestArmXTrk->at(j).GetHitT2(jj).GetHitR());
		    DXTrkHit_HRefHTestCloseVsZTest->Fill(fabs(TestArmXTrk->at(j).GetHitT2(jj).GetHitZ()),trkX1-TestArmXTrk->at(j).GetHitT2(jj).GetHitX());
		    DYTrkHit_HRefHTestCloseVsZTest->Fill(fabs(TestArmXTrk->at(j).GetHitT2(jj).GetHitZ()),trkY1-TestArmXTrk->at(j).GetHitT2(jj).GetHitY());

		    DYTrkHit_HRefHTestClose->Fill(trkY1-TestArmXTrk->at(j).GetHitT2(jj).GetHitY());  
		    DXTrkHit_HRefHTestClose->Fill(trkX1-TestArmXTrk->at(j).GetHitT2(jj).GetHitX()); 
		    DRTrkHit_HRefHTestClose->Fill(trkR1-TestArmXTrk->at(j).GetHitT2(jj).GetHitR());     
		    DPhiTrkHit_HRefHTestClose->Fill(trkphi1-TestArmXTrk->at(j).GetHitT2(jj).GetHitPhi());     
		    
		  }
		
		//double trketa2= T2CutUtils.EtaFromAverages(TestArmXTrk->at(j));
		 //std::cout<<"OKKKKKKKKK Y!!!!  phi1="<<trkphi1<<" phi2="<<trkphi2<<" eta1="<<trketa1<<" eta2="<<trketa2<<" ra:"<<trkR1-TestArmXTrk->at(j).GetHitT2(0).GetHitR()<<" rb:"<<trkR1-TestArmXTrk->at(j).GetHitT2(1).GetHitR()<<std::endl;			
		 
		internalcounter++;
		
		if(cutRefEta)
		  internalcounterEtaCut++;

		if(cutRefZImpact){
		  internalcounterZImpactCut++;
		  // internalcounterZImpactCut_forEta++;
		}
	      } //If Trk_Match
	    
   
	  }//End-For TestTrack
	
	if(internalcounter>0)
	  trackMatch++;
	if(internalcounterEtaCut>0)
	  trackMatchEtaCut++;
	if(internalcounterZImpactCut>0)
	  trackMatchZImpactCut++;
	
	
	unsigned int intmultclass=(unsigned int)floor(avgMult/5.);
	if(intmultclass>9)//High multiplicity event goes in the same bin;
	  intmultclass=9;
	
	if(numreftrksZImpactCut>0){

	  //trackMatchZImpactCut++; H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_ZimpCut
	  if((refQuarter==0)&&(internalcounterZImpactCut>0))
	    H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(trketa1,1);
	  
	  if((refQuarter==0)&&(internalcounterZImpactCut==0))
	    H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(trketa1,0);


	  if((refQuarter==1)&&(internalcounterZImpactCut>0))
	    H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(trketa1,1);
	  
	  if((refQuarter==1)&&(internalcounterZImpactCut==0))
	    H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(trketa1,0);

	  if((refQuarter==2)&&(internalcounterZImpactCut>0))
	    H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(trketa1,1);
	  
	  if((refQuarter==2)&&(internalcounterZImpactCut==0))
	    H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(trketa1,0);


	  if((refQuarter==3)&&(internalcounterZImpactCut>0))
	    H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(trketa1,1);
	  
	  if((refQuarter==3)&&(internalcounterZImpactCut==0))
	    H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(trketa1,0);

	}


	/*
	if(numreftrksZImpactCut_forEta_56Division>0){
	  //std::cout<<"so so it pass. Trk eta="<<trketa1<<std::endl;
	  //trackMatchZImpactCut_56Division++;
	  if((refQuarter==0)&&(internalcounter>0))
	    H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_56Division[intmultclass]->Fill(trketa1,1);
	  if((refQuarter==0)&&(internalcounter==0))
	    H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_56Division[intmultclass]->Fill(trketa1,0);
	}*/
	
      }//end reference Trk-for
    
    
   



    if(numreftrks>0){

      HXTrkEffi->Fill(numreftrks,(trackMatch/numreftrks));
      HXTrkEffi_AxisRawTrk->Fill(numreftrks_raw,(trackMatch/numreftrks));

      if(refQuarter==0){
	H0TrkEffi->Fill(numreftrks,(trackMatch/numreftrks));
	H0TrkEffi_AxisRawTrk->Fill(numreftrks_raw,(trackMatch/numreftrks));	
	H0TrkEffi_AxisAvgMult->Fill(avgMult,(trackMatch/numreftrks));
      }

      if(refQuarter==1){
	H1TrkEffi->Fill(numreftrks,(trackMatch/numreftrks));
	H1TrkEffi_AxisRawTrk->Fill(numreftrks_raw,(trackMatch/numreftrks));
	H1TrkEffi_AxisAvgMult->Fill(avgMult,(trackMatch/numreftrks));
      }
    
      if(refQuarter==2){
	H2TrkEffi->Fill(numreftrks,(trackMatch/numreftrks));
	H2TrkEffi_AxisRawTrk->Fill(numreftrks_raw,(trackMatch/numreftrks));
	H2TrkEffi_AxisAvgMult->Fill(avgMult,(trackMatch/numreftrks));
      }

      if(refQuarter==3){
	H3TrkEffi->Fill(numreftrks,(trackMatch/numreftrks));
	H3TrkEffi_AxisRawTrk->Fill(numreftrks_raw,(trackMatch/numreftrks));
	H3TrkEffi_AxisAvgMult->Fill(avgMult,(trackMatch/numreftrks));
      }


     
      if(numreftrksZImpactCut>0){
	if(refQuarter==0){
	  H0TrkEffi_AxisAvgMult_ZImpactCut->Fill(avgMult,(trackMatchZImpactCut/numreftrksZImpactCut));

	}
	if(refQuarter==1)
	  H1TrkEffi_AxisAvgMult_ZImpactCut->Fill(avgMult,(trackMatchZImpactCut/numreftrksZImpactCut));
	if(refQuarter==2)
	  H2TrkEffi_AxisAvgMult_ZImpactCut->Fill(avgMult,(trackMatchZImpactCut/numreftrksZImpactCut));
	if(refQuarter==3)
	  H3TrkEffi_AxisAvgMult_ZImpactCut->Fill(avgMult,(trackMatchZImpactCut/numreftrksZImpactCut));
      }

      if(numreftrksEtaCut>0){

	if(refQuarter==0){
	  H0TrkEffi_AxisRawTrk_EtaCut->Fill(numreftrks_raw,(trackMatchEtaCut/numreftrksEtaCut));
	  H0TrkEffi_AxisAvgMult_EtaCut->Fill(avgMult,(trackMatchEtaCut/numreftrksEtaCut));	 
	}

	if(refQuarter==1){
	  H1TrkEffi_AxisRawTrk_EtaCut->Fill(numreftrks_raw,(trackMatchEtaCut/numreftrksEtaCut));
	  H1TrkEffi_AxisAvgMult_EtaCut->Fill(avgMult,(trackMatchEtaCut/numreftrksEtaCut));

	}
    
	if(refQuarter==2){
	  H2TrkEffi_AxisRawTrk_EtaCut->Fill(numreftrks_raw,(trackMatchEtaCut/numreftrksEtaCut));
	  H2TrkEffi_AxisAvgMult_EtaCut->Fill(avgMult,(trackMatchEtaCut/numreftrksEtaCut));
	}

	if(refQuarter==3){
	  H3TrkEffi_AxisRawTrk_EtaCut->Fill(numreftrks_raw,(trackMatchEtaCut/numreftrksEtaCut));
	  H3TrkEffi_AxisAvgMult_EtaCut->Fill(avgMult,(trackMatchEtaCut/numreftrksEtaCut));	
	}


      }


      
    }
	
  

    //   std::cout<<"CalculateEfficiency-END"<<std::endl;
}


bool T2HalfQuarterTrkEfficiency::RecoTrk_UnfoldingCut_56Division(T1T2Track &atrk){

bool cutpassed=false;

  unsigned int counter=0; 
  T2SelectionCutUtils T2CutsUtil;
  //std::cout<<"Caso0bb RecoTrk_UnfoldingCut"<<atrk.GetHitEntries()<<std::endl;
  for(unsigned int u=0;u<atrk.GetHitEntries();u++)
    if(atrk.GetHitT2(u).GetHitClass()==1)
      counter++;
  //std::cout<<"Caso0cc RecoTrk_UnfoldingCut"<<atrk.GetHitEntries()<<std::endl;
  double C0=(atrk.GetHitT2(0).GetHitX()*atrk.GetHitT2(0).GetHitX()+atrk.GetHitT2(0).GetHitY()*atrk.GetHitT2(0).GetHitY())/(atrk.GetTx()*atrk.GetHitT2(0).GetHitX()+atrk.GetTy()*atrk.GetHitT2(0).GetHitY());
  double Z0impact=atrk.GetHitT2(0).GetHitZ()-C0;

  double eta2=T2CutsUtil.EtaFromAveragesVtxCorr(atrk,0.,0.,0.);
  if((fabs(atrk.Z_at_Rmin())<13500)||((atrk.Z_at_Rmin())*atrk.GetHitT2(0).GetHitZ()<0)){
    if((fabs(Z0impact)<5000)&&(eta2>5.6))
      cutpassed=true;

    if((fabs(Z0impact)<10000)&&(eta2<5.6))
      cutpassed=true;
  }

  //if(eta2<5.6)
  //std::cout<<"56Division cut: eta2="<<eta2<< " Returned: "<<cutpassed<<std::endl;
    //else
    

  return cutpassed;
}


bool T2HalfQuarterTrkEfficiency::RecoTrk_UnfoldingCut(T1T2Track &atrk, double &Zimp){

  bool cutpassed=false;
  unsigned int counter=0;
  
  for(unsigned int u=0;u<atrk.GetHitEntries();u++)
    if(atrk.GetHitT2(u).GetHitClass()==1)
      counter++;

  double C0=(atrk.GetHitT2(0).GetHitX()*atrk.GetHitT2(0).GetHitX()+atrk.GetHitT2(0).GetHitY()*atrk.GetHitT2(0).GetHitY())/(atrk.GetTx()*atrk.GetHitT2(0).GetHitX()+atrk.GetTy()*atrk.GetHitT2(0).GetHitY());
  double Z0impact=atrk.GetHitT2(0).GetHitZ()-C0;


  Zimp=Z0impact;
  if((fabs(atrk.Z_at_Rmin())<13500)||((atrk.Z_at_Rmin())*atrk.GetHitT2(0).GetHitZ()<0))
    if(fabs(Z0impact)<10000)//was 5000
      cutpassed=true;

       
	     
  return cutpassed;
}








T2HalfQuarterTrkEfficiency::~T2HalfQuarterTrkEfficiency()
{
}


bool T2HalfQuarterTrkEfficiency::IsVfatMapped(int symbvfat)
{
  bool toret=false;
  map<unsigned int, VFATRegisters>::const_iterator dit = testmap->readoutIdToRegisters.find(symbvfat);
	
  if(dit != testmap->readoutIdToRegisters.end()) 
    {
      toret=true;
    }
  else 
    {
      toret=false;         
    }
  
    

  return toret;
}




//--------------------------------------------------------------------------------------------------------------------------------------
// BEGIN ALIGNMENT PART FUNCTIONS
//--------------------------------------------------------------------------------------------------------------------------------------
unsigned int T2HalfQuarterTrkEfficiency::RawtoSymb(uint32_t thedet)
{
  T2DetId converter;
  unsigned int pl=converter.plane(thedet);
  unsigned int pls=converter.planeSide(thedet);
  unsigned int ht=converter.halfTelescope(thedet);
  unsigned int arm=converter.arm(thedet);
  unsigned int symbolic=pl*2+pls+ht*10+20*arm;	  
 
  return symbolic;
}







bool lookforactivepad(edm::Handle<T2HitCollection> t2hitcoll,T1T2Track TrkCit,unsigned int m)
{
 bool flag=false;
 return flag;
}
bool lookforactivestrip(edm::Handle<T2HitCollection> t2hitcoll,T1T2Track TrkCit,unsigned int m)
{
  bool flag=false;
  return flag;
}




std::vector<double>  ResiduiRPhi(std::vector<double> vpar,double x, double y,double z) //vpar[4]
{
  //double retdrdphi[2];
   std::vector<double> retdrdphi;
 for (unsigned int i=0;i<2;i++)
   retdrdphi.push_back(0.0);

  retdrdphi[0]=(vpar[0]*z+vpar[1]-x);
  retdrdphi[1]=(vpar[2]*z+vpar[3]-y);  
  return retdrdphi;
}




bool  T2HalfQuarterTrkEfficiency::HitIsInTrackColl(T2Hit trackinghit,T1T2TrackCollection trackColl) 
{
  bool HitisAlreadyInTrk=false;
  T1T2TrackCollection::const_iterator TrkCit;
  for(TrkCit=trackColl.begin(); TrkCit!=trackColl.end(); TrkCit++){
    
    for(unsigned int u=0;u<(*TrkCit).GetHitEntries();u++)

    if (fabs(trackinghit.GetHitPhi()-(*TrkCit).GetHitT2(u).GetHitPhi())<0.001)
      if (fabs(trackinghit.GetHitR()-(*TrkCit).GetHitT2(u).GetHitR())<0.001)
	if (fabs(trackinghit.GetHitZ()-(*TrkCit).GetHitT2(u).GetHitZ())<0.001)
	  {
	    HitisAlreadyInTrk=true;
	  }
  
  }
  
  return HitisAlreadyInTrk;
}



bool T2HalfQuarterTrkEfficiency::PlaneInRightHalf(uint32_t cmsswid)
{
  bool ret=false;
  if((RawtoSymb(cmsswid)/10)==0/*SelectedHalf*/)
    {
      ret=true;
      //     if(verbosity)
      //std::cout<<"PlaneInRightHalf "<<RawtoSymb(cmsswid)<<" Accepted"<<std::endl;
    }
  //  else
  //{
  //  if(verbosity)
  //std::cout<<"SYMB PLANE="<<RawtoSymb(cmsswid)<<" WHILE SELECTED-HALF="/*<<SelectedHalf*/<<std::endl;
  //}
  
  return ret; 
}




bool T2HalfQuarterTrkEfficiency::TrackInOverlappingRegion(T1T2Track trk)
 {
   bool flag=false;
   T2SelectionCutUtils T2CutUtils;
   double phitrk=T2CutUtils.PhiFitAverage(trk.GetT2TrackHits());
   unsigned int sideref=ReferenceQuarters.at(0)%2;

   double phiminup=84.;double phimaxup=96.;
   double phimindown=264.;double phimaxdown=276.;

   //It means that the reference is not going too far from the test.
   if(UseRestrictedOverlap)
     {
       if(sideref==0)
	 {
	   phiminup=87.;
	   phimaxup=96.;
	   phimindown=264.;
	   phimaxdown=273.;
	 }
        if(sideref==1)
	 {
	   phiminup=84.;
	   phimaxup=93.;
	   phimindown=267.;
	   phimaxdown=276.;
	 }
     }

 
   
   if((phitrk>phimindown)&&(phitrk<phimaxdown))
     {
       flag=true; 
     }
   
   if((phitrk>phiminup)&&(phitrk<phimaxup))
     {
       flag=true; 
     }

   return flag;
 }

bool T2HalfQuarterTrkEfficiency::HitInOverlappingRegion(T2HitCollection::const_iterator ithit)
 {
   bool flag=false;
   if((ithit->GetHitPhi()>84.)&&(ithit->GetHitPhi()<96.))
     {
       flag=true; 
     }
    if((ithit->GetHitPhi()>264.)&&(ithit->GetHitPhi()<276.))
     {
       flag=true; 
     }
   
   return flag;
 }



double T2HalfQuarterTrkEfficiency::GetMiniumumPhiDistance(T1T2Track reftrk,T1T2TrackCollection::const_iterator thebegin_, T1T2TrackCollection::const_iterator theend_)
{
  T2SelectionCutUtils T2CutUtils;
  // double minimum=-1.0;
  std::vector<T2Hit> hits;
  hits=reftrk.GetT2TrackHits();

  double phiref=T2CutUtils.PhiFitAverage(hits);
  double fabsminphidist=500.;
  double minphidist=500.;
  for(T1T2TrackCollection::const_iterator TrkCit=thebegin_; TrkCit!=theend_; TrkCit++)
    {  
      std::vector<T2Hit> hits2; T1T2Track test=(*TrkCit);
      hits2=test.GetT2TrackHits();
      double phitrk=T2CutUtils.PhiFitAverage(hits2);
      if(fabs(phitrk-phiref)<fabsminphidist)
	{
	  minphidist= (phitrk-phiref);
	  fabsminphidist=fabs(phitrk-phiref);
	}
    }
  
  return minphidist;
}



bool T2HalfQuarterTrkEfficiency::QuarterBombarda(int numplane, int multipadplane, int multstripplane,const T2PadDigiCollection* PadDigiptr,const T2StripDigiCollection* StripDigiptr)
{

  DigiContainerIterator<T2DetId, T2PadDigi> itp;   
  DigiContainerIterator<T2DetId, T2StripDigi> its;
 
  unsigned int theref=ReferenceQuarters.at(0);
  unsigned int thetest=0;
  if(theref==0)
    thetest=1;
  if(theref==1)
    thetest=0;
  if(theref==2)
    thetest=3;
  if(theref==3)
    thetest=2;
  
  
 bool bombardainquarter= false;
 int numoffullstripplane=0;
 int numoffullpadplane=0; 
 unsigned int symb;
 //Added on 15 June 2010

    for(itp= PadDigiptr->begin(); itp!=PadDigiptr->end(); ++itp)
    {

      T2DetId mydet=(*itp).first;
      symb= mydet.plane()*2 + mydet.planeSide() + mydet.arm()*20 + mydet.halfTelescope()*10;

      
      if(thetest==(symb/10))
       {
	 int numofpadsintheplane=(*itp).second.second - (*itp).second.first ;
	 if(numofpadsintheplane > multipadplane)
	   numoffullpadplane++;
       }
     }

    for(its= StripDigiptr->begin(); its!=StripDigiptr->end(); ++its)
    {

      T2DetId mydet=(*itp).first;
      symb= mydet.plane()*2 + mydet.planeSide() + mydet.arm()*20 + mydet.halfTelescope()*10;

      if(thetest==(symb/10))
       {
	 int numofstripsintheplane=(*its).second.second - (*its).second.first ;
	 if(numofstripsintheplane > multstripplane)
	   numoffullstripplane++;
       }
     }
    
    if((numoffullstripplane>numplane)||(numoffullpadplane>numplane))
	bombardainquarter=true;

    // std::cout<<"Number of plane with more the 50 pad/Strip ON: "<<numoffullpadplane<<"/"<<numoffullstripplane<<std::endl;
 
  return bombardainquarter;
}




//--------------------------------------------------------------------------------------------------------------------------------------
// ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******
//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
// ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******
//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
// ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******
//--------------------------------------------------------------------------------------------------------------------------------------
 




 



// ------------ method called once each job just before starting event loop  ------------
void T2HalfQuarterTrkEfficiency::beginJob()//const edm::EventSetup&
{

  std::cout<<"Begin Job-1"<<std::endl;
  //char *cmsswPath = getenv("CMSSW_BASE");
 
  
  std::cout<<"Begin Job-3"<<std::endl;

  TH1::AddDirectory(kFALSE);
  //std::vector<double> alldetZ;
 
  numevent=1;
  //matrixentries=0;
  trkcounter=0;
  countegood=0;
  /*
    std::auto_ptr<TProfile> H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[10];
  std::auto_ptr<TProfile> H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_ZcutOnGeant[10];
  std::auto_ptr<TProfile> H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[10];
  std::auto_ptr<TProfile> H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[10];
  std::auto_ptr<TProfile> H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[10];
  */


  for(unsigned int y=0;y<10;y++)
    {
      HITNoiseEnt.push_back(0);
      PADNoiseEnt.push_back(0);
      STRIPNoiseEnt.push_back(0);

      TotHITNoiseEnt.push_back(0);
      TotPADNoiseEnt.push_back(0);
      TotSTRIPNoiseEnt.push_back(0);      

    }
  
   int totnumberbinfordndeta=480;
   double maxetafordndetahisto=12.0; //range assumed simmetric around eta=0

 
   char sZname2[1024];
   char sZnamehist[1024];

   for(unsigned int m=0;m<10; m++)
     {
       sprintf(sZname2, "H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult %d", m); 
       sprintf(sZnamehist, "H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult-AvgPadClu= %d", m); 
       H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

       sprintf(sZname2, "H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_56Division %d", m); 
       sprintf(sZnamehist, "H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_56Division-AvgPadClu= %d", m); 
       H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_56Division[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
    
       sprintf(sZname2, "H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult %d", m); 
       sprintf(sZnamehist, "H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult-AvgPadClu= %d", m); 
       H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
       sprintf(sZname2, "H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult %d", m); 
       sprintf(sZnamehist, "H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult-AvgPadClu= %d", m); 
       H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
       sprintf(sZname2, "H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult %d", m); 
       sprintf(sZnamehist, "H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult-AvgPadClu= %d", m); 
       H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
     
       sprintf(sZname2, "H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_ZcutOnGeant %d", m); 
       sprintf(sZnamehist, "H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_ZcutOnGeant-AvgPadClu= %d", m); 
       H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_ZcutOnGeant[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

     }


  Num_RefHalfTrkH0= std::auto_ptr<TH1F>(new TH1F("Num_RefHalfTrkH0","Num_RefHalfTrkH0",60, -0.5, 59.5));
  Num_RefHalfTrkH1= std::auto_ptr<TH1F>(new TH1F("Num_RefHalfTrkH1","Num_RefHalfTrkH1",60, -0.5, 59.5));
  Num_RefHalfTrkH2= std::auto_ptr<TH1F>(new TH1F("Num_RefHalfTrkH2","Num_RefHalfTrkH2",60, -0.5, 59.5));
  Num_RefHalfTrkH3= std::auto_ptr<TH1F>(new TH1F("Num_RefHalfTrkH3","Num_RefHalfTrkH3",60, -0.5, 59.5));
  Reference_RawTrkNumHit_NoCut = std::auto_ptr<TH1F>(new TH1F("Reference_RawTrkNumHit_NoCut","Reference_RawTrkNumHit_NoCut",60, -0.5, 59.5));
  
  AvgRefMultH0= std::auto_ptr<TH1D>(new TH1D("AvgRefMultH0","AvgRefMultH0",160, -0.5, 159.5)); 
  AvgRefMultH1= std::auto_ptr<TH1D>(new TH1D("AvgRefMultH1","AvgRefMultH1",160, -0.5, 159.5)); 
  AvgRefMultH2= std::auto_ptr<TH1D>(new TH1D("AvgRefMultH2","AvgRefMultH2",160, -0.5, 159.5)); 
  AvgRefMultH3= std::auto_ptr<TH1D>(new TH1D("AvgRefMultH3","AvgRefMultH3",160, -0.5, 159.5)); 


  NumRefTrk_vs_RefAvgHitMultH0= std::auto_ptr<TProfile>(new TProfile("NumRefTrk_vs_RefAvgHitMultH0","NumRefTrk_vs_RefAvgHitMultH0",20,-0.5,19.5,""));
  NumRefTrk_vs_RefAvgHitMultH1= std::auto_ptr<TProfile>(new TProfile("NumRefTrk_vs_RefAvgHitMultH1","NumRefTrk_vs_RefAvgHitMultH1",20,-0.5,19.5,""));
  NumRefTrk_vs_RefAvgHitMultH2= std::auto_ptr<TProfile>(new TProfile("NumRefTrk_vs_RefAvgHitMultH2","NumRefTrk_vs_RefAvgHitMultH2",20,-0.5,19.5,""));
  NumRefTrk_vs_RefAvgHitMultH3= std::auto_ptr<TProfile>(new TProfile("NumRefTrk_vs_RefAvgHitMultH3","NumRefTrk_vs_RefAvgHitMultH3",20,-0.5,19.5,""));
  UnassHitFractH0= std::auto_ptr<TH1D>(new TH1D("UnassHitFractH0","UnassHitFractH0",100, -0.005, 10.05));
  UnassHitFractH1= std::auto_ptr<TH1D>(new TH1D("UnassHitFractH1","UnassHitFractH1",100, -0.005, 10.05));
  UnassHitFractH2= std::auto_ptr<TH1D>(new TH1D("UnassHitFractH2","UnassHitFractH2",100, -0.005, 10.05));
  UnassHitFractH3= std::auto_ptr<TH1D>(new TH1D("UnassHitFractH3","UnassHitFractH3",100, -0.005, 10.05));
  H0RefHitActivePlane= std::auto_ptr<TH1D>(new TH1D("H0RefHitActivePlane","H0RefHitActivePlane",10,-0.5,9.5));





  CluStripentriesGoodHit = std::auto_ptr<TH1F>(new TH1F("CluStripentriesGoodHit","Strip cluster-size (Cl1 hit, noise excluded from DR-D#phi cut)",20, -0.5, 19.5));
  CluStripentriesGoodHit->SetDirectory(0); 
  CluPadentriesGoodHit = std::auto_ptr<TH1F>(new TH1F("CluPadentriesGoodHit","Pad cluster-size (Cl1 hit, noise excluded from DR-D#phi cut)",20, -0.5, 19.5));
  CluPadentriesGoodHit->SetDirectory(0); 

 
  AllHitR= std::auto_ptr<TH1F>(new TH1F("AllHitR","AllHitR",120, 30., 150.));
  AllHitR->SetDirectory(0); 
 
  chiRredREF = std::auto_ptr<TH1F>(new TH1F("chiRredREF","Chi2ridotto riferimento",100, 0., 50.));
  clusterstripentries = std::auto_ptr<TH1F>(new TH1F("clusterstripentries","Number of strips in strip-clusters",20, -0.5, 19.5));
  clusterstripentries->SetDirectory(0); 
  clusterpadentries = std::auto_ptr<TH1F>(new TH1F("clusterpadentries","Number of pads in pad-clusters",20, -0.5, 19.5));
  clusterpadentries->SetDirectory(0); 
  diffphiGUNHIT= std::auto_ptr<TH1F>(new TH1F("diffphiGUNHIT","#Phi Generated particle - #Phi Reconstructed Hit",120, -30, 30));
  diffphiGUNHIT->SetDirectory(0); 
  RLocHRecoH= std::auto_ptr<TH1F>(new  TH1F("RLocHRecoH"," Reconstructed Hit R - GEANT4 Hit R",100,-1.5, 1.5));
  RLocHRecoH->SetDirectory(0); 
  diffRCluHit= std::auto_ptr<TH1F>(new  TH1F("diffRCluHit"," Cluster R - GEANT4 Hit R",100,-1.5, 1.5));
  diffRCluHit->SetDirectory(0); 
  Trketa= std::auto_ptr<TH1F>(new TH1F("Trketa","Reconstructed Track #eta (no cut used)",300,1.0,8.0));
  Trketa->SetDirectory(0); 
  Trkphi= std::auto_ptr<TH1F>(new TH1F("Trkphi","Reconstructed Track #phi (no cut used)",361,0,360));
  Trkphi->SetDirectory(0); 

    
  TrkphiALL= std::auto_ptr<TH1F>(new TH1F("TrkphiALL","Reconstructed Track #phi, ALL quarters",361,0,360));
  TrketaALL= std::auto_ptr<TH1F>(new TH1F("TrketaALL","Reconstructed Track #eta ",200,-20,20));
  NumhitinTrackALL= std::auto_ptr<TH1F>(new TH1F("NumhitinTrackALL","Number of hits forming the track",21,-0.5,20.5));
  TrkQuarterIdALL= std::auto_ptr<TH1F>(new TH1F("TrkQuarterIdALL","Id of the first plane in Track",40,-0.5,39.5));

  

  NumTrackALL_ONEQuarter= std::auto_ptr<TH1F>(new TH1F("NumTrackALL_ONEQuarter","NumTrackALL_ONEQuarter",151,-0.5,150.5));
  NumTrackALL_ONEQuarter->SetDirectory(0); 

  NumTrackALL_EveryQuarter= std::auto_ptr<TH1F>(new TH1F("NumTrackALL_EveryQuarter","NumTrackALL_EveryQuarter",151,-0.5,150.5));
  NumTrackALL_EveryQuarter->SetDirectory(0); 

  

  DPhiGoodTrk= std::auto_ptr<TH1F>(new TH1F("DPhiGoodTrk","#phi Track - #phi GEANT4 Track (Only primary tracks are considered)",51,-25.5,25.5));
  DPhiGoodTrk->SetDirectory(0); 
  DEtaGoodTrk= std::auto_ptr<TH1F>(new TH1F("DEtaGoodTrk","#eta Track - #phi GEANT4 Track (Only primary tracks are considered)",200,-1.0,1.0));
  DEtaGoodTrk->SetDirectory(0); 
  diffphiCluGun = std::auto_ptr<TH1F>(new TH1F("diffphiCluGun","Cluster #Phi - #Phi Generated particle",81, -40.5, 40.5));
  diffphiCluGun->SetDirectory(0); 
  diffphiCluHit= std::auto_ptr<TH1F>(new TH1F("diffphiCluHit","Cluster #Phi - GEANT4 Hit #Phi",50, -5, 5));
  diffphiCluHit->SetDirectory(0); 
  Trketagood= std::auto_ptr<TH1F>(new TH1F("Trketagood","Reconstructed Track #eta (#chi^{2}_{p} cut)",200,1.0,8.0));
  Trketagood->SetDirectory(0); 
  //Trkphigood= std::auto_ptr<TH1F>(new TH1F("Trkphigood","Reconstructed Track #phi (#chi^{2}_{p} cut)",361,0,360));
  Trkphigood= std::auto_ptr<TH1F>(new TH1F("Trkphigood","Reconstructed Track #phi",361,0,360));
  Trkphigood->SetDirectory(0); 

  Chi2RProb = std::auto_ptr<TH1F>(new TH1F("Chi2RProb","Radial #chi^{2} probability",15000,0.0,1.05));
  Chi2RProb->SetDirectory(0); 
  Chi2PhiProb = std::auto_ptr<TH1F>(new TH1F("Chi2PhiProb","Azimuthal #chi^{2} probability",15000,0.0,1.05));
  Chi2PhiProb->SetDirectory(0);

  IstoDisp= std::auto_ptr<TH1F>(new TH1F("IstoDisp","Displ.",100,0.0,10.0));
  IstoDisp->SetDirectory(0);
   
  DigiPadOccupancy= std::auto_ptr<TProfile>(new TProfile("DigiPadOccupancy","Plane Pad occupancy",10,-0.5,9.5,""));
  DigiPadOccupancy->GetYaxis()->SetTitle("Occupancy (%)");
  DigiPadOccupancy->GetXaxis()->SetTitle("Detector plane");
  DigiPadOccupancy->SetDirectory(0);

  DigiStripOccupancy= std::auto_ptr<TProfile>(new TProfile("DigiStripOccupancy","Plane Strip occupancy",10,-0.5,9.5,""));
  DigiStripOccupancy->GetYaxis()->SetTitle("Occupancy (%)");
  DigiStripOccupancy->GetXaxis()->SetTitle("Detector plane");
  DigiStripOccupancy->SetDirectory(0);

  Class1HitPadStripCLSCorrel= std::auto_ptr<TH2D>(new TH2D("Class1HitPadStripCLSCorrel","Pad-Strip Cluster size correlation (Class1 Hit)",60,-0.5,59.5,60,-0.5,59.5));
  Class1HitPadStripCLSCorrel->GetYaxis()->SetTitle("Pad CL Size");
  Class1HitPadStripCLSCorrel->GetXaxis()->SetTitle("Strip CL Size"); 
  Class1HitPadStripCLSCorrel->SetDirectory(0);  


  SingleParticleEfficiencyCut = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyCut","# Tracks reconstructed vs generated particle #eta (all cuts)",34,4.35,7.75,""));
  SingleParticleEfficiencyCut->SetDirectory(0); 
  SingleParticleEfficiency = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiency","# Tracks reconstructed vs generated particle #eta (all track)",34,4.35,7.75,""));
  SingleParticleEfficiency->SetDirectory(0); 
  MultiParticleEfficiencyCut= std::auto_ptr<TProfile>(new TProfile("MultiParticleEfficiencyCut","<# Tracks> vs # Charged Particles (all track cuts)",12,-0.5,11.5,""));
  MultiParticleEfficiencyCut->SetDirectory(0); 
  MultiParticleEfficiencyCutNorm= std::auto_ptr<TProfile>(new TProfile("MultiParticleEfficiencyCutNorm","<# Tracks / # Charged Particles > vs # Charged Particles (all track cuts)",12,-0.5,11.5,"")); 
  MultiParticleEfficiencyCutNorm->SetDirectory(0); 
   

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  


  QuarterQuarterDXEntry= std::auto_ptr<TH1D>(new TH1D("QuarterQuarterDXEntry","Ref. Trk Entry X - Test Trk Entry X",60,-15,15));
  QuarterQuarterDXExit= std::auto_ptr<TH1D>(new TH1D("QuarterQuarterDXExit","Ref. Trk Entry X - Test Trk Exit X",60,-15,15));
  QuarterQuarterDYEntry= std::auto_ptr<TH1D>(new TH1D("QuarterQuarterDYEntry","Ref. Trk Entry X - Test Trk Entry Y",60,-15,15));
  QuarterQuarterDYExit= std::auto_ptr<TH1D>(new TH1D("QuarterQuarterDYExit","Ref. Trk Entry X - Test Trk Exit Y",60,-15,15));



  
  DYTrkHit_HRefHTestClose= std::auto_ptr<TH1D>(new TH1D("DYTrkHit_HRefHTestClose","Track #Phi close to the ref. in the tested quarter",60,-15,15)) ;  
  DXTrkHit_HRefHTestClose = std::auto_ptr<TH1D>(new TH1D("DXTrkHit_HRefHTestClose","Track #Phi close to the ref. in the tested quarter",60,-15,15)) ; 
  DRTrkHit_HRefHTestClose= std::auto_ptr<TH1D>(new TH1D("DRTrkHit_HRefHTestClose","Track #Phi close to the ref. in the tested quarter",40,-10,10)) ;  
  DPhiTrkHit_HRefHTestClose= std::auto_ptr<TH1D>(new TH1D("DPhiTrkHit_HRefHTestClose","Track #Phi close to the ref. in the tested quarter",20,-10,10)) ; 

						 
  TrkRinRefQuarter=std::auto_ptr<TH1D>(new TH1D("TrkRinRefQuarter","TrkRinRefQuarter",60,30,150)) ;
  TrkPhiinRefQuarter=std::auto_ptr<TH1D>(new TH1D("TrkPhiinRefQuarter","TrkPhiinRefQuarter",120,-0.5,359.5)) ;
  TrkXinRefQuarter=std::auto_ptr<TH1D>(new TH1D("TrkXinRefQuarter","TrkXinRefQuarter",60,-15,15)) ;
  TrkYinRefQuarter=std::auto_ptr<TH1D>(new TH1D("TrkYinRefQuarter","TrkYinRefQuarter",60,-150,150)) ;

  ReferenceTrackXY=  std::auto_ptr<TH2F>(new TH2F("ReferenceTrackXY","ReferenceTrackXY",300,-150,150.,300,-150.,150.));



  TrkPhiinTestedQuarterCloseToRef= std::auto_ptr<TH1D>(new TH1D("TrkPhiinTestedQuarterCloseToRef","Track #Phi close to the ref. in the tested quarter",360,-0.5,359.5)) ; 
  TrkRinTestedQuarterCloseToRef=  std::auto_ptr<TH1D>(new TH1D("TrkRinTestedQuarterCloseToRef","TrkRinTestedQuarterCloseToRef",60,30,150.)) ;
  TrkXinTestedQuarterCloseToRef=  std::auto_ptr<TH1D>(new TH1D("TrkXinTestedQuarterCloseToRef","TrkXinTestedQuarterCloseToRef",60,-30.,30.)) ;
  TrkYinTestedQuarterCloseToRef=  std::auto_ptr<TH1D>(new TH1D("TrkYinTestedQuarterCloseToRef","TrkYinTestedQuarterCloseToRef",300,-150,150)) ;
  DPhiTrk_HRefHTest2CloseToRef= std::auto_ptr<TH1D>(new TH1D("DPhiTrk_HRefHTest2CloseToRef","",30,-15,15)) ;


  DRTrk_HRefHTest= std::auto_ptr<TH1F>(new TH1F("DRTrk_HRefHTest","DRTrk_HRefHTest",30,-15,15)) ;



  DYTrkHit_HRefHTestCloseVsZTest= std::auto_ptr<TProfile>(new TProfile("DYTrkHit_HRefHTestCloseVsZTest","DYTrkHit_HRefHTestCloseVsZTest",400,13800,14200,"")) ;
  DXTrkHit_HRefHTestCloseVsZTest= std::auto_ptr<TProfile>(new TProfile("DXTrkHit_HRefHTestCloseVsZTest","DXTrkHit_HRefHTestCloseVsZTest",400,13800,14200,"")) ;  
  DRTrkHit_HRefHTestCloseVsZTest= std::auto_ptr<TProfile>(new TProfile("DRTrkHit_HRefHTestCloseVsZTest","DRTrkHit_HRefHTestCloseVsZTest",400,13800,14200,"")) ; 
  DPhiTrkHit_HRefHTestCloseVsZTest= std::auto_ptr<TProfile>(new TProfile("DPhiTrkHit_HRefHTestCloseVsZTest","DPhiTrkHit_HRefHTestCloseVsZTest",400,13800,14200,"")) ;
  
  RefTrackXall=  std::auto_ptr<TH1D>(new TH1D("RefTrackXall","RefTrackXall",60,-30.,30.)) ;  

  HowManyPadCluCloseToRef_ArmXFar=std::auto_ptr<TH1F>(new TH1F("HowManyPadCluCloseToRef_ArmXFar","# Pad Clust. in X-Far close to the ref. track",30,-0.5,29.5)) ; 
  

  HowManyPadCluCloseToRef_ArmXNear=std::auto_ptr<TH1F>(new TH1F("HowManyPadCluCloseToRef_ArmXNear","# Pad Clust. in X-Near close to the ref. track",30,-0.5,29.5)) ;

  

  HowManyPadCluCloseToRef_PluNear =  std::auto_ptr<TH1F>(new TH1F("HowManyPadCluCloseToRef_PluNear","# Pad Clust. in +Near close to the ref. track",30,-0.5,29.5));
   
  HowManyPadCluCloseToRef_PluFar=  std::auto_ptr<TH1F>(new TH1F("HowManyPadCluCloseToRef_PluFar","# Pad Clust. in +Far close to the ref. track",30,-0.5,29.5));

  DPhiTrk_HRefHTest=  std::auto_ptr<TH1F>(new TH1F("DPhiTrk_HRefHTest","Overlap Track |#Phi reference - #Phi Test|",361,-0.5,360.5));
  DPhiDRTrk_HRefHTest=  std::auto_ptr<TH2F>(new TH2F("DPhiDRTrk_HRefHTest","Matching Track |#Phi reference - #Phi Test|,|R reference - R Test|",361,-0.5,360.5,146,-1.,145.));

  HitInOVERLTrk_ALL_Phi_Residual=  std::auto_ptr<TH1F>(new TH1F("HitInOVERLTrk_ALL_Phi_Residual","Cumulative #Delta #phi hits contained in overlap-tracks",40,-20,20));
  HitInOVERLTrk_ALL_R_Residual=  std::auto_ptr<TH1F>(new TH1F("HitInOVERLTrk_ALL_R_Residual","Cumulative #Delta R hits contained in overlap-tracks",600,-30,30));

  HitH0vsHitH1_OvrlTrk=  std::auto_ptr<TH2F>(new TH2F("HitH0vsHitH1_OvrlTrk","X=# H0 Hit   Y=# H1 Hit for a common built Ovrl Trk, clean event",12,-0.5,11.5,12,-0.5,11.5));
  
  HitH0vsHitH1_OvrlTrk->SetXTitle("H0");
  HitH0vsHitH1_OvrlTrk->SetYTitle("H1");

  ScatterPlotPhi_ROverlaps =  std::auto_ptr<TH2F>(new TH2F("ScatterPlotPhi_ROverlaps","y=R x=#Phi of the overlapping tracks",360,-0.5,359.5,220,40,150));
  ScatterPlotPhi_ROverlaps->SetXTitle("#Phi");
  ScatterPlotPhi_ROverlaps->SetXTitle("R");

  ScatterPlotPhi_ROverlapsH0 =  std::auto_ptr<TH2F>(new TH2F("ScatterPlotPhi_ROverlapsH0","y=R x=#Phi of the overlapping tracks H0",360,-0.5,359.5,220,40,150));
  ScatterPlotPhi_ROverlapsH0->SetXTitle("#Phi");
  ScatterPlotPhi_ROverlapsH0->SetXTitle("R");

  H0TrackPhiALL=  std::auto_ptr<TH1F>(new TH1F("H0TrackPhiALL","ALL Trk #phi in H0",360,-0.5,359.5));
  H0TrackPhiALL->SetXTitle("#Phi");
  
  ReferenceZImpVsEta2=  std::auto_ptr<TH2F>(new TH2F("ReferenceZImpVsEta2","y=R x=#Phi of the overlapping tracks H0",500,-10000,10000,20,5.2,6.6));

  ZIMPdistrRefH0=  std::auto_ptr<TH1F>(new TH1F("ZIMPdistrRefH0","ZIMPdistrRefH0",500,-10000,10000));
  ZIMPdistrRefH1=  std::auto_ptr<TH1F>(new TH1F("ZIMPdistrRefH1","ZIMPdistrRefH1",500,-10000,10000));
  ZIMPdistrRefH2=  std::auto_ptr<TH1F>(new TH1F("ZIMPdistrRefH2","ZIMPdistrRefH2",500,-10000,10000));
  ZIMPdistrRefH3=  std::auto_ptr<TH1F>(new TH1F("ZIMPdistrRefH3","ZIMPdistrRefH3",500,-10000,10000));

  ScatterPlotPhi_R_fullPhiH0=  std::auto_ptr<TH2F>(new TH2F("ScatterPlotPhi_R_fullPhiH0","y=R x=#Phi of the overlapping tracks H0",360,-0.5,359.5,220,40,150));
  ScatterPlotPhi_R_fullPhiH0->SetXTitle("#Phi");
  ScatterPlotPhi_R_fullPhiH0  ->SetXTitle("R");
  
  ScatterPlotPhi_ROverlapsH1 =  std::auto_ptr<TH2F>(new TH2F("ScatterPlotPhi_ROverlapsH1","y=R x=#Phi of the overlapping tracks H1",360,-0.5,359.5,220,40,150));
  ScatterPlotPhi_ROverlapsH1->SetXTitle("#Phi");
  ScatterPlotPhi_ROverlapsH1->SetXTitle("R");
  
  ScatterPlotPhi_ROverlapsH2 =  std::auto_ptr<TH2F>(new TH2F("ScatterPlotPhi_ROverlapsH2","y=R x=#Phi of the overlapping tracks H2",360,-0.5,359.5,220,40,150));
  ScatterPlotPhi_ROverlapsH2->SetXTitle("#Phi");
  ScatterPlotPhi_ROverlapsH2->SetXTitle("R");

  ScatterPlotPhi_ROverlapsH3 =  std::auto_ptr<TH2F>(new TH2F("ScatterPlotPhi_ROverlapsH3","y=R x=#Phi of the overlapping tracks H3",360,-0.5,359.5,220,40,150));
  ScatterPlotPhi_ROverlapsH3->SetXTitle("#Phi");
  ScatterPlotPhi_ROverlapsH3->SetXTitle("R");

  DPhiMinRefTrkVsTestQuarterTrkArm0=  std::auto_ptr<TH1F>(new TH1F("DPhiMinRefTrkVsTestQuarterTrkArm0","Minimum #Phi distance (arm-0) of the overlapping tracks H1",720,-359.5,359.5));
  DPhiMinRefTrkVsTestQuarterTrkArm0->SetXTitle("#Phi");

  ReasonOfMissingTrkArm0PluNear = std::auto_ptr<TH1F>(new TH1F("ReasonOfMissingTrkArm0PluNear","Reason of failure efficiency for +Near: 0=SmallNumHit 1=SmallNumHitCloseToReference 2=AlgoInefficient",3, -0.5, 2.5)); 

  ReasonOfMissingTrkArm0PluFar = std::auto_ptr<TH1F>(new TH1F("ReasonOfMissingTrkArm0PluFar","Reason of failure efficiency for +Far: 0=SmallNumHit 1=SmallNumHitCloseToReference 2=AlgoInefficient",3, -0.5, 2.5));


 
  ReasonOfMissingTrkArmXFar=std::auto_ptr<TH1F>(new TH1F("ReasonOfMissingTrkArmXFar","Reason of failure efficiency for X-Far: 0=SmallNumHit 1=SmallNumHitCloseToReference 2=AlgoInefficient",3, -0.5, 2.5))  ;  
  ReasonOfMissingTrkArmXNear=std::auto_ptr<TH1F>(new TH1F("ReasonOfMissingTrkArmXNear","Reason of failure efficiency for X-Near: 0=SmallNumHit 1=SmallNumHitCloseToReference 2=AlgoInefficient",3, -0.5, 2.5)) ;



  
  NumCl1HitInOvRegionArm0_Up= std::auto_ptr<TH1D>(new TH1D("NumCl1HitInOvRegionArm0_Up","NumCl1HitInOvRegionArm0_Up",200, -0.5, 199.5));
  NumCl1HitInOvRegionArm0_Down= std::auto_ptr<TH1D>(new TH1D("NumCl1HitInOvRegionArm0_Down","NumCl1HitInOvRegionArm0_Down",200, -0.5, 199.5)); 


  HXTrkEffi=  std::auto_ptr<TProfile>(new TProfile("HXTrkEffi","<# matching Tracks HTest / # Tracks HRef > vs # Tracks HRef. ",12,-0.5,11.5,"")); 
  HXTrkEffivsR=  std::auto_ptr<TProfile>(new TProfile("HXTrkEffivsR","HX track efficiency vs R ",4,40,150,"")); 

  HXTrkEffi_AxisRawTrk=  std::auto_ptr<TProfile>(new TProfile("HXTrkEffi_AxisRawTrk","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",15,-0.5,15.5,"")); 



  H0TrkEffi_AxisAvgMult=  std::auto_ptr<TProfile>(new TProfile("H0TrkEffi_AxisAvgMult","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 

  H1TrkEffi_AxisAvgMult=  std::auto_ptr<TProfile>(new TProfile("H1TrkEffi_AxisAvgMult","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 

  H2TrkEffi_AxisAvgMult=  std::auto_ptr<TProfile>(new TProfile("H2TrkEffi_AxisAvgMult","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 

  H3TrkEffi_AxisAvgMult=  std::auto_ptr<TProfile>(new TProfile("H3TrkEffi_AxisAvgMult","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 


  H0TrkEffi_AxisAvgMult_EtaCut=  std::auto_ptr<TProfile>(new TProfile("H0TrkEffi_AxisAvgMult_EtaCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 
  H1TrkEffi_AxisAvgMult_EtaCut=  std::auto_ptr<TProfile>(new TProfile("H1TrkEffi_AxisAvgMult_EtaCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 
  H2TrkEffi_AxisAvgMult_EtaCut=  std::auto_ptr<TProfile>(new TProfile("H2TrkEffi_AxisAvgMult_EtaCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 
  H3TrkEffi_AxisAvgMult_EtaCut=  std::auto_ptr<TProfile>(new TProfile("H3TrkEffi_AxisAvgMult_EtaCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 


  H0TrkEffi_AxisAvgMult_ZImpactCut=  std::auto_ptr<TProfile>(new TProfile("H0TrkEffi_AxisAvgMult_ZImpactCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 
  H1TrkEffi_AxisAvgMult_ZImpactCut=  std::auto_ptr<TProfile>(new TProfile("H1TrkEffi_AxisAvgMult_ZImpactCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 
  H2TrkEffi_AxisAvgMult_ZImpactCut=  std::auto_ptr<TProfile>(new TProfile("H2TrkEffi_AxisAvgMult_ZImpactCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 
  H3TrkEffi_AxisAvgMult_ZImpactCut=  std::auto_ptr<TProfile>(new TProfile("H3TrkEffi_AxisAvgMult_ZImpactCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 




  H0TrkEffi_AxisRawTrk=  std::auto_ptr<TProfile>(new TProfile("H0TrkEffi_AxisRawTrk","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",16,-0.5,15.5,"")); 

  H1TrkEffi_AxisRawTrk=  std::auto_ptr<TProfile>(new TProfile("H1TrkEffi_AxisRawTrk","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",16,-0.5,15.5,"")); 

  H2TrkEffi_AxisRawTrk=  std::auto_ptr<TProfile>(new TProfile("H2TrkEffi_AxisRawTrk","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",16,-0.5,15.5,"")); 

  H3TrkEffi_AxisRawTrk=  std::auto_ptr<TProfile>(new TProfile("H3TrkEffi_AxisRawTrk","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",16,-0.5,15.5,"")); 


  
  H0TrkEffi_AxisRawTrk_EtaCut=  std::auto_ptr<TProfile>(new TProfile("H0TrkEffi_AxisRawTrk_EtaCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",16,-0.5,15.5,"")); 

  H1TrkEffi_AxisRawTrk_EtaCut=  std::auto_ptr<TProfile>(new TProfile("H1TrkEffi_AxisRawTrk_EtaCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",16,-0.5,15.5,"")); 

  H2TrkEffi_AxisRawTrk_EtaCut=  std::auto_ptr<TProfile>(new TProfile("H2TrkEffi_AxisRawTrk_EtaCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",16,-0.5,15.5,"")); 

  H3TrkEffi_AxisRawTrk_EtaCut=  std::auto_ptr<TProfile>(new TProfile("H3TrkEffi_AxisRawTrk_EtaCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",16,-0.5,15.5,"")); 










  HXTrkEffi_Yplus=  std::auto_ptr<TProfile>(new TProfile("HXTrkEffi_Yplus","<# matching Tracks H0 / # Tracks H1 > vs # Tracks HRef. (Y>0)",12,-0.5,11.5,"")); 
  HXTrkEffi_Yminus=  std::auto_ptr<TProfile>(new TProfile("HXTrkEffi_Yminus","<# matching Tracks H0 / # Tracks H1 > vs # Tracks HRef. (Y<0)",12,-0.5,11.5,"")); 

  HXTrkEffivsRplus=  std::auto_ptr<TProfile>(new TProfile("HXTrkEffivsRplus","HX track efficiency vs R. (Y>0)",4,40,150,"")); 
  HXTrkEffivsRminus=  std::auto_ptr<TProfile>(new TProfile("HXTrkEffivsRminus","HX track efficiency vs R. (Y<0)",4,40,150,""));

  HX_ClusterPhi= std::auto_ptr<TH1D>(new TH1D("HX_ClusterPhi","HX_ClusterPhi",120, -0.5, 359.5));
  
  
   



  H0TrkEffi =  std::auto_ptr<TProfile>(new TProfile("H0TrkEffi","<# matching Tracks H0 / # Tracks H0 > vs # Tracks H0. ",12,-0.5,11.5,""));
  H1TrkEffi =  std::auto_ptr<TProfile>(new TProfile("H1TrkEffi","<# matching Tracks H1 / # Tracks H1> vs # Tracks  H1 ",12,-0.5,11.5,""));
  H2TrkEffi =  std::auto_ptr<TProfile>(new TProfile("H2TrkEffi","<# matching Tracks H2 / # Tracks H2> vs # Tracks  H2 ",12,-0.5,11.5,""));
  H3TrkEffi=  std::auto_ptr<TProfile>(new TProfile("H3TrkEffi","<# matching Tracks H3 / # Tracks H3> vs # Tracks  H3 ",12,-0.5,11.5,""));

  H0TrkEffi->SetDirectory(0);  H0TrkEffi->SetXTitle("# Tracks H1 ");
  H1TrkEffi->SetDirectory(0); H1TrkEffi->SetXTitle("# Tracks H0 ");
  H2TrkEffi->SetDirectory(0); H2TrkEffi->SetXTitle("# Tracks H3 ");
  H3TrkEffi->SetDirectory(0); H3TrkEffi->SetXTitle("# Tracks H2 ");

H0TrkEffivsR=  std::auto_ptr<TProfile>(new TProfile("H0TrkEffivsR","H0 single ov. track efficiency vs R ",4,40,150,""));
H0TrkEffivsR->SetXTitle("<R> mm ");
H1TrkEffivsR=  std::auto_ptr<TProfile>(new TProfile("H1TrkEffivsR","H1 single ov. track efficiency vs R. ",4,40,150,""));
H1TrkEffivsR->SetXTitle("<R> mm ");

 H0TrkEffi_Yplus =  std::auto_ptr<TProfile>(new TProfile("H0TrkEffi_Yplus","<# matching Tracks H0 / # Tracks H1 > vs # Tracks H1. (Y>0)",12,-0.5,11.5,""));
 H1TrkEffi_Yplus =  std::auto_ptr<TProfile>(new TProfile("H1TrkEffi_Yplus","<# matching Tracks H1 / # Tracks H0> vs # Tracks  H0 (Y>0)",12,-0.5,11.5,""));
 H2TrkEffi_Yplus =  std::auto_ptr<TProfile>(new TProfile("H2TrkEffi_Yplus","<# matching Tracks H2 / # Tracks H3> vs # Tracks  H3 (Y>0)",12,-0.5,11.5,""));
 H3TrkEffi_Yplus =  std::auto_ptr<TProfile>(new TProfile("H3TrkEffi_Yplus","<# matching Tracks H3 / # Tracks H2> vs # Tracks  H2 (Y>0)",12,-0.5,11.5,""));
  
 H0TrkEffi_Yminus =  std::auto_ptr<TProfile>(new TProfile("H0TrkEffi_Yminus","<# matching Tracks H0 / # Tracks H0> vs # Tracks H1 (Y<0)",12,-0.5,11.5,""))  ;
 H1TrkEffi_Yminus =  std::auto_ptr<TProfile>(new TProfile("H1TrkEffi_Yminus","<# matching Tracks H0 / # Tracks H0> vs # Tracks H0 (Y<0)",12,-0.5,11.5,""))   ;
 H2TrkEffi_Yminus = std::auto_ptr<TProfile>(new TProfile("H2TrkEffi_Yminus","<# matching Tracks H0 / # Tracks H0> vs # Tracks H3 (Y<0)",12,-0.5,11.5,""))  ;
 H3TrkEffi_Yminus = std::auto_ptr<TProfile>(new TProfile("H3TrkEffi_Yminus","<# matching Tracks H0 / # Tracks H0> vs # Tracks H2 (Y<0)",12,-0.5,11.5,""))  ;

H0TrkEffi_Yplus->SetDirectory(0);
H0TrkEffi_Yplus->SetXTitle("# Tracks H1 (Y>0)");

H1TrkEffi_Yplus->SetDirectory(0);
H1TrkEffi_Yplus->SetXTitle("# Tracks H0 (Y>0)");

H2TrkEffi_Yplus->SetDirectory(0);
H2TrkEffi_Yplus->SetXTitle("# Tracks H3 (Y>0)");

H3TrkEffi_Yplus->SetDirectory(0);
H3TrkEffi_Yplus->SetXTitle("# Tracks H2 (Y>0)");


H0TrkEffi_Yminus->SetDirectory(0);
H0TrkEffi_Yminus->SetXTitle("# Tracks H1 (Y<0)");

H1TrkEffi_Yminus->SetDirectory(0);
H1TrkEffi_Yminus->SetXTitle("# Tracks H0 (Y<0)");

H2TrkEffi_Yminus->SetDirectory(0);
H2TrkEffi_Yminus->SetXTitle("# Tracks H3 (Y<0)");

H3TrkEffi_Yminus->SetDirectory(0);
H3TrkEffi_Yminus->SetXTitle("# Tracks H2 (Y<0)");

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      

 }







// ------------ method called once each job just after ending the event loop  ------------

void T2HalfQuarterTrkEfficiency::endJob()
{
 
  std::cout<<"End Job. Number of non-corrupted event: "<<countegood<<std::endl;



  std::cout<<"A"<<std::endl;
 /* ----------- Save the output ------------ */
 /* ----------- Save the output ------------ */
 /* ----------- Save the output ------------ */


  
TFile *f = TFile::Open(outputFileName.c_str(), "recreate");
 if( !f || !f->IsWritable() ){
   std::cout << "Output file not opened correctly !!" << std::endl;
 }
 


 gDirectory = f->mkdir("Global");



 /*------------------------------------------------- ----------- Efficiency  ------------ */
 /*------------------------------------------------- ----------- Efficiency ------------ */
 /*------------------------------------------------- ----------- Efficiency ------------ */

 std::cout<<"B"<<std::endl;

 gDirectory = f->mkdir("Efficiency");


 /*
 ReferenceTrackXY->Write("");   
 QuarterQuarterDYExit->Write("");   
 QuarterQuarterDYEntry->Write("");   
 QuarterQuarterDXExit->Write("");   
 QuarterQuarterDXEntry->Write("");  
 
 std::cout<<"C"<<std::endl;

 DYTrkHit_HRefHTestClose->Write("");  
 DXTrkHit_HRefHTestClose->Write(""); 
 DRTrkHit_HRefHTestClose->Write("");     
 DPhiTrkHit_HRefHTestClose->Write("");     
 DRTrk_HRefHTest->Write("");     

 std::cout<<"D"<<std::endl;
 
 RefTrackXall->Write("");   
 DPhiTrkHit_HRefHTestCloseVsZTest->Write("");  
 DRTrkHit_HRefHTestCloseVsZTest->Write("");  
 DXTrkHit_HRefHTestCloseVsZTest->Write("");  
 DYTrkHit_HRefHTestCloseVsZTest->Write(""); 

 std::cout<<"E"<<std::endl;

 TrkPhiinTestedQuarterCloseToRef->Write(""); 
 TrkRinTestedQuarterCloseToRef->Write(""); 
 TrkXinTestedQuarterCloseToRef->Write(""); TrkYinTestedQuarterCloseToRef->Write(""); DPhiTrk_HRefHTest2CloseToRef->Write("");
 TrkPhiinRefQuarter->Write("");
 TrkXinRefQuarter->Write(""); TrkYinRefQuarter->Write("");

 std::cout<<"F"<<std::endl;

 DPhiTrkHit_HRefHTestCloseVsZTest->Write("");
 DRTrkHit_HRefHTestCloseVsZTest->Write("");
 DXTrkHit_HRefHTestCloseVsZTest->Write("");
 DYTrkHit_HRefHTestCloseVsZTest->Write("");
  
 std::cout<<"G"<<std::endl;
 */ 
 TrkRinRefQuarter->Write(""); 
 DPhiTrk_HRefHTest->Write("");
 DRTrk_HRefHTest->Write("");
 ReferenceZImpVsEta2->Write("");
 Num_RefHalfTrkH0->Write("");
 Num_RefHalfTrkH1->Write("");
 Num_RefHalfTrkH2->Write("");
 Num_RefHalfTrkH3->Write("");
 
Reference_RawTrkNumHit_NoCut->Write("");

 
H0RefHitActivePlane->Write("");
UnassHitFractH0->Write("");
NumRefTrk_vs_RefAvgHitMultH0->Write("");
UnassHitFractH1->Write("");
NumRefTrk_vs_RefAvgHitMultH1->Write("");
UnassHitFractH2->Write("");
NumRefTrk_vs_RefAvgHitMultH2->Write("");
UnassHitFractH3->Write("");
NumRefTrk_vs_RefAvgHitMultH3->Write("");

AvgRefMultH0->Write("");
AvgRefMultH1->Write("");
AvgRefMultH2->Write("");
AvgRefMultH3->Write("");

 gDirectory = f->mkdir("BinnedPrimaryEfficiency");
 for(unsigned int m=0;m<10; m++)
    {
      H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]->Write();
      H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_56Division[m]->Write();
      H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_ZcutOnGeant[m]->Write();
      H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]->Write();
      H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]->Write();
      H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]->Write();
    }


 std::cout<<"I"<<std::endl;
 

 HXTrkEffi->Write("");
 HXTrkEffivsR->Write("");
 HXTrkEffi_AxisRawTrk->Write("");
 //std::cout<<"L"<<std::endl;

 H0TrkEffi->Write("");
 H0TrkEffi_AxisRawTrk->Write("");
 // H0TrkEffivsR->Write("");

 H1TrkEffi->Write("");
 H1TrkEffi_AxisRawTrk->Write("");
 // H1TrkEffivsR->Write("");

 H2TrkEffi->Write("");
 H2TrkEffi_AxisRawTrk->Write("");
 // H2TrkEffivsR->Write("");

 H3TrkEffi->Write("");
 H3TrkEffi_AxisRawTrk->Write("");
 // H3TrkEffivsR->Write("");

 H0TrkEffi_AxisRawTrk_EtaCut->Write("");
 H1TrkEffi_AxisRawTrk_EtaCut->Write("");
 H2TrkEffi_AxisRawTrk_EtaCut->Write("");
 H3TrkEffi_AxisRawTrk_EtaCut->Write("");

 H0TrkEffi_AxisAvgMult_EtaCut->Write("");
 H1TrkEffi_AxisAvgMult_EtaCut->Write("");
 H2TrkEffi_AxisAvgMult_EtaCut->Write("");
 H3TrkEffi_AxisAvgMult_EtaCut->Write(""); 


 H0TrkEffi_AxisAvgMult->Write("");
 H1TrkEffi_AxisAvgMult->Write("");
 H2TrkEffi_AxisAvgMult->Write("");
 H3TrkEffi_AxisAvgMult->Write("");


 H0TrkEffi_AxisAvgMult_ZImpactCut->Write("");
 H1TrkEffi_AxisAvgMult_ZImpactCut->Write("");
 H2TrkEffi_AxisAvgMult_ZImpactCut->Write("");
 H3TrkEffi_AxisAvgMult_ZImpactCut->Write("");
 /*
 HXTrkEffi_Yplus->Write(""); 
   std::cout<<"M"<<std::endl;
 HXTrkEffi_Yminus->Write(""); 
 HXTrkEffivsRplus->Write(""); 
 HXTrkEffivsRminus->Write("");
 HX_ClusterPhi->Write("");
 */

ZIMPdistrRefH0->Write(""); ZIMPdistrRefH1->Write("");
ZIMPdistrRefH2->Write(""); ZIMPdistrRefH3->Write("");



 f->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(T2HalfQuarterTrkEfficiency);
