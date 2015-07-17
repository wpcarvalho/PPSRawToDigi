#include "TotemT1T2Validation/GeneratorAnalysis/interface/T2GeneratorAnalyzer.h"

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

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

#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"

#include "CLHEP/Vector/LorentzVector.h"


T2GeneratorAnalyzer::T2GeneratorAnalyzer(const edm::ParameterSet& iConfig){
  Chicut = iConfig.getParameter<double>("Chicut");
  PhiChiProbCut = iConfig.getParameter<double>("PhiChiProbCut");
  RChiProbCut = iConfig.getParameter<double>("RChiProbCut");  
  energy=iConfig.getParameter<double>("FitdzEnergy");
  DZScale=iConfig.getParameter<double>("DZScale");
  tracketamin=iConfig.getParameter<double>("TrkEtamin");
  tracketamax=iConfig.getParameter<double>("TrkEtaMAX");
  singleparticle=iConfig.getParameter<bool>("singleparticle");

  outputFileName = iConfig.getUntrackedParameter<std::string>("OutputFile");
  HepMCProductLabel = iConfig.getParameter<std::string>("HepMCProductLabel");
  CluLabel = iConfig.getParameter<std::string>("CluLabel");
  HitLabel = iConfig.getParameter<std::string>("HitLabel");
  RoadLabel = iConfig.getParameter<std::string>("RoadLabel");
  TrackLabel= iConfig.getParameter<std::string>("TrackLabel");
  FullInfos= iConfig.getParameter<bool>("FullInfos");
}


T2GeneratorAnalyzer::~T2GeneratorAnalyzer()
{
}




double PartCharge(int pdgcode){
  double carica=0;
 switch(pdgcode){
 case(11):
   carica=-1; //electron
   break;
case(-11):
   carica=1; 
   break;

 case(13):
	carica=-1; //muon
	break;
case(-13):
	carica=+1; 
	break;
 case(15):
	carica=-1; //tau
	break;
case(-15):
	carica=+1; //tau+
	break;

 case(211):
	carica=1; //picharg+
	break;
 case(-211):
	carica=-1; //picharg+
	break;

 case(213):
	carica=1; //rhocharg+
	break;
 case(-213):
	carica=-1; //rhocharg-
	break;

 case(321):
	carica=1; //kcharg+
	break;

 case(-321):
	carica=-1; //kcharg-
	break;

 case(2212):
	carica=+1; //proton
	break;
 case(-2212):
	carica=-1; //antiproton
	break;

 case(1114):
	carica=-1; //delta-
	break;
 case(2214):
	carica=-1; //delta+
	break;
 case(2224):
	carica=-1; //delta++
	break;
 case(411):
	carica=1; //D+
	break;
case(-411):
	carica=-1; //D-
	break;
case(521):
	carica=1; //B+
	break;
case(-521):
	carica=-1; //B-
	break;
 }
 return carica;
}


//DZ curve obtained with 50 GeV pion
double sigmaZfunctE50(double eta){
  double theDZ=0.0;
  double nepero=10.0;

  if(fabs(eta)<5.2)
    theDZ=3000.0;

  if(fabs(eta)>6.6)
    theDZ=4000.0;

  if((fabs(eta)>5.2)&&(fabs(eta)<6.6))
    theDZ= 7.79425*pow(nepero,6.0)-  5.16635*pow(nepero,6.0)*eta  +  1.28288*pow(nepero,6.0)*eta*eta  -1.41458*pow(nepero,5.0)*eta*eta*eta+ 5.84678*pow(nepero,3.0)*eta*eta*eta*eta;

  if((fabs(eta)>5.5)&&(fabs(eta)<5.6))
    theDZ=2000.0;

  return theDZ;
}

//DZ curve obtained with 10-50  GeV pion
double sigmaZfunctE30(double eta){

  double theDZ=0.0;
  double nepero=10.0;

  
  if(fabs(eta)<5.37)
     theDZ=2200.0;
  
   if(fabs(eta)>6.4)
     theDZ=4000.0;

   if((fabs(eta)>5.37)&&(fabs(eta)<6.4)){
       theDZ=1.11966*pow(nepero,5.0)-3.91045*pow(nepero,4.0)*eta+3.47461*pow(nepero,3.0)*eta*eta;
     }
   
 
   if((fabs(eta)>5.45)&&(fabs(eta)<5.7))  
     theDZ=2400.0;

  return theDZ;
}

unsigned int RawtoSymb(uint32_t thedet)
{
  T2DetId converter;
  unsigned int pl=converter.plane(thedet);
  unsigned int pls=converter.planeSide(thedet);
  unsigned int ht=converter.halfTelescope(thedet);
  unsigned int arm=converter.arm(thedet);
  unsigned int symbolic=pl*2+pls+ht*10+20*arm;	  
 
  return symbolic;
}



bool IsinT2(double eta, double tracketamin, double tracketamax){
  bool flag;
  if ((fabs(eta)>tracketamin)&&(fabs(eta)<tracketamax))
    flag=true;
  else
    flag=false;
  return flag;
}

//
// member functions
//
// ------------ method called to for each event  ------------
void T2GeneratorAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace HepMC;

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  /* LOADING OF ALL THE RECORDS FROM THE EVENT */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */


  // edm::ESHandle <ParticleDataTable> pdt;
  // if(!(pdt.isValid())) 
  // iSetup.getData(pdt);
  //std::cout<<"CC"<<std::endl;
  
  Handle<T2StripClusterCollection> t2strclcoll;
  
  Handle<T2PadClusterCollection> t2padclcoll;
  
  /*::::::Take particle generated with Gun::::::*/
  Handle<HepMCProduct> EvtHandle ;
  iEvent.getByLabel(HepMCProductLabel, EvtHandle ) ;
  const GenEvent* evt = EvtHandle->GetEvent();

  /*::::::Take Geant local hit on T2::::::*/
  /*
  Handle<CrossingFrame<PSimHit> > cFrame;
  iEvent.getByLabel("mix", "g4SimHitsTotemHitsT2Gem", cFrame);
  // get hits from G4Sim
  const string nameOfHits("TotemHitsT2Gem");
  auto_ptr<MixCollection<PSimHit> >  T2SimHits( new MixCollection<PSimHit>( cFrame.product() ) );


  // map hits to T2 plane
  map<int, PSimHitContainer> hitMap;
  for(MixCollection<PSimHit>::MixItr hitItr = T2SimHits->begin(); hitItr != T2SimHits->end(); ++hitItr) {
    hitMap[hitItr->detUnitId()].push_back(*hitItr);
  }
  */

  /*
  // :::::::::::::Take The Clusters::::::::::::

  Handle<T2StripClusterCollection> t2strclcoll;
  iEvent.getByLabel(CluLabel,"T2StripClusters",t2strclcoll);
  Handle<T2PadClusterCollection> t2padclcoll;
  iEvent.getByLabel(CluLabel,"T2PadClusters",t2padclcoll);

  //:::::::Take  T2  Hits::::::
  Handle<T2HitCollection> t2hitcoll;

  iEvent.getByLabel(HitLabel,"T2Hits",t2hitcoll);
  //::::::Take  T2  Roads:::::

  Handle<T2RoadCollection> t2roadcoll;
  iEvent.getByLabel(RoadLabel,"T2RoadColl",t2roadcoll);
  */

  //:::::: Take T2 tracks ::::::
   unsigned int numTrk=0;
  /*
  Handle<T1T2TrackCollection> trackCollection;
  
  if(FullInfos==true){
    iEvent.getByLabel(TrackLabel,"T2TrackColl",trackCollection);  
    iEvent.getByLabel(CluLabel,"T2PadClusters",t2padclcoll);
    iEvent.getByLabel(CluLabel,"T2StripClusters",t2strclcoll);
    numTrk=trackCollection->size();
  }
  */
   vectorMBALLMCGenerator.clear();vectorMBALLMCGeneratorPtCut.clear();
   vectorMBALLMCGenerator_K0s.clear();
   vectorMBALLMCGenerator_gammaEInLhcf.clear();
   for(unsigned int iu=0;iu<totnumberbinfordndeta+1;iu++)
    {
      vectorMBALLMCGenerator.push_back(0);vectorMBALLMCGeneratorPtCut.push_back(0);
      vectorMBALLMCGenerator_K0s.push_back(0);
    }
   
   for(unsigned int i=0;i<=35;i++)
     {
       vectorMBALLMCGenerator_gammaEInLhcf.push_back(0);
     }
  int chpAlice=0;
  int chp12pos=0;
  int chp12neg=0;int nump04575=0;
  bool  lhcbaccept=false; bool  t2accept=false; bool  aliceaccept=false;
  bool  CMSaccept=false; 
  //std::cout<<"A"<<std::endl;

  for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){

    double etag=(*p)->momentum().eta();


//if((*p)->pdg_id()==22) 
// if((fabs(etag)<8.99)&&(fabs(etag)>8.5))
     //std::cout<<"GAMMA FOUNDAA!"<<etag<<std::endl;

    //std::cout<<"Particle PID: "<<(*p)->pdg_id()<<"  charge:"<<charge<<" status:"<<(*p)->status()<<std::endl;

     if((*p)->status()==1)
       EtaParticleVsNumTrk->Fill(((*p)->momentum().eta()),numTrk);
      
   

    if(((*p)->pdg_id()==111)&&((*p)->momentum().eta()<7.5)&&((*p)->momentum().eta()>4.5))
      {
	nump04575++;
	if(((*p)->momentum().e())>1.)
	  P0EtaEnergy->Fill(((*p)->momentum().eta()),((*p)->momentum().e()));
      }

    
    //const ParticleData * part = pdt->particle((*p)->pdg_id());
    //int charge = part->charge();
   
    int charge=PartCharge((*p)->pdg_id());  //OBSOLETE WAY
 
    if((*p)->pdg_id()==(-211))
      if((*p)->status()==1){
	
	PmenoEtaEnergy2->Fill(((*p)->momentum().eta()),((*p)->momentum().e()));
	//  std::cout<<etag<<std::endl;
	if((fabs(etag)<6.5)&&(fabs(etag)>5.3)){
	  EnergyPiMenoPlusInT2->Fill((*p)->momentum().e());
	  
	}
      }

    if((*p)->pdg_id()==111){
      EnergyPi0InT2->Fill((*p)->momentum().e());
      
      //  std::cout<<etag<<std::endl;
      if((fabs(etag)<6.5)&&(fabs(etag)>5.3)){
	
	StatusPi0InT2->Fill((*p)->status());
      }
    }
    
    if((*p)->pdg_id()==22)
      if((fabs(etag)<8.99)&&(fabs(etag)>8.81)){
	//std::cout<<"GAMMA FOUND!"<<std::endl;
	int binE=((*p)->momentum().e())/100;
	if(binE<35)
	  vectorMBALLMCGenerator_gammaEInLhcf.at(binE)=vectorMBALLMCGenerator_gammaEInLhcf.at(binE)+1;
	else
	  std::cout<<"Warn Egamma="<<((*p)->momentum().e())<<std::endl;

      }


    if((*p)->pdg_id()==310){	
      
      int trkbinfordndeta_=(int)((etag-(-maxetafordndetahisto))/etabinsize); 
      if((trkbinfordndeta_<(int)totnumberbinfordndeta)&&(trkbinfordndeta_>0))
	vectorMBALLMCGenerator_K0s.at(trkbinfordndeta_)= vectorMBALLMCGenerator_K0s.at(trkbinfordndeta_)+1;      
    }

    if((abs(charge)>0)&&((*p)->status()==1))
      {	      	    	    
	
	if((etag<2.4)&&(etag>-2.4))
	  CMSaccept=true;

	if((etag<4.5)&&(etag>2))
	  lhcbaccept=true;

	if((fabs(etag)<6.5)&&(fabs(etag)>5.3))
	  t2accept=true;

	if(fabs(etag)<1){
	  aliceaccept=true;
	  chpAlice++;}

	Eta2_vs_Energy->Fill(etag,((*p)->momentum().e()));


	if(fabs(etag)<12.)
	  {
	    if((fabs(etag)<6.5)&&(fabs(etag)>5.3))
	      {
		StableChPartilceEnergy->Fill((*p)->momentum().e());
		if(etag>0)
		  chp12pos++; 
		else
		  chp12neg++;
	      }
	   int trkbinfordndeta_=(int)((etag-(-maxetafordndetahisto))/etabinsize); 
	   if((trkbinfordndeta_<(int)totnumberbinfordndeta)&&(trkbinfordndeta_>0)){
	     vectorMBALLMCGenerator.at(trkbinfordndeta_)= vectorMBALLMCGenerator.at(trkbinfordndeta_)+1;
	     if((*p)->momentum().perp()>0.03)
	        vectorMBALLMCGeneratorPtCut.at(trkbinfordndeta_)= vectorMBALLMCGeneratorPtCut.at(trkbinfordndeta_)+1;
	   }
	    //vectorMBALLMCGenerator_RecoAccepted.at(trkbinfordndeta_)= vectorMBALLMCGenerator_RecoAccepted.at(trkbinfordndeta_)+1;
	  }
      }
  }          

  // std::cout<<"B"<<std::endl;

  double etacentre=0.;
  for(unsigned int bin=0; bin<vectorMBALLMCGenerator_gammaEInLhcf.size(); bin++){
    double ecentre=bin*100+50;
    DNDE_gamma_NoEvtRequirement->Fill(ecentre,vectorMBALLMCGenerator_gammaEInLhcf.at(bin)*(1.0/100.));
  } 
  // std::cout<<"C"<<std::endl;

  for(unsigned int bin=0; bin<vectorMBALLMCGenerator.size(); bin++)
    {
      etacentre=-maxetafordndetahisto+bin*etabinsize+ etabinsize/2.0;    

      
      DNDetaMBALLMCGenerator_NoEvtRequirement->Fill(etacentre,vectorMBALLMCGenerator.at(bin)*(1.0/etabinsize));
      DNDetaMBALLMCGenerator_K0s_NoEvtRequirement->Fill(etacentre,vectorMBALLMCGenerator_K0s.at(bin)*(1.0/etabinsize));
      if(lhcbaccept)
	DNDetaMBALLMCGenerator_LHCbRequirement->Fill(etacentre,vectorMBALLMCGenerator.at(bin)*(1.0/etabinsize));

      if(t2accept){
	DNDetaMBALLMCGenerator_T2Requirement->Fill(etacentre,vectorMBALLMCGenerator.at(bin)*(1.0/etabinsize));
	DNDetaMBALLMCGenerator_GeneratorTriggered->Fill(etacentre,vectorMBALLMCGeneratorPtCut.at(bin)*(1.0/etabinsize));
      }

      if(aliceaccept)
	DNDetaMBALLMCGenerator_ALICE->Fill(etacentre,vectorMBALLMCGenerator.at(bin)*(1.0/etabinsize));

      if(CMSaccept){
	DNDetaMBALLMCGenerator_CMS->Fill(etacentre,vectorMBALLMCGenerator.at(bin)*(1.0/etabinsize));
	DNDetaMBALLMCGenerator_K0s_CMS->Fill(etacentre,vectorMBALLMCGenerator_K0s.at(bin)*(1.0/etabinsize));
      }
    }


  //std::cout<<"C"<<std::endl;



  




  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //



  /*






 //////////////////////Begin test///////////////



//Simulate a Trigger Condition
  if(numTrk>0){
    




  std::vector<unsigned int> idexplored;
    double rr=0.;double ff=0.;double xx=0.;double yy=0.;

 for(T2StripClusterCollection::const_iterator itstrip = t2strclcoll->begin(); itstrip != t2strclcoll->end(); itstrip++){
    

   vector<T2Cluster> stripClv = itstrip->second;
    T2DetId *detID =new T2DetId(itstrip->first);
    
    uint32_t cmsswdId= detID->calculateRawId(detID->arm(),detID->halfTelescope(),detID->plane(),detID->planeSide());
    unsigned int symbol=RawtoSymb(cmsswdId);

    unsigned int quarter=symbol/10;
    unsigned int plane=symbol%10;

   
    
    unsigned int countCutted0=0;
    unsigned int countCutted1=0;
    unsigned int countCutted2=0;
    unsigned int countCutted3=0;
    unsigned int countCutted0Plus=0;unsigned int countCutted0Minus=0;
    unsigned int countCutted1Plus=0;unsigned int countCutted1Minus=0;
    unsigned int countCutted2Plus=0;unsigned int countCutted2Minus=0;
    unsigned int countCutted3Plus=0;unsigned int countCutted3Minus=0;
    idexplored.push_back(symbol);
    
    for(unsigned int k=0;k<stripClv.size();k++){
      
      rr=stripClv[k].GetClusterR();
      ff=stripClv[k].GetClusterPhi()*3.14159/180.;    
      xx=rr*cos(ff);
      yy=rr*sin(ff);

      if(quarter==0)
	{
	  StripCluSizeVsPlaneAll3H0->Fill(plane,stripClv[k].GetNoOfEntries());
	  //CumulativeStripCluSizeAll3H0->Fill(stripClv[k].GetNoOfEntries());
	  if(stripClv[k].GetNoOfEntries()<=5){
	    countCutted0++;
	    if(yy>0)
	      countCutted0Plus++;
	    else
	      countCutted0Minus++;
	  }
	}
      
      if(quarter==1)
	{
	  StripCluSizeVsPlaneAll3H1->Fill(plane,stripClv[k].GetNoOfEntries());
	  if(stripClv[k].GetNoOfEntries()<=5){
	    countCutted1++;
	    if(yy>0)
	      countCutted1Plus++;
	    else
	      countCutted1Minus++;
	  }
	}

      if(quarter==2)
	{
	  StripCluSizeVsPlaneAll3H2->Fill(plane,stripClv[k].GetNoOfEntries());
	  if(stripClv[k].GetNoOfEntries()<=5){
	    countCutted2++;
	    if(yy>0)
	      countCutted2Plus++;
	    else
	      countCutted2Minus++;
	  }
	}

      if(quarter==3)
	{
	  StripCluSizeVsPlaneAll3H3->Fill(plane,stripClv[k].GetNoOfEntries());
	  if(stripClv[k].GetNoOfEntries()<=5){
	    countCutted3++;
	  
	    if(yy>0)
	      countCutted3Plus++;
	    else
	      countCutted3Minus++;
	  }
	}

      //if(symbol<10)
      //if(verbosity)
      //  std::cout<<" ev300Cluster Strip at "<<stripClv[k].GetClusterR()<<" plane: "<<symbol<<std::endl;
    }


    if(quarter==0){
      NumStripCluVsPlaneAll3H0->Fill(plane,stripClv.size());
      //NumStripCluVsPlaneAll3H0_Cutted->Fill(plane,countCutted0);
      //NumStripCluVsPlaneAll3H0_CuttedMinus->Fill(plane,countCutted0Minus);
      //NumStripCluVsPlaneAll3H0_CuttedPlus->Fill(plane,countCutted0Plus);	
      // CumulativeNumStripCluAll3H0->Fill(stripClv.size());
    }
    
    if(quarter==1){
      NumStripCluVsPlaneAll3H1->Fill(plane,stripClv.size());
      //NumStripCluVsPlaneAll3H1_Cutted->Fill(plane,countCutted1);
      //NumStripCluVsPlaneAll3H1_CuttedMinus->Fill(plane,countCutted1Minus);
     // NumStripCluVsPlaneAll3H1_CuttedPlus->Fill(plane,countCutted1Plus);
    }

    if(quarter==2){
      NumStripCluVsPlaneAll3H2->Fill(plane,stripClv.size());
      //NumStripCluVsPlaneAll3H2_Cutted->Fill(plane,countCutted2);
      //NumStripCluVsPlaneAll3H2_CuttedMinus->Fill(plane,countCutted2Minus);
     // NumStripCluVsPlaneAll3H2_CuttedPlus->Fill(plane,countCutted2Plus);
    }

    if(quarter==3){
      NumStripCluVsPlaneAll3H3->Fill(plane,stripClv.size());
      //NumStripCluVsPlaneAll3H3_Cutted->Fill(plane,countCutted3);
      //NumStripCluVsPlaneAll3H3_CuttedMinus->Fill(plane,countCutted3Minus);
     // NumStripCluVsPlaneAll3H3_CuttedPlus->Fill(plane,countCutted3Plus);
    }


  }

 for(unsigned int i=0; i<40;i++)
   {
     if (std::find(idexplored.begin(),idexplored.end(),i)==idexplored.end())
       {
	 unsigned int qq=i/10;
	 unsigned int plane=i%10;
	 if(qq==0){
	   NumStripCluVsPlaneAll3H0->Fill(plane,0);
	   
	 }
	 if(qq==1){
	   NumStripCluVsPlaneAll3H1->Fill(plane,0);
	 }
	 if(qq==2){
	   NumStripCluVsPlaneAll3H2->Fill(plane,0);
	 }
	 if(qq==3){
	   NumStripCluVsPlaneAll3H3->Fill(plane,0);              
	 }
       }
   }
 
 idexplored.clear();

 for(T2PadClusterCollection::const_iterator itpad= t2padclcoll->begin(); itpad != t2padclcoll->end(); itpad++){
    vector<T2Cluster> padClv = itpad->second;
    T2DetId *detID =new T2DetId(itpad->first);
    
    uint32_t cmsswdId= detID->calculateRawId(detID->arm(),detID->halfTelescope(),detID->plane(),detID->planeSide());
    unsigned int symbol=RawtoSymb(cmsswdId);
    unsigned int quarter=symbol/10;
    unsigned int plane=symbol%10;

    unsigned int countCutted0=0;
    unsigned int countCutted1=0;
    unsigned int countCutted2=0;
    unsigned int countCutted3=0;
  unsigned int countCutted0Plus=0;unsigned int countCutted0Minus=0;
    unsigned int countCutted1Plus=0;unsigned int countCutted1Minus=0;
    unsigned int countCutted2Plus=0;unsigned int countCutted2Minus=0;
    unsigned int countCutted3Plus=0;unsigned int countCutted3Minus=0;
    idexplored.push_back(symbol);

      for(unsigned int k=0;k<padClv.size();k++){
	rr=padClv[k].GetClusterR();
	ff=padClv[k].GetClusterPhi()*3.14159/180.;    
	xx=rr*cos(ff);
	yy=rr*sin(ff);

	
	//	PadClusterR_AllvsPlane[symbol]->Fill(padClv[k].GetClusterR());
	//	PadClusterSize_AllvsPlane[symbol]->Fill(padClv[k].GetNoOfEntries());

	if(quarter==0)
	  {
	    PadCluSizeVsPlaneAll3H0->Fill(plane,padClv[k].GetNoOfEntries());
	    //CumulativePadCluSizeAll3H0->Fill(padClv[k].GetNoOfEntries());
	    if(padClv[k].GetNoOfEntries()<=5){
	      countCutted0++;
	      if(yy>0)
		countCutted0Plus++;
	      else
		countCutted0Minus++;
	    }
	  }
      
	if(quarter==1)
	  {
	    PadCluSizeVsPlaneAll3H1->Fill(plane,padClv[k].GetNoOfEntries());
	    if(padClv[k].GetNoOfEntries()<=5){
	      countCutted1++;
	       if(yy>0)
		countCutted1Plus++;
	      else
		countCutted1Minus++;
	    }
	  }

	if(quarter==2)
	  {
	    PadCluSizeVsPlaneAll3H2->Fill(plane,padClv[k].GetNoOfEntries());
	    if(padClv[k].GetNoOfEntries()<=5){
	      countCutted2++;
	       if(yy>0)
		countCutted2Plus++;
	       else
		countCutted2Minus++;
	    }
	  }
	
	if(quarter==3)
	  {
	    PadCluSizeVsPlaneAll3H3->Fill(plane,padClv[k].GetNoOfEntries());
	    if(padClv[k].GetNoOfEntries()<=5){
	      countCutted3++;
	      if(yy>0)
		countCutted3Plus++;
	      else
		countCutted3Minus++;
	    }
	  }
		
    }


      if(quarter==0)
	{
	  NumPadCluVsPlaneAll3H0->Fill(plane,padClv.size());
	 // CumulativeNumPadCluAll3H0->Fill(padClv.size());
	 // NumPadCluVsPlaneAll3H0_Cutted->Fill(plane,countCutted0);
	 // NumPadCluVsPlaneAll3H0_CuttedMinus->Fill(plane,countCutted0Minus);
	 // NumPadCluVsPlaneAll3H0_CuttedPlus->Fill(plane,countCutted0Plus);
	  
	 // if(countCutted0<5)
	  //  NumPadCluVsPlaneAll3H0_Cutted_LOWMultipl->Fill(plane,countCutted0);
	}
      
      if(quarter==1)
	{
	  NumPadCluVsPlaneAll3H1->Fill(plane,padClv.size());
	  //NumPadCluVsPlaneAll3H1_Cutted->Fill(plane,countCutted1);
	  //NumPadCluVsPlaneAll3H1_CuttedMinus->Fill(plane,countCutted1Minus);
	  //NumPadCluVsPlaneAll3H1_CuttedPlus->Fill(plane,countCutted1Plus);

	  //if(countCutted1<5)
	   // NumPadCluVsPlaneAll3H1_Cutted_LOWMultipl->Fill(plane,countCutted1);
	}

      if(quarter==2)
	{
	  NumPadCluVsPlaneAll3H2->Fill(plane,padClv.size());
	  //NumPadCluVsPlaneAll3H2_Cutted->Fill(plane,countCutted2);
	  //NumPadCluVsPlaneAll3H2_CuttedMinus->Fill(plane,countCutted2Minus);
	  //NumPadCluVsPlaneAll3H2_CuttedPlus->Fill(plane,countCutted2Plus);

	  //if(countCutted2<5)
	   // NumPadCluVsPlaneAll3H2_Cutted_LOWMultipl->Fill(plane,countCutted2);
	}

      if(quarter==3)
	{
	   NumPadCluVsPlaneAll3H3->Fill(plane,padClv.size());   
	   //NumPadCluVsPlaneAll3H3_Cutted->Fill(plane,countCutted3);
	   //NumPadCluVsPlaneAll3H3_CuttedMinus->Fill(plane,countCutted3Minus);
	   //NumPadCluVsPlaneAll3H3_CuttedPlus->Fill(plane,countCutted3Plus);
	   //if(countCutted3<5)
	   //NumPadCluVsPlaneAll3H3_Cutted_LOWMultipl->Fill(plane,countCutted3);
	}

  }

 for(unsigned int i=0; i<40;i++)
   {
     if (std::find(idexplored.begin(),idexplored.end(),i)==idexplored.end())
       {
	 unsigned int qq=i/10;
	 unsigned int plane=i%10;
	 if(qq==0)
	   NumPadCluVsPlaneAll3H0->Fill(plane,0);
	 if(qq==1)
	   NumPadCluVsPlaneAll3H1->Fill(plane,0);
	 if(qq==2)
	   NumPadCluVsPlaneAll3H2->Fill(plane,0);
	 if(qq==3)
	   NumPadCluVsPlaneAll3H3->Fill(plane,0);              
       }
   }


  }//If NumTrk>0


*/


 //////////////////////end test///////////////



































  /*

    //--***************************************************         Reco TRACK        ************************************************************--/
  //--***************************************************        Reco TRACK        ************************************************************--/
  //--***************************************************        Reco TRACK        ************************************************************--/
  //std::cout<<"Track reco Only"<<std::endl;




    
  T1T2TrackCollection::const_iterator TrkCit;
  double trketa=0.;
  double trkphi=0.;
  unsigned int trackcounter=0;
  unsigned int trackcountergood=0;
  double DZgood;
  double numrecotrackcutmag0=0.;
  double chiRProb=0.;
  double chiPhiProb=0.;
  bool chi2condition;
  for(TrkCit=trackCollection->begin(); TrkCit!=trackCollection->end(); TrkCit++){
    trketa= (*TrkCit).Eta();
    trkphi=(*TrkCit).Phi()*180/3.14159265;
    Trketa->Fill(trketa);
    Trkphi->Fill(trkphi);
    chiRProb=TMath::Prob((*TrkCit).ChiSquaredR(),((*TrkCit).GetHitEntries()-2));
    chiPhiProb=TMath::Prob((*TrkCit).ChiSquaredPhi(),((*TrkCit).GetHitEntries()-1));  
    Chi2RProb->Fill(chiRProb);
    Chi2PhiProb->Fill(chiPhiProb);
    
    if((chiPhiProb<PhiChiProbCut)&&(chiRProb<RChiProbCut))  
      chi2condition=false;
    else
      chi2condition=true;

    DZgood=1000.;
    if (energy==30.)
      DZgood=sigmaZfunctE30(fabs(trketa));
    if (energy==50.)
      DZgood=sigmaZfunctE50(fabs(trketa));

    if((chi2condition)&&(IsinT2(trketa,tracketamin,tracketamax))&&(fabs((*TrkCit).Z_at_Rmin())<DZgood)){
      trackcountergood++;
      Trketagood->Fill(trketa);
      Trkphigood->Fill(trkphi);
      if(trketa>0)	       
	numrecotrackcutmag0= numrecotrackcutmag0+1;
      
      
    }
    trackcounter++;
  }

  double lastParticleeta=0.;
  double chp12mag0=0.;
  // double chp12min0=0.;
  for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
    lastParticleeta=(*p)->momentum().eta();
    if ((PartCharge((*p)->pdg_id())!=0)&&((*p)->status()==1)&&((*p)->momentum().eta()<6.5)&&((*p)->momentum().eta()>5.3))
      chp12mag0++;
    //  if ((PartCharge((*p)->pdg_id())!=0)&&((*p)->status()==1)&&((*p)->momentum().eta()<-5.3)&&((*p)->momentum().eta()>-6.5))
    //chp12min0++;

  }

  SingleParticleEfficiency->Fill(lastParticleeta,trackcounter);
  SingleParticleEfficiencyCut->Fill(lastParticleeta,trackcountergood);


  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  GEANT TRACK @@@@@@@@@@@@@@@@@@@@@@@@@@
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  GEANT TRACK @@@@@@@@@@@@@@@@@@@@@@@@@@
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  GEANT TRACK @@@@@@@@@@@@@@@@@@@@@@@@@@
 
  double SingleTrackEvent_EtaRec=0;
  double PI = 3.1415927;
  double SingleTrackEvent_PhiRec=0;
  

  //get G4 tracks
  
//std::cout<<"Track Reco + Geant"<<std::endl;
  Handle<SimTrackContainer> G4TrkContainer;
  iEvent.getByType(G4TrkContainer);
  if (!G4TrkContainer.isValid()) {
    LogError("TrackerHitAnalyzer::analyze") << "Unable to find SimTrack in event!";
    return;
  }
  
  
  int numsimtrk=fabs(G4TrkContainer->begin(); itTrk != G4TrkContainer->end());
  G4TrkContainerSize

 if((trackcountergood==1)&&(singleparticle)) 
  {
    vector< pair<double,double> > PrimarySimTracks;
    
    for (SimTrackContainer::const_iterator itTrk = G4TrkContainer->begin(); itTrk != G4TrkContainer->end(); ++itTrk) {
      double eta =0, phi =0, p =0;
      const HepLorentzVector G4Trk(itTrk->momentum().x(),itTrk->momentum().y(),itTrk->momentum().z(), itTrk->momentum().e() ) ;
      p = sqrt(G4Trk[0]*G4Trk[0]+G4Trk[1]*G4Trk[1]+G4Trk[2]*G4Trk[2]);
      if ( p == 0){
	LogError("TrackerHitAnalyzer::analyze") << "TrackerTest::INFO: Primary has p = 0 ";
      } else {
	double costheta  = G4Trk[2]/p;
	double theta = acos(TMath::Min(TMath::Max(costheta, -1.),1.));
	eta = -log(tan(theta/2));
	if ( G4Trk[0] != 0 || G4Trk[1] != 0)
	  phi = atan2(G4Trk[1],G4Trk[0]);

	if(phi<0)
	  phi = 2*PI + phi;
	phi=phi*180/PI;

	pair<double,double> *myprimarytrack = new pair<double,double>(eta,phi);
	PrimarySimTracks.push_back(*myprimarytrack);
      }
    }

    // comparison reconstructed track  - GEANT track    
    for(T1T2TrackCollection::const_iterator TC_it=trackCollection->begin(); TC_it!=trackCollection->end(); TC_it++)
      {
	DZgood=1000.;
	if (energy==30.){
	  DZgood=sigmaZfunctE30(fabs((*TC_it).Eta()));
	}
	if (energy==50.) {
	  DZgood=sigmaZfunctE50(fabs((*TC_it).Eta()));
	}

	SingleTrackEvent_EtaRec = (*TC_it).Eta();
	SingleTrackEvent_PhiRec = (*TC_it).Phi()*180.0/3.14159265;

	for(vector< pair<double,double> >::iterator pair_it=PrimarySimTracks.begin(); pair_it!=PrimarySimTracks.end(); pair_it++){
	  double DE_temp = SingleTrackEvent_EtaRec - (*pair_it).first;
	  double DF_temp = SingleTrackEvent_PhiRec - (*pair_it).second;

	  ////std::cout << "Track DE, DF = " << DE_temp<< " , " << DF_temp << std::endl;
	  ////std::cout << "Phi reco: "<< SingleTrackEvent_PhiRec << "Phi G4: "  << (*pair_it).second << std::endl;
	  chiRProb=TMath::Prob((*TC_it).ChiSquaredR(),((*TC_it).GetHitEntries()-2));
	  chiPhiProb=TMath::Prob((*TC_it).ChiSquaredPhi(),((*TC_it).GetHitEntries()-1));          
	  if((chiPhiProb<PhiChiProbCut)&&(chiRProb<RChiProbCut))  
	    chi2condition=false;
	  else
	    chi2condition=true;

	  if((chi2condition)&&(IsinT2((*TC_it).Eta(),tracketamin,tracketamax))){
	    double DZgood=1000.;

	    if (energy==30.)
	      DZgood=sigmaZfunctE30(fabs((*TC_it).Eta()));
	    if (energy==50.)
	      DZgood=sigmaZfunctE50(fabs((*TC_it).Eta()));
	    
	    if(fabs((*TC_it).Z_at_Rmin())<DZgood){
	      DEtaGoodTrk->Fill(DE_temp);
	      DPhiGoodTrk->Fill(DF_temp);
	    }

	  }
	}
      }
  }

  //::::::::::Local variable Declaration:::::::::
  double lmyx=0.0;
  double lmyy=0.0;
  double lmyr=0.0;
  double phipart=0.0;

  T2DetId myT2Det;
  int myhalftele;
  int myplaneside;
  int myplane;
  double zdetshift;
  double zglobhit;

  //int zmm[2][10]; //2 row for halftelescope. 10 coloumn for planes
  double z1=13828.3;          //14035.605; //first Gem first drift gas zone (mm)
  double planedist= 86.0;
  double btbdist=24.6;   //25.0;
  double ovdist=43.0;
  double zinsidedet=1.5;           //4.5; From first Drift zone to RO board = 9mm


  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
//std::cout<<"strip Only"<<std::endl;

  for(T2StripClusterCollection::const_iterator itstrip = t2strclcoll->begin(); itstrip != t2strclcoll->end(); itstrip++){
    vector<T2Cluster> stripClv = itstrip->second;
    for(unsigned int k=0;k<stripClv.size();k++){
      clusterstripentries->Fill(stripClv[k].GetNoOfEntries());
    }
  }


//std::cout<<"strip+geant"<<std::endl;
 if((trackcountergood==1)&&(singleparticle))  
 for(T2StripClusterCollection::const_iterator itstrip = t2strclcoll->begin(); itstrip != t2strclcoll->end(); itstrip++){
      myT2Det = itstrip->first;
      myhalftele= myT2Det.halfTelescope();
      myplaneside= myT2Det.planeSide();
      myplane= myT2Det.plane();
      vector<T2Cluster> stripClv = itstrip->second;
      for(unsigned int k=0;k<stripClv.size();k++){
    

      for(map<int, PSimHitContainer>::const_iterator hitMapItr = hitMap.begin(); hitMapItr != hitMap.end(); ++hitMapItr){

          const PSimHitContainer & planeSimHits = hitMapItr->second;
          T2DetId *theT2DetId =new T2DetId(hitMapItr->first);

          if((theT2DetId->plane()==myplane)&&(theT2DetId->planeSide()==myplaneside)&&(planeSimHits.size()==1)&&(stripClv.size()==1)){
            for (unsigned int l=0; l<planeSimHits.size(); l++){

              PSimHit myplanehits= planeSimHits[l];
              lmyx= myplanehits.localPosition().x();
              lmyy= myplanehits.localPosition().y();
              lmyr= sqrt(lmyx*lmyx + lmyy*lmyy);
	      diffRCluHit->Fill(stripClv[k].GetClusterR()-lmyr);
	      
	    }
	  }
      }

    }
  }



  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
//std::cout<<"pad+geant"<<std::endl;

for(T2PadClusterCollection::const_iterator itpad = t2padclcoll->begin(); itpad != t2padclcoll->end(); itpad++){
  vector<T2Cluster> padClv = itpad->second;
  for(unsigned int k=0;k<padClv.size();k++){
    clusterpadentries->Fill(padClv[k].GetNoOfEntries());
  }
}

  if((trackcountergood==1)&&(singleparticle))                       //Comupte the differences only if there is no secondaries
    for(T2PadClusterCollection::const_iterator itpad = t2padclcoll->begin(); itpad != t2padclcoll->end(); itpad++){

      myT2Det = itpad->first;
      myhalftele= myT2Det.halfTelescope();
      myplaneside= myT2Det.planeSide();
      myplane= myT2Det.plane();

      vector<T2Cluster> padClv = itpad->second;

      for(unsigned int k=0;k<padClv.size();k++){

	// clusterpadentries->Fill(padClv[k].GetNoOfEntries());

        // now we run for each plane the digi-simulation
        for(map<int, PSimHitContainer>::const_iterator hitMapItr = hitMap.begin(); hitMapItr != hitMap.end(); ++hitMapItr){

          const PSimHitContainer & planeSimHits = hitMapItr->second;
          T2DetId *theT2DetId =new T2DetId(hitMapItr->first);

          if((theT2DetId->plane()==myplane)&&(theT2DetId->planeSide()==myplaneside)&&(planeSimHits.size()==1)&&(padClv.size()==1)){
            for (unsigned int l=0; l<planeSimHits.size(); l++){

              PSimHit myplanehits= planeSimHits[l];
              lmyx= myplanehits.localPosition().x();
              lmyy= myplanehits.localPosition().y();
              lmyr= sqrt(lmyx*lmyx + lmyy*lmyy);

              //strange GEANT x-y sign convention
              if((myplaneside==1)&&(myhalftele==0)){
                phipart=atan2(lmyy,lmyx);
                phipart=phipart*180.0/3.14159265;
                if (phipart<0)
                  phipart= 360.0 - fabs(phipart);
              }
              if((myplaneside==0)&&(myhalftele==0)){
                phipart=atan2(-lmyy,lmyx);
                phipart=phipart*180.0/3.14159265;
                if (phipart<0)
                  phipart= 360.0 - fabs(phipart);
              }
              if((myplaneside==1)&&(myhalftele==1)){
                phipart=atan2(-lmyy,-lmyx);
                phipart=phipart*180.0/3.14159265;
                if (phipart<0)
                  phipart= 360.0 - fabs(phipart);
              }
              if((myplaneside==0)&&(myhalftele==1)){
                phipart=atan2(lmyy,-lmyx);
                phipart=phipart*180.0/3.14159265;
                if (phipart<0)
                  phipart= 360.0 - fabs(phipart);
                ////std::cout<<"PS0 HT1  Reco: "<<padClv[k].GetClusterPhi()<<"; Geant: "<<phipart<<";    X-Y: "<<lmyx<<"-"<<lmyy<<std::endl;
              }
              diffphiCluHit->Fill(padClv[k].GetClusterPhi()-phipart);
            }
          }
        }                 // RPhiHitGeant

        for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
          phipart=(*p)->momentum().phi()*180.0/3.14159265;

          if (phipart<0)
            phipart= 360.0 - fabs(phipart);  // From particleGun Phi to Clusterizator Phi
          diffphiCluGun->Fill(padClv[k].GetClusterPhi()-phipart);
        }
      }
    }


  // **************************************************** HIT STUDIES *********************************************** //
  // **************************************************** HIT STUDIES ************************************************ //
  // **************************************************** HIT STUDIES ************************************************ //
  // **************************************************** HIT STUDIES ************************************************ //
  // **************************************************** HIT STUDIES ************************************************ //

  // COMPARISON RECOHIT LOCALHIT
  // now we run for each plane the digi-simulation
//std::cout<<"hit+geant"<<std::endl;
 if((trackcountergood==1)&&(singleparticle))
    for(map<int, PSimHitContainer>::const_iterator hitMapItr = hitMap.begin();  hitMapItr != hitMap.end(); ++hitMapItr){

      const PSimHitContainer & planeSimHits = hitMapItr->second;
      T2DetId *theT2DetId =new T2DetId(hitMapItr->first);

      myhalftele= theT2DetId->halfTelescope();
      myplaneside= theT2DetId->planeSide();
      myplane= theT2DetId->plane();
      zdetshift=myplane*planedist+myplaneside*btbdist;
      if (myhalftele==0)
        zdetshift=zdetshift+ovdist;
      zglobhit=zdetshift+z1+zinsidedet;

      for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++){
        if (fabs(zglobhit-(ithit->GetHitZ()))<0.5)
          for (unsigned int l=0; l<planeSimHits.size(); l++){
            PSimHit myplanehits= planeSimHits[l];
            lmyx= myplanehits.localPosition().x();
            lmyy= myplanehits.localPosition().y();  //prima entryPoint()
            lmyr= sqrt(lmyx*lmyx + lmyy*lmyy);

            RLocHRecoH->Fill((ithit->GetHitR())-lmyr);
          }
      }
    }

  //Studies on Reco Hit

  double hitr=0.0;
  double hitphi=0.0;
  double hitdr=0.0;
  double hitdphi=0.0;
  unsigned int hnumstrip=0;
  unsigned int hnumpad=0;
  double phipar;
//std::cout<<"hit onlyt"<<std::endl;
  if(trackcountergood==1)
    for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++) {
      hitr=ithit->GetHitR();
      hitphi=ithit->GetHitPhi();
      hitdr=ithit->GetHitDR();
      hitdphi=ithit->GetHitDPhi();
      hnumstrip=ithit->GetHitNumStrip();
      hnumpad=ithit->GetHitNumPad();

      for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
        phipar=(*p)->momentum().phi()*180.0/3.14159265;

        if (phipar<0)
          phipar= 360.0 + phipar;  // From particleGun Phi to Clusterizator Phi

        if ((fabs((*p)->momentum().eta())<6.5)&&(fabs((*p)->momentum().eta())>5.3)){
          if (fabs(phipar-hitphi)<25.0){
            //cout<<"fabs(phipar-hitphi)  ="<<fabs(phipar-hitphi)<<"   particle gen phi: "<<phipar<<"   Hit phi"<< hitphi<<endl;
            diffphiGUNHIT->Fill(phipar-hitphi);
          }
        }
      }
    }



if(chp12mag0>0)
  {
 MultiParticleEfficiencyCutNorm->Fill(chp12mag0,((double)numrecotrackcutmag0)/((double)chp12mag0));
 //cout<<" eta>0 Ch -- Tr:   "<<chp12mag0 <<" --  "<<numrecotrackcutmag0<<endl;
 MultiParticleEfficiencyCut->Fill(((double)chp12mag0),((double)numrecotrackcutmag0));
}

  */

  numevent++;


  if((numevent%2000)==0)
    std::cout<<"Evt: "<<numevent<<std::endl;
}







// ------------ method called once each job just before starting event loop  ------------
void T2GeneratorAnalyzer::beginJob()
{
  TH1::AddDirectory(kFALSE);
  numevent=0;
std::cout<<"Abeg"<<std::endl;

  totnumberbinfordndeta=480;
  maxetafordndetahisto=12.0; //range assumed simmetric around eta=0
  etabinsize=(maxetafordndetahisto*2)/totnumberbinfordndeta;

  EtaParticleVsNumTrk= std::auto_ptr<TProfile> (new TProfile("EtaParticleVsNumTrk","#eta vs NumTrk",50,4.4,5.6));
  EtaParticleVsNumTrk->SetYTitle("NumTrk");

  DNDetaMBALLMCGenerator_LHCbRequirement = std::auto_ptr<TProfile> (new TProfile("DNDetaMBALLMCGenerator_LHCbRequirement","MB Charged Particle dN/d#eta (at IP) (Pythia MC)",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

  DNDetaMBALLMCGenerator_T2Requirement = std::auto_ptr<TProfile> (new TProfile("DNDetaMBALLMCGenerator_T2Requirement","MB Charged Particle dN/d#eta (at IP) (Pythia MC)",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

  DNDetaMBALLMCGenerator_K0s_CMS= std::auto_ptr<TProfile> (new TProfile("DNDetaMBALLMCGenerator_K0s_CMS","MB K0s dN/d#eta (at IP) (Pythia MC)",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
DNDetaMBALLMCGenerator_K0s_NoEvtRequirement= std::auto_ptr<TProfile> (new TProfile("DNDetaMBALLMCGenerator_K0s_NoEvtRequirement","MB K0s dN/d#eta (at IP) (Pythia MC)",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 
 DNDetaMBALLMCGenerator_ALICE= std::auto_ptr<TProfile> (new TProfile("DNDetaMBALLMCGenerator_ALICE","MB Charged Particle dN/d#eta (at IP) (Pythia MC)",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

  DNDetaMBALLMCGenerator_CMS= std::auto_ptr<TProfile> (new TProfile("DNDetaMBALLMCGenerator_CMS","MB Charged Particle dN/d#eta (at IP) (Pythia MC)",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

 
  DNDetaMBALLMCGenerator_GeneratorTriggered= std::auto_ptr<TProfile> (new TProfile("DNDetaMBALLMCGenerator_GeneratorTriggered","DNDetaMBALLMCGenerator_GeneratorTriggered",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

  DNDetaMBALLMCGenerator_NoEvtRequirement = std::auto_ptr<TProfile> (new TProfile("DNDetaMBALLMCGenerator_NoEvtRequirement","MB Charged Particle dN/d#eta (at IP) (Pythia MC)",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
  DNDetaMBALLMCGenerator_NoEvtRequirement->SetYTitle("dN/d#eta");
  DNDetaMBALLMCGenerator_NoEvtRequirement->SetXTitle("#eta");
  
  DNDE_gamma_NoEvtRequirement= std::auto_ptr<TProfile> (new TProfile("DNDE_gamma_NoEvtRequirement","MB #gamma dN/dE (at IP) 8.81-8.99 360deg (Pythia MC)",35,0,3500));

  // Chi2PhiProbLogy= std::auto_ptr<TCanvas>(new TCanvas("Chi2PhiProblog","Azimuthal #chi^{2} probability",400,400));
  // Chi2PhiProbLogy->SetDirectory(0);
  // Chi2RProbLogy= std::auto_ptr<TCanvas>(new TCanvas("Chi2RProblog","Radial #chi^{2} probability",400,400));
  //Chi2RProbLogy->SetDirectory(0);

  clusterstripentries = std::auto_ptr<TH1F>(new TH1F("clusterstripentries","Number of strips in strip-clusters",10, -0.5, 9.5));
  clusterstripentries->SetDirectory(0); 

  clusterpadentries = std::auto_ptr<TH1F>(new TH1F("clusterpadentries","Number of pads in pad-clusters",10, -0.5, 9.5));
  clusterpadentries->SetDirectory(0);
 
  diffphiGUNHIT= std::auto_ptr<TH1F>(new TH1F("diffphiGUNHIT","#Phi Generated particle - #Phi Reconstructed Hit",120, -30, 30));
  diffphiGUNHIT->SetDirectory(0); 

  RLocHRecoH= std::auto_ptr<TH1F>(new  TH1F("RLocHRecoH"," Reconstructed Hit R - GEANT4 Hit R",100,-1.5, 1.5));
  RLocHRecoH->SetDirectory(0); 

  diffRCluHit= std::auto_ptr<TH1F>(new  TH1F("diffRCluHit"," Cluster R - GEANT4 Hit R",100,-1.5, 1.5));
  diffRCluHit->SetDirectory(0); 

  Trketa= std::auto_ptr<TH1F>(new TH1F("Trketa","Reconstructed Track #eta (no cut used)",300,3.0,8.0));
  Trketa->SetDirectory(0); 

  StatusPi0InT2= std::auto_ptr<TH1F>(new TH1F("StatusPi0InT2","StatusPi0InT2",5,-0.5,4.5)); 
  EnergyPi0InT2= std::auto_ptr<TH1F>(new TH1F("EnergyPi0InT2","EnergyPi0InT2",300,0.,600.)); 
  EnergyPiMenoPlusInT2= std::auto_ptr<TH1F>(new TH1F("EnergyPiMenoPlusInT2","EnergyPiMenoPlusInT2",300,0.,600.)); 

  Trkphi= std::auto_ptr<TH1F>(new TH1F("Trkphi","Reconstructed Track #phi (no cut used)",361,0,360));
  Trkphi->SetDirectory(0); 

  DPhiGoodTrk= std::auto_ptr<TH1F>(new TH1F("DPhiGoodTrk","#phi Track - #phi GEANT4 Track (Only primary tracks are considered)",51,-25.5,25.5));
  DPhiGoodTrk->SetDirectory(0); 

  DEtaGoodTrk= std::auto_ptr<TH1F>(new TH1F("DEtaGoodTrk","#eta Track - #phi GEANT4 Track (Only primary tracks are considered)",200,-1.0,1.0));
  DEtaGoodTrk->SetDirectory(0); 

  diffphiCluGun = std::auto_ptr<TH1F>(new TH1F("diffphiCluGun","Cluster #Phi - #Phi Generated particle",81, -40.5, 40.5));
  diffphiCluGun->SetDirectory(0);
 
  diffphiCluHit= std::auto_ptr<TH1F>(new TH1F("diffphiCluHit","Cluster #Phi - GEANT4 Hit #Phi",50, -5, 5));
  diffphiCluHit->SetDirectory(0); 

  Trketagood= std::auto_ptr<TH1F>(new TH1F("Trketagood","Reconstructed Track #eta (all cut used)",300,3.0,8.0));
  Trketagood->SetDirectory(0); 

  Trkphigood= std::auto_ptr<TH1F>(new TH1F("Trkphigood","Reconstructed Track #phi (all cut used)",361,0,360));
  Trkphigood->SetDirectory(0); 

  Chi2RProb = std::auto_ptr<TH1F>(new TH1F("Chi2RProb","Radial #chi^{2} probability",15000,0.0,1.05));
  Chi2RProb->SetDirectory(0); 

  Chi2PhiProb = std::auto_ptr<TH1F>(new TH1F("Chi2PhiProb","Azimuthal #chi^{2} probability",15000,0.0,1.05));
  Chi2PhiProb->SetDirectory(0);

  SingleParticleEfficiencyCut = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyCut","# Tracks reconstructed vs generated particle #eta (all cuts)",34,4.35,7.75,""));
  SingleParticleEfficiencyCut->SetDirectory(0);
 
  SingleParticleEfficiency = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiency","# Tracks reconstructed vs generated particle #eta (all track)",34,4.35,7.75,""));
  SingleParticleEfficiency->SetDirectory(0); 

  MultiParticleEfficiencyCut= std::auto_ptr<TProfile>(new TProfile("MultiParticleEfficiencyCut","<# Tracks> vs # Charged Particles (all track cuts)",12,-0.5,11.5,""));
  MultiParticleEfficiencyCut->SetDirectory(0); 

  MultiParticleEfficiencyCutNorm= std::auto_ptr<TProfile>(new TProfile("MultiParticleEfficiencyCutNorm","<# Tracks / # Charged Particles > vs # Charged Particles (all track cuts)",12,-0.5,11.5,"")); 
  MultiParticleEfficiencyCutNorm->SetDirectory(0); 
  
  P0EtaEnergy= std::auto_ptr<TH2F>(new TH2F("P0EtaEnergy","#eta vs Energy #pi^{0}",30,4.5,7.5,40,1.,400.));
  P0EtaEnergy->SetDirectory(0);

  PmenoEtaEnergy2= std::auto_ptr<TH2F>(new TH2F("PmenoEtaEnergy2","#eta vs Energy #pi^{-}",30,4.5,7.5,40,1.,400.));









NumPadCluVsPlaneAll3H0=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H0","Num Pad Cluster  vs Plane -cluinfo- H0",12,-1.5,10.5));
 PadCluSizeVsPlaneAll3H0=std::auto_ptr<TProfile> (new TProfile("PadCluSizeVsPlaneAll3H0","Pad Cluster Size vs Plane -cluinfo- H0",12,-1.5,10.5));
NumPadCluVsPlaneAll3H1=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H1","Num Pad Cluster  vs Plane -cluinfo- H1",12,-1.5,10.5));
 PadCluSizeVsPlaneAll3H1=std::auto_ptr<TProfile> (new TProfile("PadCluSizeVsPlaneAll3H1","Pad Cluster Size vs Plane -cluinfo- H1",12,-1.5,10.5));
NumPadCluVsPlaneAll3H2=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H2","Num Pad Cluster  vs Plane -cluinfo- H2",12,-1.5,10.5));
 PadCluSizeVsPlaneAll3H2=std::auto_ptr<TProfile> (new TProfile("PadCluSizeVsPlaneAll3H2","Pad Cluster Size vs Plane -cluinfo- H2",12,-1.5,10.5));
NumPadCluVsPlaneAll3H3=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H3","Num Pad Cluster  vs Plane -cluinfo- H3",12,-1.5,10.5));
 PadCluSizeVsPlaneAll3H3=std::auto_ptr<TProfile> (new TProfile("PadCluSizeVsPlaneAll3H3","Pad Cluster Size vs Plane -cluinfo- H3",12,-1.5,10.5));
NumStripCluVsPlaneAll3H0=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H0","Num Strip Cluster  vs Plane -cluinfo- H0",12,-1.5,10.5));
 StripCluSizeVsPlaneAll3H0=std::auto_ptr<TProfile> (new TProfile("StripCluSizeVsPlaneAll3H0","Strip Cluster Size vs Plane -cluinfo- H0",12,-1.5,10.5));
NumStripCluVsPlaneAll3H1=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H1","Num Strip Cluster  vs Plane -cluinfo- H1",12,-1.5,10.5));
 StripCluSizeVsPlaneAll3H1=std::auto_ptr<TProfile> (new TProfile("StripCluSizeVsPlaneAll3H1","Strip Cluster  vs Plane -cluinfo- H1",12,-1.5,10.5));
NumStripCluVsPlaneAll3H2=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H2","Num Strip Cluster Size vs Plane -cluinfo- H2",12,-1.5,10.5));
 StripCluSizeVsPlaneAll3H2=std::auto_ptr<TProfile> (new TProfile("StripCluSizeVsPlaneAll3H2","Strip Cluster  vs Plane -cluinfo- H2",12,-1.5,10.5));
NumStripCluVsPlaneAll3H3=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H3","Num Strip Cluster Size vs Plane -cluinfo- H3",12,-1.5,10.5));
 StripCluSizeVsPlaneAll3H3=std::auto_ptr<TProfile> (new TProfile("StripCluSizeVsPlaneAll3H3","Strip Cluster  vs Plane -cluinfo- H3",12,-1.5,10.5));
















  clusterstripentries->SetXTitle("# strips in one strip-cluster");
  clusterstripentries->SetYTitle("# cluster-strip");
  clusterpadentries->SetXTitle("# pad in one pad-cluster");
  clusterpadentries->SetYTitle("# cluster-pad");
  clusterpadentries->GetYaxis()->SetTitleOffset(1.25);
  clusterstripentries->GetYaxis()->SetTitleOffset(1.25);
  diffphiGUNHIT->SetXTitle("#Delta #Phi (deg)");
  RLocHRecoH->SetXTitle("#Delta R (mm)");
  diffRCluHit->SetXTitle("#Delta R (mm)");
  Trketa->SetXTitle("Track #eta");
  Trkphi->SetXTitle("Track #phi (deg)");
  DEtaGoodTrk->SetXTitle("#Delta #eta");
  diffphiCluGun->SetXTitle("#Delta #Phi (deg)");
  diffphiCluHit->SetXTitle("#Delta #Phi (deg)");
  DPhiGoodTrk->SetXTitle("#Delta #Phi (deg)");
  Trkphigood->SetXTitle("Track #phi (deg)");
  Trketagood->SetXTitle("Track #eta");
  SingleParticleEfficiency->SetXTitle("Particle #eta");
  SingleParticleEfficiencyCut->SetXTitle("Particle #eta");
  SingleParticleEfficiency->SetYTitle("Pseudo-Efficiency");
  SingleParticleEfficiencyCut->SetYTitle("Efficiency");
  MultiParticleEfficiencyCutNorm->SetYTitle("Efficiency");
  MultiParticleEfficiencyCutNorm->SetXTitle("ch. particle in T2");  
  MultiParticleEfficiencyCut->SetYTitle("<Track>");
  MultiParticleEfficiencyCut->SetYTitle("ch. particle in T2");    
  Chi2RProb->SetXTitle("Probability");
  Chi2PhiProb->SetXTitle("Probability");




  StableChPartilceEnergy= std::auto_ptr<TH1D>(new TH1D("StableChPartilceEnergy","Stable Ch Partilce Energy at the IP",60,0.,300.));
  StableChPartilceEnergy->SetYTitle("E (GeV)");
  StableChPartilceEnergy->SetDirectory(0);
  
  Eta2_vs_Energy= std::auto_ptr<TProfile>(new TProfile("Eta2_vs_Energy","Eta2_vs_Energy",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
  Eta2_vs_Energy->SetYTitle("E (GeV)"); 
  Eta2_vs_Energy->SetXTitle("eta (GeV)"); 
  Eta2_vs_Energy->SetDirectory(0); 
  //Chi2RProbLogy->cd();
  //Chi2RProb->Draw("");
  //Chi2RProbLogy->SetLogy();
  //Chi2PhiProbLogy->cd();
  //Chi2PhiProb->Draw("");
  //Chi2PhiProbLogy->SetLogy();
  //std::cout<<"Abeg-end"<<std::endl;

}

// ------------ method called once each job just after ending the event loop  ------------

void T2GeneratorAnalyzer::endJob()
{
  //std::cout << " Aend-beg"<<std::endl;
TFile *f = TFile::Open(outputFileName.c_str(), "recreate");
 if( !f || !f->IsWritable() ){
   std::cout << "Output file not opened correctly !!" << std::endl;
 }
  clusterstripentries->Write("");
  clusterpadentries->Write("");
  diffphiGUNHIT->Write("");
  DNDetaMBALLMCGenerator_NoEvtRequirement->Write("");
  DNDetaMBALLMCGenerator_T2Requirement->Write(""); 
  DNDetaMBALLMCGenerator_LHCbRequirement->Write("");
  DNDetaMBALLMCGenerator_K0s_CMS->Write("");
  DNDetaMBALLMCGenerator_K0s_NoEvtRequirement->Write("");
  DNDE_gamma_NoEvtRequirement->Write("");
  DNDetaMBALLMCGenerator_CMS->Write("");
  DNDetaMBALLMCGenerator_GeneratorTriggered->Write("");
  DNDetaMBALLMCGenerator_ALICE->Write("");
  DPhiGoodTrk->Write("");
  DEtaGoodTrk->Write("");
  RLocHRecoH->Write("");
  diffRCluHit->Write("");
  diffphiCluGun->Write("");
  diffphiCluHit->Write("");
  SingleParticleEfficiencyCut->Write("");
  SingleParticleEfficiency->Write("");

  EtaParticleVsNumTrk->Write("");

  StableChPartilceEnergy->Write(""); 
  StatusPi0InT2->Write("");  
  Eta2_vs_Energy->Write("");
  EnergyPi0InT2->Write("");  
  EnergyPiMenoPlusInT2->Write(""); 
  Trketagood->Write("");
  Trkphigood->Write("");
  Trketa->Write("");
  Trkphi->Write("");
  MultiParticleEfficiencyCutNorm->Write("");
  MultiParticleEfficiencyCut->Write("");
  Chi2PhiProb->Write("");
  Chi2RProb->Write(""); 
  //Chi2PhiProbLogy->Write("");
  //  Chi2RProbLogy->Write("");
  
  PmenoEtaEnergy2->Write(""); 
  P0EtaEnergy->Write(""); 






  NumPadCluVsPlaneAll3H0->Write("");
 PadCluSizeVsPlaneAll3H0->Write("");
NumPadCluVsPlaneAll3H1->Write("");
 PadCluSizeVsPlaneAll3H1->Write("");
NumPadCluVsPlaneAll3H2->Write("");
 PadCluSizeVsPlaneAll3H2->Write("");
NumPadCluVsPlaneAll3H3->Write("");
 PadCluSizeVsPlaneAll3H3->Write("");
NumStripCluVsPlaneAll3H0->Write("");
 StripCluSizeVsPlaneAll3H0->Write("");
NumStripCluVsPlaneAll3H1->Write("");
 StripCluSizeVsPlaneAll3H1->Write("");
NumStripCluVsPlaneAll3H2->Write("");
 StripCluSizeVsPlaneAll3H2->Write("");
NumStripCluVsPlaneAll3H3->Write("");
 StripCluSizeVsPlaneAll3H3->Write("");




 // std::cout << " Aend-end"<<std::endl;



  f->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(T2GeneratorAnalyzer);
