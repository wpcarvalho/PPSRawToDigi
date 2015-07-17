// -*- C++ -*-
//

// system include files

#include <cmath>
#include <iostream>
#include <list>
#include <map>

// user include files
//#include "TotemT1T2Validation/T2BackgroundAnalysis/interface/T2BackgroundAn.h"
//#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TotemT1T2Validation/T2BackgroundAnalysis/interface/T2BackgroundAn.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "TRandom.h"
//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

/////////////// WARNING //////////////////////////

//Look to the tag:
//Very-Important-Parameter in order to see some key 
//parameter of the effi truth analysis 

//////////////////////////////////////////////////

T2BackgroundAn::T2BackgroundAn(const edm::ParameterSet& iConfig)
{
 
  T2StripDigiCollectionLabel = iConfig.getParameter<edm::InputTag>("T2StripDigiCollectionLabel");
  T2PadDigiCollectionLabel = iConfig.getParameter<edm::InputTag>("T2PadDigiCollectionLabelLabel");
  SimVertexContainerLabel = iConfig.getParameter<edm::InputTag>("SimVertexContainerLabel");
  SimTrackContainerLabel = iConfig.getParameter<edm::InputTag>("SimTrackContainerLabel");



  //Read parameters
  verbosity = iConfig.getParameter<int>("verbosity");
  //TrackLabel2  = iConfig.getParameter<std::string>("TrackModuleLabel2");
  outputFileName = iConfig.getParameter<std::string>("outputFileName");
 
  CluLabel = iConfig.getParameter<std::string>("CluLabel");
  HitLabel = iConfig.getParameter<std::string>("HitLabel");
  RoadLabel = iConfig.getParameter<std::string>("RoadLabel");
  TrackLabel= iConfig.getParameter<std::string>("TrackLabel");
  RoadInstanceLabel= iConfig.getParameter<std::string>("RoadInstanceLabel");
  VtxClusterDistance= iConfig.getParameter<double>("VtxClusterDistance");
  fastSimulation= iConfig.getParameter<bool>("fastSimulation");
  //T2CutsUtil.SetCuts(T2_TrkEtamin,T2_TrkEtaMAX,T2_trkMultimin,T2_trkMultiMAX,T2_DZMultiplicator,T2_PhiChiProbCut,T2_RChiProbCut,T2_QuarterUsed,IgnoredSmallAnglePrimarySlopeCut,_T2_usesXYtracking);
  selected_event= iConfig.getParameter<int>("selected_event");

  PadRoadFinderAnalysis= iConfig.getParameter<bool>("PadRoadFinderAnalysis");
  numhitRequiredFormMatching= iConfig.getParameter<int>("numhitRequiredFormMatching");
  fastAnalysis= iConfig.getParameter<bool>("fastAnalysis");
  T2_QuarterUsed=iConfig.getParameter<std::vector<int> >("T2_QuarterUsed");
  EnergyCutinPrimaryEfficiency= iConfig.getParameter<double>("EnergyCutinPrimaryEfficiency");
  PtCutinPrimaryEfficiency= iConfig.getParameter<double>("PtCutinPrimaryEfficiency");
  ZEffiCutImpact=iConfig.getParameter<std::vector<double> >("ZEffiCutImpact");
  MaxPadCluOfFittingCurves= iConfig.getParameter<int>("MaxPadCluOfFittingCurves");
  UseselectedplanesforAPM= iConfig.getParameter<bool>("UseselectedplanesforAPM");
  NameOfGenerator= iConfig.getParameter<std::string>("NameOfGenerator");
  // qused.push_back(0);qused
  T2CutsUtil.SetCuts(4.5,7.5,4,11,2.,0.01,0.01,T2_QuarterUsed,0.001,true);
  
  // T2CutsUtil.SetCuts(T2_TrkEtamin,T2_TrkEtaMAX,T2_trkMultimin,T2_trkMultiMAX,T2_DZMultiplicator,T2_PhiChiProbCut,T2_RChiProbCut,T2_QuarterUsed,IgnoredSmallAnglePrimarySlopeCut,_T2_usesXYtracking);
  
  //VTX simu template
  //t2_quarter_used="0"
  std::cout<<"Done!"<<std::endl;
}


T2BackgroundAn::~T2BackgroundAn()
{

}



// ------------ method called to produce the data  ------------


void T2BackgroundAn::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //LOAD TRACK PRODUCER 2 TRACKS
  //Container for all the event tracks  (TRACK PRODUCER 2)  
  if(verbosity>0)
    std::cout<<"Begin Analyse"<<std::endl;
  std::auto_ptr<T1T2TrackCollection> theT2Tracks2 (new T1T2TrackCollection());
  loadTrackproducer2Tracks(iEvent, theT2Tracks2);  
  if(verbosity>0)
    std::cout<<"Begin Analyse-2"<<std::endl;
  std::auto_ptr<T1T2TrackCollection> theT2TracksRaw (new T1T2TrackCollection());
  loadTrackproducerRaw(iEvent, theT2TracksRaw);  
  

  
  numevent++;

  //LOAD SIMULATED TRACKS
  //Container for all simtracks
  std::auto_ptr<edm::SimTrackContainer> theSimTracks (new edm::SimTrackContainer());
  loadSimulatedTracks(iEvent,theSimTracks);

  //LOAD SIMULATED VERTICES
  //Container for all simvertices
  std::auto_ptr<edm::SimVertexContainer> theSimVertices (new edm::SimVertexContainer());
  loadSimulatedVertices(iEvent,theSimVertices);

  //LOAD SIMULATED HITS
  //Container for all simhits
  std::auto_ptr<edm::PSimHitContainer> theSimHits (new edm::PSimHitContainer());
  //Container for unknown tracks. Only store the trackId and particle type.
  //Unknown track means that we don't have information about the simulated track and vertice for the corresponding track ID.
  std::map<unsigned int, int> unknownTrackList;
  //Container for all hits corresponding to track.
  //Unsigned int tells track ID and list corresponding simulated hit coordinates for the track.
  std::map<unsigned int, std::list<Local3DPoint> > trackHitList;
  loadSimulatedHits(iEvent,theSimVertices,theSimTracks,theSimHits,trackHitList,unknownTrackList);

  //std::vector<unsigned int> idTrkalreadyconsidered;
  std::vector<int> idVtxAlreadyconsidered;

  edm::Handle<edm::PSimHitContainer> psimHitCollection;
  iEvent.getByLabel("g4SimHits","TotemHitsT2Gem",psimHitCollection);

  edm::Handle<edm::HepMCProduct> EvtHandle ;
  //edm::Handle<HepMCProduct> EvtHandle;

  iEvent.getByLabel("generator", EvtHandle ) ;
  const HepMC::GenEvent* evt = EvtHandle->GetEvent();

  edm::Handle<T2HitCollection> t2hitcoll; 
  iEvent.getByLabel("T2Hits","T2Hits",t2hitcoll);

   /* :::::::::::::TakeDigi::::::::::::*/

  edm::Handle<T2PadDigiCollection> t2paddigicoll;
  iEvent.getByLabel(T2PadDigiCollectionLabel, t2paddigicoll);

  t2paddigicoll.product();
  edm::Handle<T2StripDigiCollection> t2stripdigicoll;
  iEvent.getByLabel(T2StripDigiCollectionLabel, t2stripdigicoll);

  t2stripdigicoll.product();
  DigiContainerIterator<T2DetId, T2PadDigi> itp;
  DigiContainerIterator<T2DetId, T2StripDigi> its;
  
  //  process.T2TrackColl3.RoadModuleLabel="T2RoadPadFinder"
  //  process.T2TrackColl3.RoadInstanceLabel="NewRoadFinderRELOAD"
  //  event.getByLabel(RoadModuleLabel,RoadInstanceLabel,myRoadColl);

  
  edm::Handle<T2RoadCollection> NewT2RoadsPad; 
  edm::Handle<T2RoadCollection> NewT2Roads;
  if(PadRoadFinderAnalysis)
    iEvent.getByLabel(RoadLabel,RoadInstanceLabel,NewT2RoadsPad);


   edm::Handle<T2StripClusterCollection> Strclcoll;
  iEvent.getByLabel("T2MCl","T2StripClusters",Strclcoll);

  edm::Handle<T2PadClusterCollection> Padclcoll;
  iEvent.getByLabel("T2MCl","T2PadClusters",Padclcoll);

  bool showerinH1=false;
  bool sometrkexcluding64=false;

  unsigned int trkCount=0; unsigned int trkCountless64 = 0;
 unsigned int trkCountZimp7m=0;

  for(T1T2TrackCollection::const_iterator ittrack = theT2Tracks2->begin();ittrack!=theT2Tracks2->end();++ittrack){ 
   
   if(ittrack->GetHitEntries()>=4){
     
     int symb=RawtoSymb(ittrack->GetHitT2(0).GetHitDetRawId());
     symb=symb/10;
     double eta2=T2CutsUtil.EtaFromAveragesVtxCorr((*ittrack),0.,0.,0.);
     if((fabs(eta2)<6.4)&&(fabs(eta2)>5.3))
       trkCountless64++;

     double C0=(ittrack->GetHitT2(0).GetHitX()*ittrack->GetHitT2(0).GetHitX()+ittrack->GetHitT2(0).GetHitY()*ittrack->GetHitT2(0).GetHitY())/(ittrack->GetTx()*ittrack->GetHitT2(0).GetHitX()+ittrack->GetTy()*ittrack->GetHitT2(0).GetHitY());
     double Z0impact=ittrack->GetHitT2(0).GetHitZ()-C0;
     trkCount++;
     if(fabs(Z0impact)<7000.)
       trkCountZimp7m++;
   
   }
  
  }

  if(trkCountless64>0)
    sometrkexcluding64=true;

  NumTrkInT2PerEvt->Fill(trkCount);
  NumTrkInT2PerEvt7mcut->Fill(trkCountZimp7m);
 
  if(trkCount>0)
    NumEventTrackTrigger->Fill(1);
  else
    NumEventTrackTrigger->Fill(0);



  //edm::ESHandle <ParticleDataTable> pdt;
  // if(!(pdt.isValid())) 
  //iSetup.getData(pdt);

  // std::cout<<"HERE0: multiplicity calc"<<std::endl;


  int ActivePlH0[10]= {0,0,0,0,0,0,0,0,0,0};
  int ActivePlH1[10]= {0,0,0,0,0,0,0,0,0,0};
  int ActivePlH2[10]= {0,0,0,0,0,0,0,0,0,0}; 
  int ActivePlH3[10]= {0,0,0,0,0,0,0,0,0,0};

  unsigned int intplane=0;
  unsigned int totHitAllH0=0;unsigned int totHitAllH1=0;unsigned int totHitAllH2=0;unsigned int totHitAllH3=0;  
  unsigned int totHitAllH0_sel=0;unsigned int totHitAllH1_sel=0;unsigned int totHitAllH2_sel=0;unsigned int totHitAllH3_sel=0;
  unsigned int numplaneactiveH0=0;unsigned int numplaneactiveH1;
  unsigned int numplaneactiveStripH0=0; unsigned int numplaneactiveStripH1=0;
   for(T2StripClusterCollection::const_iterator itstrip = Strclcoll->begin(); itstrip != Strclcoll->end(); itstrip++){
    vector<T2Cluster> stripClv = itstrip->second;
    T2DetId *detID =new T2DetId(itstrip->first);
    uint32_t cmsswdId= detID->calculateRawId(detID->arm(),detID->halfTelescope(),detID->plane(),detID->planeSide());
    unsigned int symb=RawtoSymb(cmsswdId);
    if((symb/10)==0){
      if(itstrip->second.size()>0)
	numplaneactiveStripH0++;
    }
    if((symb/10)==1){
      if(itstrip->second.size()>0)
	numplaneactiveStripH1++;
    }
    delete detID;
  }


  for(T2PadClusterCollection::const_iterator itpad= Padclcoll->begin(); itpad != Padclcoll->end(); itpad++){//Reft2padclcoll
    vector<T2Cluster> padClv = itpad->second;
    T2DetId *detID =new T2DetId(itpad->first);
    
    uint32_t cmsswdId= detID->calculateRawId(detID->arm(),detID->halfTelescope(),detID->plane(),detID->planeSide());
    unsigned int symb=RawtoSymb(cmsswdId);
    //unsigned int quarter=symbol/10;
    //unsigned int plane=symbol%10;
    intplane=55;

    //To avoid possible bias due to dead vfat
    // H0: 2-4-7 H1:12-14-17 H2:23-24-27 H3: 32-34-37
   
    PlaneClusterActiveCounter->Fill(symb);
    if((symb/10)==0){
	totHitAllH0=totHitAllH0+padClv.size();
	intplane=(symb%10);
	ActivePlH0[intplane]=ActivePlH0[intplane]+padClv.size();

	if(((symb%10)==2)||((symb%10)==4)||((symb%10)==7))
	  totHitAllH0_sel=totHitAllH0_sel+padClv.size();

	if(padClv.size()>0)
	  numplaneactiveH0++;
      }
      if((symb/10)==1){	
	
	if(((symb%10)==2)||((symb%10)==4)||((symb%10)==7))
	  totHitAllH1_sel=totHitAllH0_sel+padClv.size();

	totHitAllH1=totHitAllH1+padClv.size();
	intplane=(symb%10);
	ActivePlH1[intplane]=ActivePlH1[intplane]+padClv.size();
	

	if(padClv.size()>0)
	  numplaneactiveH1++;
      }
      if((symb/10)==2){
	totHitAllH2=totHitAllH2+padClv.size();
	intplane=(symb%10);
	ActivePlH2[intplane]=ActivePlH2[intplane]+padClv.size();

	if(((symb%10)==3)||((symb%10)==4)||((symb%10)==7))
	  totHitAllH2_sel=totHitAllH0_sel+padClv.size();

      }
      if((symb/10)==3){
	totHitAllH3=totHitAllH3+padClv.size();     
	intplane=(symb%10);
	ActivePlH3[intplane]=ActivePlH3[intplane]+padClv.size();
	if(((symb%10)==2)||((symb%10)==4)||((symb%10)==7))
	  totHitAllH3_sel=totHitAllH0_sel+padClv.size();
      }  
      delete detID;
 }

  AveragePadCLSDistrH0->Fill(totHitAllH0/10.);
  AveragePadCLSDistrH1->Fill(totHitAllH1/10.);


  if(trkCount>0){
    numtotevent++;
    Histonumtoeventgen->Fill(1);
    if((totHitAllH1/10.)>66.){
      showerinH1=true;
      numshowevent++;
      Histonumshowevent->Fill(1);
    }else{
      numcleanevent++;
      Histonumcleanevent->Fill(1);
    }
  }


  // if(PadRoadFinderAnalysis)
  //iEvent.getByLabel(RoadLabel,"T2RoadColl",NewT2Roads); 
  std::vector<double> PrimaryTrkenergy_ForTheRecoHits;
  //  std::cout<<"Numevent:"<<numevent<<std::endl;
  
  
  T2GeometryUtil conv;     

  
  std::map<int, unsigned int> countsinsamerowH0U;countsinsamerowH0U.clear();
  std::map<int, unsigned int> countsinsamerowH0D;countsinsamerowH0D.clear();
  std::map<int, unsigned int> countsinsamerowH1U;countsinsamerowH1U.clear();
  std::map<int, unsigned int> countsinsamerowH1D;countsinsamerowH1D.clear();  
  // int minicount=0;



 



  // BEGIN of the SIMULATION TRACKs
  vectorEta2CounterFakePrimary.clear();
  vectorEta2CounterDnDeta.clear();
  vectorMBALLMCGenerator.clear();
  vectorEta2CounterRecoDnDeta2UnfoldCut2.clear();  
  vectorEta2CounterRecoDnDeta2UnfoldCut2_Primary.clear();
  vectorEta2CounterRecoDnDeta2UnfoldCut.clear();  
  vectorEta2CounterRecoDnDeta2UnfoldCut_Primary.clear();
  
  vectorEta2CounterDnDetaPrimary.clear();
  vectorEta2SimulatedPrimaryTrkInT2_H0.clear(); 
  vectorEta2SimulatedPrimaryTrkInT2_H1.clear();
  vectorEta2CounterDnDetaNoCut.clear();
  vectorEta2CounterDnDetaPrimaryNoCut.clear();
  //  std::cout<<"A"<<std::endl;
  vectorEta2GeantPrimaryTrkInH0_ZCut.clear();
  vectorEta2PrimaryGeantDNDeta2_IfOneRecoZImpInH0.clear();
  for(unsigned int iu=0;iu<totnumberbinfordndeta+1;iu++)
    {
      vectorEta2GeantPrimaryTrkInH0_ZCut.push_back(0);
      vectorEta2CounterFakePrimary.push_back(0);
      vectorEta2CounterDnDeta.push_back(0);
      vectorEta2CounterDnDetaPrimary.push_back(0);vectorEta2CounterDnDetaNoCut.push_back(0);
      vectorEta2SimulatedPrimaryTrkInT2_H0.push_back(0);
      vectorEta2SimulatedPrimaryTrkInT2_H1.push_back(0);
      vectorEta2CounterDnDetaPrimaryNoCut.push_back(0);

      vectorEta2CounterRecoDnDeta2UnfoldCut.push_back(0);  
      vectorEta2CounterRecoDnDeta2UnfoldCut_Primary.push_back(0); 
      vectorEta2CounterRecoDnDeta2UnfoldCut2.push_back(0);
      vectorEta2CounterRecoDnDeta2UnfoldCut2_Primary.push_back(0);
      vectorMBALLMCGenerator.push_back(0);
      vectorEta2PrimaryGeantDNDeta2_IfOneRecoZImpInH0.push_back(0);
    }

  

   double GeantTrkEffCountPlus=0.;

   
   double GeantTrkEffCountMinus=0.;

   unsigned int Garm=22; unsigned int GQuarter=22;
   unsigned int Garm0=22; 
   unsigned int Garm1=22; 
   unsigned int GQuarter0=22;unsigned int GQuarter1=22;  unsigned int GQuarter2=22;unsigned int GQuarter3=22;
   std::map<unsigned int, std::vector<int> > trkIdPadCol_map;
   std::map<unsigned int, std::vector<int> > trkIdPadRow_map;
   std::map<unsigned int, std::vector<int> > trkIdPlane_map;

  
 


   




   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   unsigned int   chp12pos=0;
   unsigned int   chp12neg=0;
   unsigned int   chp12_5364=0;
   
   
   for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
	   
     //const ParticleData * part = pdt->particle((*p)->pdg_id());
     //int charge = part->charge();

    //std::cout<<"Particle PID: "<<(*p)->pdg_id()<<"  charge:"<<charge<<std::endl;
     int charge=PartCharge((*p)->pdg_id());  //OBSOLETE WAY
    
    
      //      if((*p)->momentum().e()>EnergyCutinPrimaryEfficiency)
     if((abs(charge)>0)&&((*p)->status()==1))
      if((*p)->momentum().perp()>PtCutinPrimaryEfficiency)
	{	      	    	    
	double etag=(*p)->momentum().eta();

	if(trkCount>0){
	  DNDetaT2ALL->Fill(etag);
	  
	  if(showerinH1==false)
	    DNDetaT2Clean->Fill(etag);
	  if(showerinH1)
	    DNDetaT2Show->Fill(etag);
	}

	if(fabs(etag)<maxetafordndetahisto)
	  {
	    if((fabs(etag)<6.5)&&(fabs(etag)>5.3))
	      {
		//StableChPartilceEnergy->Fill((*p)->momentum().e());
		if(etag>0)
		  chp12pos++; 
		else
		  chp12neg++;

		if(fabs(etag)<6.4)
		  chp12_5364++;
	      }
	    int trkbinfordndeta_=(int)((etag-(-maxetafordndetahisto))/etabinsize); 
	    vectorMBALLMCGenerator.at(trkbinfordndeta_)= vectorMBALLMCGenerator.at(trkbinfordndeta_)+1;
	    
	    //vectorMBALLMCGenerator_RecoAccepted.at(trkbinfordndeta_)= vectorMBALLMCGenerator_RecoAccepted.at(trkbinfordndeta_)+1;
	  }
      }
  }
   
   bool generator_Inelastic5364=false;
   if(chp12_5364>0)
     generator_Inelastic5364=true;

   bool generator_Inelastic=false;
   if((chp12pos>0)||(chp12neg>0)){
     generator_Inelastic=true;
     NumEventGeneratorTrigger->Fill(1);
   }else
     NumEventGeneratorTrigger->Fill(0);

   NumChPartInT2PerEvt->Fill(chp12pos+chp12neg);

   bool AtLeastOneTrkZimpInH0=false;
 

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //   std::cout<<"----------New EVENT: Begin of main SimTrk loop"<<std::endl;

   
	  double theX=0.;
	  double theY=0.;
	  double theR=0.;
	  double theZ=0.;

   
   unsigned int Trkhitcounter=0;
   unsigned int TrkCounter=0;int trackId=0;bool trkassociated=false;unsigned int particleType=0;
   double trkenergy=0.1;  double trkpt=0.01;


    std::vector<std::pair<double,double> > SimHitRZ_ofGeantTrk_H0;
    std::vector<std::pair<double,double> > SimHitRZ_ofGeantTrk_H1;
    std::vector<std::pair<double,double> > SimHitRZ_ofGeantTrk_ForEff;
    std::vector<unsigned int> PadsRowForSameTrk;
    std::vector<unsigned int> PadsColForSameTrk;
    std::vector<unsigned int> StripRowForSameTrk;
    std::vector<unsigned int> StripColForSameTrk;
    T2GeometryUtil::T2DetInfo detInfo2;
    unsigned int NumT2HitinGeantTrkPrimary=0;



    SimTrack PythiaIPChParticleTrack;
   int motherbarcode=99;

    
    bool trkprimcond=false;

    int GeneratorBarcode=0;    
    int ParenttrackId=0; 
    unsigned int detUnitId =0; 
    T2Geometry t2geom;Local3DPoint localPosition ;

for(edm::SimTrackContainer::const_iterator ittrack = theSimTracks->begin();ittrack != theSimTracks->end(); ++ittrack){
  
  trackId = ittrack->trackId();
  //std::map<unsigned int, std::list<Local3DPoint> >::iterator checkIsT2Trk;
  Trkhitcounter=0;
  TrkCounter++;
  Garm=22; GQuarter=22;
  Garm0=0; 
  Garm1=0; 
  GQuarter0=0; GQuarter1=0; GQuarter2=0; GQuarter3=0;


  // std::cout<<"New Gtrack Analysis:"<<std::endl;
  // std::cout<<"a)Geant Trk PID Printing: "<<std::endl;

  
    //This Trk Id has some PSimHit Found and saved in trackHitList
    //If It is not true, like in the SimHit generated in Pi0 Gun, you cannot do more

    trkassociated=false;
    particleType=0;
    trkenergy=0.1; 
    trkpt=0.01;
    trkpt=ittrack->momentum().pt();
    trkenergy=ittrack->momentum().e();
    //  std::cout<<"Gtrk E:"<<trkenergy<<" PID:"<<ittrack->type()<<std::endl;
    SimHitRZ_ofGeantTrk_H0.clear();
    SimHitRZ_ofGeantTrk_H1.clear();
    SimHitRZ_ofGeantTrk_ForEff.clear();
    PadsRowForSameTrk.clear();
    PadsColForSameTrk.clear();
    StripRowForSameTrk.clear();
    StripColForSameTrk.clear();
    
    NumT2HitinGeantTrkPrimary=0;

    //  std::cout<<"B"<<std::endl;
    // std::cout<<" trackId: "<<trackId<<std::endl;


    
    
     PythiaIPChParticleTrack=(*ittrack);
      motherbarcode=99;

    
    // thePidGenerator=GetTrkOlderMotherPid(PythiaIPChParticleTrack,theSimVertices,theSimTracks,evt,motherbarcode);
    //bool isHitPrimary = isTrackPrimary(trackId, theSimTracks, theSimVertices);     	        
    //HoSimHitActiveplane
   
     trkprimcond=PrimaryTrackCondition(PythiaIPChParticleTrack,theSimVertices,theSimTracks,evt,motherbarcode);
    // std::cout<<k:"<<trkprimcond<<std::endl;

    if(trkprimcond){
      // std::cout<<"New Prim Track analysis"<<std::endl;

       GeneratorBarcode=PythiaIPChParticleTrack.genpartIndex();
       ParenttrackId=0; 
       detUnitId =0; 
      
      
    for(edm::PSimHitContainer::const_iterator ithit = psimHitCollection->begin();ithit != psimHitCollection->end(); ++ithit){
	        
       ParenttrackId = ithit->trackId();  
      // if(TrkCounter==1)
      //std::cout<<"ALL HitTrkID vs first TrkID:"<<ParenttrackId<<"-"<<  trackId<<std::endl;
      //      std::cout<<"New Prim Track analysis-A"<<std::endl;
      if(ParenttrackId == trackId)
	{
	  // T2GeometryUtil::T2DetInfo detInfo2;
 
	  Trkhitcounter++; 
	  NumT2HitinGeantTrkPrimary++;

	  localPosition = ithit->localPosition();   
	  detUnitId = ithit->detUnitId();
	  t2geom.setPlane(detUnitId);
	  detInfo2=conv.GetT2Info(detUnitId);
	  particleType = ithit->particleType();//
	  //std::cout<<" Hit Type: "<<particleType<<std::endl;
	  SimVertex simuvertex; 
	   
	
	      
	  PrimTrkTestedPlaneQuarter->Fill(detInfo2.symb);

	  // std::cout<<"GHit in plane:"<<detInfo2.symb<<std::endl;
	  
	  if(detInfo2.symb>=20){
	    Garm=1;
	    Garm1++;
	    if(detInfo2.symb>=30){
	      GQuarter=3;
	      GQuarter3++;
	    }
	    else{
	      GQuarter=2;
	      GQuarter2++;
	    }
	  }
	  else{
	    Garm0++;
	    Garm=0;
	    if(detInfo2.symb>=10){
	      H1PrimaryTrkPdgID->Fill(particleType);
	      H1PrimaryTrkE->Fill(trkenergy);
	      GQuarter=1;
	      GQuarter1++;
	    }
	    else{
	      H0PrimaryTrkPdgID->Fill(particleType);
	      H0PrimaryTrkE->Fill(trkenergy);
	      GQuarter=0;
	      GQuarter0++;
	    }
	  }

	  
	  Local3DPoint corPos = CalculateCorrectXYZPos(Local3DPoint(localPosition.x(),localPosition.y(),localPosition.z()),detUnitId);

	  //int padCol = t2geom.getNearestPadCo(&corPos);	  
	  //  int padRow = t2geom.getNearestPadRo(&corPos);
	  //int stripCol = t2geom.getNearestStripCo(&corPos); 
	  //int stripRow = t2geom.getNearestStripRo(&corPos);
	  // std::pair<double, double> Padrc(padRow,padCol);
	  
	  if(GQuarter==0){
	    corposprimaryXYQ0->Fill(corPos.x(),corPos.y());
	  }
	  if(GQuarter==1){
	    corposprimaryXYQ1->Fill(corPos.x(),corPos.y());
	  }
	  //std::cout<<"New Prim Track analysis-B1"<<std::endl;
	  theX=localPosition.x();
	  theY=localPosition.y();
	  theR=sqrt(theX*theX+theY*theY);
	  //std::cout<<"Hit xyr:"<<theX<<" "<<theY<<" "<<theR<<std::endl;
	  
	  theZ=detInfo2.Zdet;
	  std::pair<double, double> SimHitRZ_PrimTrk_ForEff(theR,theZ);
	  SimHitRZ_ofGeantTrk_ForEff.push_back(SimHitRZ_PrimTrk_ForEff);
	  //if (getVertex(ittrack->vertIndex(), theSimVertices,simuvertex) == true) //This vertex exists
	  //if (simuvertex.parentIndex() < 0){
	  // Local3DPoint localPositione = ithit->localPosition();		  		  	       
	  std::pair<double, double> SimHitRZ_PrimTrk(theR,theZ);
	  if((theZ>0)&&((detInfo2.symb/10)==0))
	    SimHitRZ_ofGeantTrk_H0.push_back(SimHitRZ_PrimTrk);
	  if((theZ>0)&&((detInfo2.symb/10)==1))
	    SimHitRZ_ofGeantTrk_H1.push_back(SimHitRZ_PrimTrk);
	  	
	  // std::cout<<"New Prim Track analysis-B2"<<std::endl;
	  
	     	      	      	      
	}//if(ParenttrackId == trackId)

      //std::cout<<"New Prim Track analysis-C"<<std::endl;  
   	  
    }//end simHit loop Now you are again in the simTrk loop
	
	
    //  std::cout<<"SimHitLoopEnd"<<std::endl;

    if((trkassociated==false)&&(psimHitCollection->size()>0))
      ParticlePdgvsENotinT2->Fill(particleType,trkenergy);
  
    //idTrkalreadyconsidered.push_back(trackId);
    //  std::cout<<"GTrack associated in Q: "<<GQuarter<<std::endl;

    double eta2geant=0;
    //if(trkenergy>EnergyCutinPrimaryEfficiency)
      
    if(trkpt>PtCutinPrimaryEfficiency){
    if((SimHitRZ_ofGeantTrk_H0.size()>=4)&&(NumT2HitinGeantTrkPrimary>=4))
      {
      	eta2geant=GetEta2FromGeantHits(SimHitRZ_ofGeantTrk_H0);
	PrimaryGeanteta2->Fill(eta2geant);     
	double trkbinfordndeta2=(eta2geant-(-maxetafordndetahisto))/etabinsize; 
	//  std::cout<<"Geant eta:"<<eta2geant<<std::endl;
	vectorEta2SimulatedPrimaryTrkInT2_H0.at(trkbinfordndeta2)=vectorEta2SimulatedPrimaryTrkInT2_H0.at(trkbinfordndeta2)+1;	
	//This Geant DN/Deta is for one quarter only 
      }

    if((SimHitRZ_ofGeantTrk_H1.size()>=4)&&(NumT2HitinGeantTrkPrimary>=4))
      {
      	eta2geant=GetEta2FromGeantHits(SimHitRZ_ofGeantTrk_H1);
	PrimaryGeanteta2->Fill(eta2geant);     
	double trkbinfordndeta2=(eta2geant-(-maxetafordndetahisto))/etabinsize; 
	//  std::cout<<"Geant eta:"<<eta2geant<<std::endl;
	vectorEta2SimulatedPrimaryTrkInT2_H1.at(trkbinfordndeta2)=vectorEta2SimulatedPrimaryTrkInT2_H1.at(trkbinfordndeta2)+1;	
	//This Geant DN/Deta is for one quarter only 
      }
    }

    int  howmanyRecoTrk=0; 

    double eta2geant_forEff=GetEta2FromGeantHits(SimHitRZ_ofGeantTrk_ForEff);
    eta2geant_forEffHisto1->Fill(eta2geant_forEff);
//    double trkbinfordndeta2=(eta2geant_forEff-(-maxetafordndetahisto))/etabinsize;
    
    GeneratorEtaFromTrkBarcode(GeneratorBarcode,evt);

    //Warning: you are still inside the SimTrk loop
    //If the simtrack is primary and have at least 4 hit in T2

    
    bool trkEffiEnergyCondition=false;
    // if(trkenergy>EnergyCutinPrimaryEfficiency)
    //trkEffiEnergyCondition=true;
    
    if(trkpt>PtCutinPrimaryEfficiency)
      trkEffiEnergyCondition=true;
    //std::cout<<"PrimTRK?:"<<trkprimcond<<" Etrk:"<<trkenergy<<" Ecut:"<<EnergyCutinPrimaryEfficiency<<" #PrimGeantHit:"<<NumT2HitinGeantTrkPrimary<<std::endl;
    if(trkEffiEnergyCondition)
    if(NumT2HitinGeantTrkPrimary>=4) //Calculate trk Reco Effi;
      {
	//Criteria: There will be 3 or more reco hit which is mathcing Geant
	if(Garm==0){
	  GeantTrkEffCountPlus+=1;
	}
	else{
	  GeantTrkEffCountMinus+=1;
	}
	double lastpassingZImpact=15000.;
	double lastPassingZmin=15000.;
	unsigned int NumGeantTrk_PassingCuts=0;	
	unsigned int NumGeantTrk_PassingCuts_56Division=0;
	unsigned int NumGeantTrk_PassingCuts_2=0;
	unsigned int NumGeantTrk_PassingCutsDiscrepancyZCutCond=0; unsigned int NumGeantTrk_PassingCutsDiscrepancy=0;
	T1T2Track TheGeantTrkReconstr2;double effiasRecRec=-1.0;bool GeantTrkFound2=false;
	//std::cout<<"GeantTest-Begin "<<std::endl;
	//
	int BestTrkIdMatchedPosition=-1;
	eta2geant_forEffHisto2->Fill(eta2geant_forEff);
	//This are the number of track associated to the primaryGeant, while if the cut is passed is written in NumGeantTrk_PassingCuts 
	howmanyRecoTrk=GeantTrkReconstructed(iEvent,theSimHits,theSimTracks,theSimVertices,trackHitList,theT2Tracks2,trackId,false,NumGeantTrk_PassingCuts,NumGeantTrk_PassingCuts_56Division,NumGeantTrk_PassingCutsDiscrepancyZCutCond,NumGeantTrk_PassingCutsDiscrepancy,NumGeantTrk_PassingCuts_2/*,ThePrimaryGeantHitPos*/,TheGeantTrkReconstr2,GeantTrkFound2,effiasRecRec,BestTrkIdMatchedPosition,lastPassingZmin,lastpassingZImpact,true); //Last one are for discrepancy studies    
	double Phigeant_forEff=0.;
	if(howmanyRecoTrk!=-1)
	  Phigeant_forEff=TheGeantTrkReconstr2._phiRZ*180./3.14159;
	
	if(verbosity>1)
	  std::cout<<"After GeantTest: howmanyrecotrk :  "<<howmanyRecoTrk<<std::endl;

	//	std::cout<<"TheGeantTr"<<std::endl;

	/*	
	//-------------RESOLUTION STUDIES BEGIN---------------
	if(BestTrkIdMatchedPosition>=0){
	  
	  T1T2TrackCollection::const_iterator recTraccia =(theT2Tracks2->begin()+BestTrkIdMatchedPosition);
	  double eta2MatchingTrk=T2CutsUtil.EtaFromAveragesVtxCorr((*recTraccia),0.,0.,0.);	 
	  EtaResolGeant->Fill(etaGenerator,fabs(eta2MatchingTrk-etaGenerator));
	  EtaResolGenerator->Fill(eta2geant_forEff,fabs(eta2MatchingTrk-eta2geant_forEff));
	}
	//Just for PGun studies: 1 simTrk, 1 reco trk. See the eta2 difference. 
	if(abs(theT2Tracks2->begin()-theT2Tracks2->end())==1)
	  if(abs(theSimTracks->begin()-theSimTracks->end())==1){
	    T1T2TrackCollection::const_iterator recTraccia =(theT2Tracks2->begin());
	    double eta2MatchingTrk=T2CutsUtil.EtaFromAveragesVtxCorr((*recTraccia),0.,0.,0.);
	    EtaResolFreeGeant->Fill(etaGenerator,fabs(eta2MatchingTrk-etaGenerator));
	    EtaResolFreeGenerator->Fill(eta2geant_forEff,fabs(eta2MatchingTrk-eta2geant_forEff));
	}
	//-------------RESOLUTION STUDIES END---------------
	*/
	bool unambigous_Quarter_association=true;
	if(((GQuarter0>0)&&(GQuarter1>0))||((GQuarter2>0)&&(GQuarter3>0)))
	  unambigous_Quarter_association=false;

	if(unambigous_Quarter_association){
	if(GQuarter0>0)
	  eta2geant_forEffQ0->Fill(eta2geant_forEff);	    
	
	if(GQuarter1>0)
	  eta2geant_forEffQ1->Fill(eta2geant_forEff);	    
	
	if(GQuarter2>0)
	  eta2geant_forEffQ2->Fill(eta2geant_forEff);	    
	
	if(GQuarter3>0)
	  eta2geant_forEffQ3->Fill(eta2geant_forEff);
	}


	if(unambigous_Quarter_association){
	  eta2geant_forEffHisto3->Fill(eta2geant_forEff);
	  if(howmanyRecoTrk>=0){//Could be -1 in that case you didn't search for the matching reco trk
	    if(verbosity>1)
	      std::cout<<"Here00"<<std::endl;
	    
	    double ZImpGeant=ZimpactFromRecTrack(TheGeantTrkReconstr2);
	    ZImpOfAPrimaryGeantTrk_Cumulative->Fill(ZImpGeant);
	    double ZMin=TheGeantTrkReconstr2.Z_at_Rmin();
	    if((fabs(ZMin)<13500)||(ZMin*eta2geant_forEff<0))
	      ZImpOfAPrimaryGeantTrk_CumulativeZMinCut->Fill(ZImpGeant);	    
	  }
	}
	
	if(verbosity>1)
	  std::cout<<"HereA"<<std::endl;
	
	PrimTrkTestedQuarter->Fill(GQuarter); 
	if(unambigous_Quarter_association)
	  PrimTrkTestedQuarterUnamb->Fill(GQuarter);

	if(verbosity>1)
	  std::cout<<"HereB"<<std::endl;

	double AvgmultclassPlus=((totHitAllH0+totHitAllH1)/20.);
	double AvgmultclassMinus=((totHitAllH2+totHitAllH3)/20.);
	
	if(howmanyRecoTrk>0){
	 /* 
	  if((unambigous_Quarter_association)){
	  T1T2TrackCollection::const_iterator recTraccia =(theT2Tracks2->begin()+BestTrkIdMatchedPosition);
	  double eta2MatchingRecTrk=T2CutsUtil.EtaFromAveragesVtxCorr((*recTraccia),0.,0.,0.);
	  T1T2Track MyRecTrk=(*recTraccia);
	  double ZImpRec=ZimpactFromRecTrack(MyRecTrk);
	  ZImpOfAPrimaryGeantTrk_Cumulative_match->Fill(ZImpRec);
	  double ZMinRec=(*recTraccia).Z_at_Rmin();
	   if((fabs(ZMinRec)<13500)||(ZMinRec*eta2MatchingRecTrk<0))
	    ZImpOfAPrimaryGeantTrk_CumulativeZMinCut_match->Fill(ZImpRec);
	  }
	*/

	  if(Garm0>0)
	    ArmPlus_TrkEtaEfficiencyCumulative->Fill(AvgmultclassPlus,1);
	  if(Garm1>0)
	    ArmMinus_TrkEtaEfficiencyCumulative->Fill(AvgmultclassMinus,1);
	}

	if(howmanyRecoTrk==0){
	  if(verbosity>1)
	    std::cout<<"HereB1"<<std::endl;
	   if(Garm0>0)
	     ArmPlus_TrkEtaEfficiencyCumulative->Fill(AvgmultclassPlus,0);
	   if(Garm1>0)
	    ArmMinus_TrkEtaEfficiencyCumulative->Fill(AvgmultclassMinus,0);
	}


	//This means a trk is associated to the geant pad but you don't know if it will pass the cuts.
	if(howmanyRecoTrk>0){
	T1T2TrackCollection::const_iterator recTraccia =(theT2Tracks2->begin()+BestTrkIdMatchedPosition);
	if(verbosity>1)
	  std::cout<<"HereD"<<std::endl;
	double eta2MatchingTrk=T2CutsUtil.EtaFromAveragesVtxCorr((*recTraccia),0.,0.,0.);

	 if(verbosity>1)
	  std::cout<<"Heredfr"<<std::endl;
	if((fabs(eta2MatchingTrk))>5.35)
	  if((fabs(eta2MatchingTrk))<6.4)
	    {
	      int BinAccordingTo_eta2 = (int)((fabs(eta2MatchingTrk)-5.35)/0.05);//21 Bin
	      if((BinAccordingTo_eta2>=0)&&(BinAccordingTo_eta2<=21)){
		if(GQuarter0>0)
		  Associated_PrimaryTrkZImpact_H0[BinAccordingTo_eta2]->Fill(lastpassingZImpact);
		
		if(GQuarter1>0)
		  Associated_PrimaryTrkZImpact_H1[BinAccordingTo_eta2]->Fill(lastpassingZImpact);
		
		if(GQuarter2>0)
		  Associated_PrimaryTrkZImpact_H2[BinAccordingTo_eta2]->Fill(lastpassingZImpact);
		
		if(GQuarter3>0)
		  Associated_PrimaryTrkZImpact_H3[BinAccordingTo_eta2]->Fill(lastpassingZImpact);		     	  		     
	      }
	    }
	}

	unsigned int intmultclass=0;
	double multclass=0.;
	if(verbosity>1)
	  std::cout<<"Here3. Mults:"<<totHitAllH0<<" "<<totHitAllH1<<" "<<totHitAllH2<<" "<<totHitAllH3<<std::endl;

	if((unambigous_Quarter_association)&&(howmanyRecoTrk!=-1)){//-1 means does not hav enogh hits or there is a bug
	if(NumGeantTrk_PassingCuts>=1)
	  {
	    ZImpOfAPrimaryGeantTrk_Cumulative_match->Fill(lastpassingZImpact);
	    if((fabs(lastPassingZmin)<13500)||(lastPassingZmin*eta2geant_forEff<0))
	      ZImpOfAPrimaryGeantTrk_CumulativeZMinCut_match->Fill(lastpassingZImpact);

	    if(GQuarter0>0)
	      eta2geant_forEffQ0eff->Fill(eta2geant_forEff);	    
	    
	    if(GQuarter1>0)
	      eta2geant_forEffQ1eff->Fill(eta2geant_forEff);	    
	    
	    if(GQuarter2>0)
	      eta2geant_forEffQ2eff->Fill(eta2geant_forEff);	    
	    
	    if(GQuarter3>0)
	      eta2geant_forEffQ3eff->Fill(eta2geant_forEff);
	    
	    // To avoid bias due to missing vfat you can also do:
	    // H0: 2-4-7 H1:12-14-17 H2:23-24-27 H3: 33-34-37

	    if(Garm0>0){	     
	      //TrkRecoSuccessPlus+=1;
	      if(GQuarter0>0){
		H0_TrkEtaEfficiencyCumulative_UnfoldCuts->Fill(eta2geant_forEff,1);
		
		multclass=(totHitAllH0/10.);
		
		if(UseselectedplanesforAPM)
		  multclass=(totHitAllH0_sel/3.0);
		
		//class are in steps of 5 pad-clusters
		intmultclass=(unsigned int)floor(multclass/5.);
		H0_TrkPhiEfficiencyCumulative_UnfoldCutsCumul->Fill(Phigeant_forEff,1);

		if(intmultclass>9)//High multiplicity event goes in the same bin;
		  intmultclass=9;

		H0_BinMultStat->Fill(intmultclass);
		H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(eta2geant_forEff,1); 
		H0_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(Phigeant_forEff,1);
	
	      }

	      if(GQuarter1>0){
		H1_TrkEtaEfficiencyCumulative_UnfoldCuts->Fill(eta2geant_forEff,1);
		multclass=(totHitAllH1/10.);
		if(UseselectedplanesforAPM)
		  multclass=(totHitAllH1_sel/3.0);

		//class are in steps of 5 pad-clusters
		intmultclass=(unsigned int)floor(multclass/5.);
		H1_TrkPhiEfficiencyCumulative_UnfoldCutsCumul->Fill(Phigeant_forEff,1);


		if(intmultclass>9)//High multiplicity event goes in the same bin;
		  intmultclass=9;		
		
		H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(eta2geant_forEff,1);
		H1_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(Phigeant_forEff,1);
		H1_BinMultStat->Fill(intmultclass);
		/*
		if(NumGeantTrk_PassingCutsDiscrepancyZCutCond>=1){
		  H1TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut->Fill(multclass,1);
		}
		if(NumGeantTrk_PassingCutsDiscrepancy>=1){
		  H1TrkEffi_AxisAvgMult_ZImpactCut->Fill(multclass,1);
		}
		*/
	      }

	      ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts->Fill(eta2geant_forEff,1);
	    }

	   if(Garm1>0){
	      // TrkRecoSuccessMinus+=1;
	      if(GQuarter2>0){
		H2_TrkEtaEfficiencyCumulative_UnfoldCuts->Fill(eta2geant_forEff,1);
		multclass=(totHitAllH2/10.);
		if(UseselectedplanesforAPM)
		  multclass=(totHitAllH2_sel/3.0);

		//class are in steps of 5 pad-clusters
		intmultclass=(unsigned int)floor(multclass/5.);
		if(intmultclass>9)//High multiplicity event goes in the same bin;
		  intmultclass=9;		
		H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(eta2geant_forEff,1);
		if((Phigeant_forEff>=37.)&&(Phigeant_forEff<=329.))
		  H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_PhiCut[intmultclass]->Fill(eta2geant_forEff,1);

		H2_BinMultStat->Fill(intmultclass);
		/*
		if(NumGeantTrk_PassingCutsDiscrepancyZCutCond>=1){
		  H2TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut->Fill(multclass,1);
		}
		if(NumGeantTrk_PassingCutsDiscrepancy>=1){
		  H2TrkEffi_AxisAvgMult_ZImpactCut->Fill(multclass,1);
		}
		*/
	      }
	      if(GQuarter3>0){
		H3_TrkEtaEfficiencyCumulative_UnfoldCuts->Fill(eta2geant_forEff,1);
		multclass=(totHitAllH3/10.);
		if(UseselectedplanesforAPM)
		  multclass=(totHitAllH3_sel/3.0);
		//class are in steps of 5 pad-clusters
		intmultclass=(unsigned int)floor(multclass/5.);
		if(intmultclass>9)//High multiplicity event goes in the same bin;
		  intmultclass=9;		
		H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(eta2geant_forEff,1);
		H3_BinMultStat->Fill(intmultclass);
		/*
		if(NumGeantTrk_PassingCutsDiscrepancyZCutCond>=1){
		  H3TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut->Fill(multclass,1);
		}
		if(NumGeantTrk_PassingCutsDiscrepancy>=1){
		  H3TrkEffi_AxisAvgMult_ZImpactCut->Fill(multclass,1);
		}
		*/
	      }

	      ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts->Fill(eta2geant_forEff,1);	     	      
	    }
	    
	  }//If NumGeantTrk_PassingCuts>=1 end

	if(NumGeantTrk_PassingCuts==0){
	  if(Garm0>0){
	    ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts->Fill(eta2geant_forEff,0);

	    if(GQuarter0>0){
	      H0_TrkEtaEfficiencyCumulative_UnfoldCuts->Fill(eta2geant_forEff,0);
	      H0_TrkPhiEfficiencyCumulative_UnfoldCutsCumul->Fill(Phigeant_forEff,0);

	       multclass=(totHitAllH0/10.);
	       	
	      //class are in steps of 5 pad-clusters
	       intmultclass=(unsigned int)floor(multclass/5.);
	      if(intmultclass>9)//High multiplicity event goes in the same bin;
		intmultclass=9;		
	      
	    

	      H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(eta2geant_forEff,0);
	      H0_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(Phigeant_forEff,0);
	      H0_BinMultStat->Fill(intmultclass);
	      H0_BinMultStatFail->Fill(intmultclass);
	     
	      // if(NumGeantTrk_PassingCutsDiscrepancyZCutCond>=1)
	      //H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_ZcutOnGeant[intmultclass]->Fill(eta2geant_forEff,0);
	      
	      //if(NumGeantTrk_PassingCutsDiscrepancyZCutCond>=1)
	      /*
	      if(intmultclass<14)
	      {
		//	H0TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut->Fill(multclass,0);

		//BEGIN DEBUG THE PROCESS
		unsigned int numtrkinh0passingCuts=0;
		unsigned int numtrkinh0=0;
		double minDx=400.;double minDy=400.;double minDistance=400.;	      
		
		  {

		    double Zimpgeanteffifail=ZimpactFromRecTrack(TheGeantTrkReconstr2);
		      if(verbosity>1)
			std::cout<<"!!!!!!!!!!!!!! Inefficiency, filling H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult:"<<intmultclass<<" GeantTrK ZImp "<<Zimpgeanteffifail<< "#Strip: "<<numplaneactiveStripH0<<std::endl;
		    unsigned int numtrkinh0=0; 
		    std::vector<int> OneHalf;OneHalf.push_back(0);
		    ZImpOfAGeantTrkWhenEffiFailH0->Fill(Zimpgeanteffifail);
		    for(T1T2TrackCollection::const_iterator recTrack = theT2Tracks2->begin(); recTrack!=theT2Tracks2->end();++recTrack){
		      T1T2Track debtrk=(*recTrack);

		      double zimprectrackeffifail=ZimpactFromRecTrack(debtrk);
		     if(verbosity>1)
		       std::cout<<"RecTrk with ZImp"<<zimprectrackeffifail<<std::endl;
		      if(RecoTrk_UnfoldingCut(debtrk))
			if(T2CutsUtil.TrkAlsoInQuarter(debtrk,OneHalf)){		  
			  numtrkinh0++;
			  ZImpOfARecTrkWhenEffiFail->Fill(zimprectrackeffifail);
			  if(RecoTrk_UnfoldingCut(debtrk))
			    numtrkinh0passingCuts++;
			  
			  double thisdx=(TheGeantTrkReconstr2.GetHitT2(0).GetHitX()-recTrack->GetHitT2(0).GetHitX());
			  double thisdy=(TheGeantTrkReconstr2.GetHitT2(0).GetHitY()-recTrack->GetHitT2(0).GetHitY());
			  double Distance=sqrt(thisdx*thisdx+thisdy*thisdy);
			  if(Distance<minDistance){
			    minDx=thisdx;
			    minDy=thisdy;
			    minDistance=Distance;
			    
			  }
			}
		    }
		    Numtrk_vs_NumtrkpasCutsWhenEffiFail->Fill(numtrkinh0,numtrkinh0passingCuts);
		    
		    minDxDyWhenGeantEffiFail->Fill(minDx,minDy);//minPhi  minR
		    NumH0TrkWhenGeantEffiFail->Fill(numtrkinh0);		  
		    NumH0PlaneActiveWhenGeantEffiFail->Fill(numplaneactiveH0);
		    NumH0PlaneActive_vsMultipl_WhenGeantEffiFail->Fill(numplaneactiveH0,totHitAllH0);
		  }
		  //END DEBUG THE PROCESS
		  }
	      
	      */
	      
		if(intmultclass<5){
		 H0_EtaPhiGeantTrkWhenEffiFailMultLess20->Fill(eta2geant_forEff,Phigeant_forEff);
		 //	 H0_EtaGeantTrkWhenEffiFailMultLess20->Fill(eta2geant_forEff);
		}


	      if(NumGeantTrk_PassingCutsDiscrepancy>=1){
		  H0TrkEffi_AxisAvgMult_ZImpactCut->Fill(multclass,0);
	      }
	    }
	    
	  

	    /*
	    if((GQuarter1>0)||(GQuarter0>0)){
	      H0H1_PhiGeantTrkWhenEffiFail->Fill(Phigeant_forEff);
	      H0H1_EtaPhiGeantTrkWhenEffiFail->Fill(Phigeant_forEff,eta2geant_forEff);	    
	    }
	    */

	    if(GQuarter1>0){
	      H1_TrkEtaEfficiencyCumulative_UnfoldCuts->Fill(eta2geant_forEff,0);
	      H1_TrkPhiEfficiencyCumulative_UnfoldCutsCumul->Fill(Phigeant_forEff,0);
	      
	       multclass=(totHitAllH1/10.);
	      //class are in steps of 5 pad-clusters
	      intmultclass=(unsigned int)floor(multclass/5.);
	      if(intmultclass>9)//High multiplicity event goes in the same bin;
		intmultclass=9;	
	


	      /*
	      if(intmultclass<14)
	      {
		//	H0TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut->Fill(multclass,0);

		//BEGIN DEBUG THE PROCESS
		unsigned int numtrkinh0passingCuts=0;
		unsigned int numtrkinh0=0;
		double minDx=400.;double minDy=400.;double minDistance=400.;	      
		
		  {

		    double Zimpgeanteffifail=ZimpactFromRecTrack(TheGeantTrkReconstr2);
		      if(verbosity>1)
			std::cout<<"!!!!!!!!!!!!!! Inefficiency in H1 , filling H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult:"<<intmultclass<<" GeantTrK ZImp "<<Zimpgeanteffifail<<"# PadCLuinH1:"<<numplaneactiveH1<<"  #Strip: "<<numplaneactiveStripH1<<" Event:"<<numevent<<std::endl;
		    unsigned int numtrkinh0=0; 
		    std::vector<int> OneHalf;OneHalf.push_back(0);
		    ZImpOfAGeantTrkWhenEffiFailH0->Fill(Zimpgeanteffifail);
		    for(T1T2TrackCollection::const_iterator recTrack = theT2Tracks2->begin(); recTrack!=theT2Tracks2->end();++recTrack){
		      T1T2Track debtrk=(*recTrack);

		      double zimprectrackeffifail=ZimpactFromRecTrack(debtrk);
		     if(verbosity>1)
		       std::cout<<"RecTrk with ZImp"<<zimprectrackeffifail<<std::endl;
		      if(RecoTrk_UnfoldingCut(debtrk))
			if(T2CutsUtil.TrkAlsoInQuarter(debtrk,OneHalf)){		  
			  numtrkinh0++;
			  ZImpOfARecTrkWhenEffiFail->Fill(zimprectrackeffifail);
			  if(RecoTrk_UnfoldingCut(debtrk))
			    numtrkinh0passingCuts++;
			  
			  double thisdx=(TheGeantTrkReconstr2.GetHitT2(0).GetHitX()-recTrack->GetHitT2(0).GetHitX());
			  double thisdy=(TheGeantTrkReconstr2.GetHitT2(0).GetHitY()-recTrack->GetHitT2(0).GetHitY());
			  double Distance=sqrt(thisdx*thisdx+thisdy*thisdy);
			  if(Distance<minDistance){
			    minDx=thisdx;
			    minDy=thisdy;
			    minDistance=Distance;
			    
			  }
			}
		    }
		    Numtrk_vs_NumtrkpasCutsWhenEffiFail->Fill(numtrkinh0,numtrkinh0passingCuts);
		    
		    minDxDyWhenGeantEffiFail->Fill(minDx,minDy);//minPhi  minR
		    NumH0TrkWhenGeantEffiFail->Fill(numtrkinh0);		  
		    NumH0PlaneActiveWhenGeantEffiFail->Fill(numplaneactiveH0);
		    NumH0PlaneActive_vsMultipl_WhenGeantEffiFail->Fill(numplaneactiveH0,totHitAllH0);
		  }
		  //END DEBUG THE PROCESS
		  }

	      */



	      H1_BinMultStat->Fill(intmultclass);
	      H1_BinMultStatFail->Fill(intmultclass);
	      H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(eta2geant_forEff,0);
	      H1_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(Phigeant_forEff,0);
	    }
	  }
	  
	  if(Garm1>0){
	    ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts->Fill(eta2geant_forEff,0);
	    if(GQuarter2>0){
	      H2_TrkEtaEfficiencyCumulative_UnfoldCuts->Fill(eta2geant_forEff,0);
	       multclass=(totHitAllH2/10.);
	      //class are in steps of 5 pad-clusters
	       intmultclass=(unsigned int)floor(multclass/5.);
	      if(intmultclass>9)//High multiplicity event goes in the same bin;
		intmultclass=9;		

	      /*
	      {
		//	H0TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut->Fill(multclass,0);

		//BEGIN DEBUG THE PROCESS
		unsigned int numtrkinh0passingCuts=0;
		unsigned int numtrkinh0=0;
		double minDx=400.;double minDy=400.;double minDistance=400.;	      
		
		  {

		    double Zimpgeanteffifail=ZimpactFromRecTrack(TheGeantTrkReconstr2);
		      if(verbosity>1)
			std::cout<<"!!!!!!!!!!!!!! Inefficiency H3, filling H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult:"<<intmultclass<<" GeantTrK ZImp "<<Zimpgeanteffifail<< "#Strip: "<<numplaneactiveStripH0<<std::endl;
		    unsigned int numtrkinh0=0; 
		    std::vector<int> OneHalf;OneHalf.push_back(0);
		    ZImpOfAGeantTrkWhenEffiFailH0->Fill(Zimpgeanteffifail);
		    for(T1T2TrackCollection::const_iterator recTrack = theT2Tracks2->begin(); recTrack!=theT2Tracks2->end();++recTrack){
		      T1T2Track debtrk=(*recTrack);

		      double zimprectrackeffifail=ZimpactFromRecTrack(debtrk);
		     if(verbosity>1)
		       std::cout<<"RecTrk with ZImp"<<zimprectrackeffifail<<std::endl;
		      if(RecoTrk_UnfoldingCut(debtrk))
			if(T2CutsUtil.TrkAlsoInQuarter(debtrk,OneHalf)){		  
			  numtrkinh0++;
			  ZImpOfARecTrkWhenEffiFail->Fill(zimprectrackeffifail);
			  if(RecoTrk_UnfoldingCut(debtrk))
			    numtrkinh0passingCuts++;
			  
			  double thisdx=(TheGeantTrkReconstr2.GetHitT2(0).GetHitX()-recTrack->GetHitT2(0).GetHitX());
			  double thisdy=(TheGeantTrkReconstr2.GetHitT2(0).GetHitY()-recTrack->GetHitT2(0).GetHitY());
			  double Distance=sqrt(thisdx*thisdx+thisdy*thisdy);
			  if(Distance<minDistance){
			    minDx=thisdx;
			    minDy=thisdy;
			    minDistance=Distance;
			    
			  }
			}
		    }
		    Numtrk_vs_NumtrkpasCutsWhenEffiFail->Fill(numtrkinh0,numtrkinh0passingCuts);
		    
		    minDxDyWhenGeantEffiFail->Fill(minDx,minDy);//minPhi  minR
		    NumH0TrkWhenGeantEffiFail->Fill(numtrkinh0);		  
		    NumH0PlaneActiveWhenGeantEffiFail->Fill(numplaneactiveH0);
		    NumH0PlaneActive_vsMultipl_WhenGeantEffiFail->Fill(numplaneactiveH0,totHitAllH0);
		  }
		  //END DEBUG THE PROCESS
		  }
	      */

	      H2_BinMultStat->Fill(intmultclass);
	      H2_BinMultStatFail->Fill(intmultclass);
	      H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(eta2geant_forEff,0);
	      if((Phigeant_forEff>=37.)&&(Phigeant_forEff<=329.))
		H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_PhiCut[intmultclass]->Fill(eta2geant_forEff,0);
	    }
	   if(GQuarter3>0){
	      H3_TrkEtaEfficiencyCumulative_UnfoldCuts->Fill(eta2geant_forEff,0);	      
	       multclass=(totHitAllH3/10.);
	      //class are in steps of 5 pad-clusters
	      intmultclass=(unsigned int)floor(multclass/5.);
	      if(intmultclass>9)//High multiplicity event goes in the same bin;
		intmultclass=9;	
	

	 

	      H3_BinMultStat->Fill(intmultclass);
	      H3_BinMultStatFail->Fill(intmultclass);
	      H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[intmultclass]->Fill(eta2geant_forEff,0);
	      
	    }
	  }


	}
	
	}
	//std::cout<<"Here4"<<std::endl;
  	
      }//if(NumT2HitinGeantTrkPrimary>=4) 



    }//if(trkprimcond)
    if(verbosity>1)
      std::cout<<"Next geant trk:"<<std::endl;


  
 }// END of the SIMULATION TRACK Loop

 if(verbosity>1)
   std::cout<<"Here5: end main SimTrk loop"<<std::endl;

//std::cout<<"G"<<std::endl;
// if(ghostfound==false)
  
   

 // if(wantedconditionfordndeta)

 for(unsigned int bin=0; bin<vectorEta2CounterFakePrimary.size(); bin++)
   {
     double etacentre=-maxetafordndetahisto+bin*etabinsize+ etabinsize/2.0;

     FakeSecondaryRecoDnDeta2FromVtx->Fill(etacentre,vectorEta2CounterFakePrimary.at(bin)*(1.0/etabinsize)); 

     PrimaryRecoDnDeta2FromVtx->Fill(etacentre,vectorEta2CounterDnDetaPrimary.at(bin)*(1.0/etabinsize));

     PrimaryRecoDnDeta2NoCut->Fill(etacentre,vectorEta2CounterDnDetaPrimaryNoCut.at(bin)*(1.0/etabinsize));

     RecoDnDeta2FromVtx->Fill(etacentre,vectorEta2CounterDnDeta.at(bin)*(1.0/etabinsize)); 

     //     if(verbosity>1)
     //std::cout<<"In counter"<<std::endl;
     
     //Put on 13-8-2011
     PrimaryGeantDNDeta2->Fill(etacentre,vectorEta2SimulatedPrimaryTrkInT2_H0.at(bin)*(1.0/etabinsize));
     
     PrimaryGeantDNDeta2H0_ZCut->Fill(etacentre,vectorEta2GeantPrimaryTrkInH0_ZCut.at(bin)*(1.0/etabinsize));
     
     
      if(AtLeastOneTrkZimpInH0){
	PrimaryGeantDNDeta2_IfOneRecoZImpInH0->Fill(etacentre,vectorEta2SimulatedPrimaryTrkInT2_H0.at(bin)*(1.0/etabinsize));	
      }
      //if(verbosity>1)
      // std::cout<<"In countera"<<std::endl;
      if(t2trackVectorBothSide.size()>0){
	PrimaryGeantDNDeta2_IfOneTrkInT2_H1->Fill(etacentre,vectorEta2SimulatedPrimaryTrkInT2_H1.at(bin)*(1.0/etabinsize));
	PrimaryGeantDNDeta2_IfOneTrkInT2_H0->Fill(etacentre,vectorEta2SimulatedPrimaryTrkInT2_H0.at(bin)*(1.0/etabinsize));
	if(sometrkexcluding64)
	  PrimaryGeantDNDeta2_IfOneTrkInT2_H1sometrkexcluding64->Fill(etacentre,vectorEta2SimulatedPrimaryTrkInT2_H1.at(bin)*(1.0/etabinsize));
      }
      //if(verbosity>1)
      // std::cout<<"In counter1"<<std::endl;
     
      if(generator_Inelastic){	    
	//it uses EnergyCutinPrimaryEfficiency          
	DNDetaMBALLMCGenerator_GeneratorTriggered->Fill(etacentre,vectorMBALLMCGenerator.at(bin)*(1.0/etabinsize));	
      }

      // if(verbosity>1)
      //	std::cout<<"In counter2"<<std::endl;
      
      if(generator_Inelastic5364){
	DNDetaMBALLMCGenerator_GeneratorTriggeredAtLeastOne5364->Fill(etacentre,vectorMBALLMCGenerator.at(bin)*(1.0/etabinsize));
      }
      
      // if(verbosity>1)
      //	std::cout<<"In counter3"<<std::endl;
      
      RecoDnDeta2NoCut->Fill(etacentre,vectorEta2CounterDnDetaNoCut.at(bin)*(1.0/etabinsize)); 
      
      RecoDnDeta2UnfoldCut2->Fill(etacentre,vectorEta2CounterRecoDnDeta2UnfoldCut2.at(bin)*(1.0/etabinsize)); 
      
      RecoDnDeta2UnfoldCut2_Primary->Fill(etacentre,vectorEta2CounterRecoDnDeta2UnfoldCut2_Primary.at(bin)*(1.0/etabinsize)); 
      //if(verbosity>1)
      //std::cout<<"In counter4"<<std::endl;
      RecoDnDeta2UnfoldCut->Fill(etacentre,vectorEta2CounterRecoDnDeta2UnfoldCut.at(bin)*(1.0/etabinsize)); 
      
      RecoDnDeta2UnfoldCut_Primary->Fill(etacentre,vectorEta2CounterRecoDnDeta2UnfoldCut_Primary.at(bin)*(1.0/etabinsize)); 

   }

 //if(verbosity>1)
 //  std::cout<<"End analyze"<<std::endl;

  //CLASSIFY RECONSTRUCTED TRACKS
  //MatchRecotracksToSimulatedTracks(iEvent,theSimHits,theSimTracks,theSimVertices,trackHitList,theT2Tracks2);

 // return;
}









//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//     ------------------      PRODUCE  END ---------------------     
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////




void T2BackgroundAn::TrkDeltaTheta(T1T2Track &atrk,double &DeltaThetaX,double &DeltaThetaY){

  double trk_Tx=atrk.GetTx();
  double trk_Ty=atrk.GetTy();
  double hitx=0;double hity=0;double hitz=0;
  double hit_Tx=0.;
  double hit_Ty=0.;
  
  for(unsigned int u=0;u<atrk.GetHitEntries();u++){
    if(atrk.GetHitT2(u).GetHitClass()==1)
      if((atrk.GetHitT2(u).GetHitNumPad()<5)&&(atrk.GetHitT2(u).GetHitNumStrip()<5))
	{
	  hitz=atrk.GetHitT2(u).GetHitZ();
	  hitx=atrk.GetHitT2(u).GetHitX();
	  hity=atrk.GetHitT2(u).GetHitY();	  	    
	  continue;
	}
  }
  
  if(hitz!=0){
    hit_Tx=hitx/hitz;
    hit_Ty=hity/hitz;

    DeltaThetaX=trk_Tx-hit_Tx;
    DeltaThetaY=trk_Ty-hit_Ty;
    
  }
  
  
}
	  

bool T2BackgroundAn::RecoTrk_UnfoldingCut_2(T1T2Track &atrk){

  bool cutpassed=false;
  unsigned int counter=0;
  double DeltaThetaX=0.;
  double DeltaThetaY=0.;
  TrkDeltaTheta(atrk,DeltaThetaX,DeltaThetaY);
  for(unsigned int u=0;u<atrk.GetHitEntries();u++)
    if(atrk.GetHitT2(u).GetHitClass()==1)
      counter++;



  if((counter>=3)&&(atrk.GetHitEntries()>=4))
    if(T2CutsUtil.ChiCutCond(atrk, true, 0.01, 0.01))
      if(fabs(DeltaThetaX)<0.035)
	if(fabs(DeltaThetaY)<0.035)
	  cutpassed=true;
	     
	     
  return cutpassed;
}
          

bool T2BackgroundAn::RecoTrk_UnfoldingCut_56Division(T1T2Track &atrk){
bool cutpassed=false;
  unsigned int counter=0;
  //std::cout<<"Caso0bb RecoTrk_UnfoldingCut"<<atrk.GetHitEntries()<<std::endl;
  for(unsigned int u=0;u<atrk.GetHitEntries();u++)
    if(atrk.GetHitT2(u).GetHitClass()==1)
      counter++;
  //std::cout<<"Caso0cc RecoTrk_UnfoldingCut"<<atrk.GetHitEntries()<<std::endl;

  double C0=(atrk.GetHitT2(0).GetHitX()*atrk.GetHitT2(0).GetHitX()+atrk.GetHitT2(0).GetHitY()*atrk.GetHitT2(0).GetHitY())/(atrk.GetTx()*atrk.GetHitT2(0).GetHitX()+atrk.GetTy()*atrk.GetHitT2(0).GetHitY());
  double Z0impact=atrk.GetHitT2(0).GetHitZ()-C0;

  double eta2=T2CutsUtil.EtaFromAveragesVtxCorr(atrk,0.,0.,0.);
  if((fabs(atrk.Z_at_Rmin())<5000)||((atrk.Z_at_Rmin())*atrk.GetHitT2(0).GetHitZ()<0)){  
    if((fabs(Z0impact)<5000)&&(eta2>5.6))
      cutpassed=true;

    if((fabs(Z0impact)<10000)&&(eta2<5.6))
      cutpassed=true;
    
  }
  return cutpassed;
}



bool T2BackgroundAn::RecoTrk_UnfoldingCut(T1T2Track &atrk){

  
  //std::cout<<"Caso0aa RecoTrk_UnfoldingCut"<<std::endl; 
  bool cutpassed=false;
  unsigned int countercl1Hit=0;
  //std::cout<<"Caso0bb RecoTrk_UnfoldingCut"<<atrk.GetHitEntries()<<std::endl;
  for(unsigned int u=0;u<atrk.GetHitEntries();u++)
    if(atrk.GetHitT2(u).GetHitClass()==1)
      countercl1Hit++;
  //std::cout<<"Caso0cc RecoTrk_UnfoldingCut"<<atrk.GetHitEntries()<<std::endl;
  
  double C0=(atrk.GetHitT2(0).GetHitX()*atrk.GetHitT2(0).GetHitX()+atrk.GetHitT2(0).GetHitY()*atrk.GetHitT2(0).GetHitY())/(atrk.GetTx()*atrk.GetHitT2(0).GetHitX()+atrk.GetTy()*atrk.GetHitT2(0).GetHitY());
  double Z0impact=atrk.GetHitT2(0).GetHitZ()-C0;


  int symb=RawtoSymb(atrk.GetHitT2(0).GetHitDetRawId());
  symb=symb/10;
  double ZImp_Left=0.;double ZImp_Right=0.;
  
  
  double the_eta2=T2CutsUtil.EtaFromAveragesVtxCorr(atrk,0.,0.,0.);
  ZImpactRangeFromFit(symb,the_eta2,ZImp_Left, ZImp_Right);
  
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW Associated_PrimaryTrkZImpact_H0
  //ZImpact studies according to eta2

  if(symb==0){
    CutZimpLeftH0->Fill(the_eta2,ZImp_Left);
    CutZimpRightH0->Fill(the_eta2,ZImp_Right);
  }
   
  if(symb==1){
    CutZimpLeftH1->Fill(the_eta2,ZImp_Left);  
    CutZimpRightH1->Fill(the_eta2,ZImp_Right);
  }

  if(symb==2){
    CutZimpLeftH2->Fill(fabs(the_eta2),ZImp_Left);
    CutZimpRightH2->Fill(fabs(the_eta2),ZImp_Right);
  }

  if(symb==3){
    CutZimpLeftH3->Fill(fabs(the_eta2),ZImp_Left);
    CutZimpRightH3->Fill(fabs(the_eta2),ZImp_Right);
  }



  
  //7TeV case
  /*
    double ZimpCutThr=5000;
    ZimpCutThr=ZEffiCutImpact.at(symb);//7 TeV analysis
  if((fabs(atrk.Z_at_Rmin())<13500)||((atrk.Z_at_Rmin())*atrk.GetHitT2(0).GetHitZ()<0))   
    if(fabs(Z0impact)<ZimpCutThr)
      cutpassed=true;
  */



 //8TeV case

double ProbChi2_XY=TMath::Prob(atrk.ChiSquared(),(countercl1Hit*2-4));
if(ProbChi2_XY>0.01)
  if((fabs(atrk.Z_at_Rmin())<13500)||((atrk.Z_at_Rmin())*atrk.GetHitT2(0).GetHitZ()<0))   
    if(Z0impact<ZImp_Right)
      if(Z0impact>ZImp_Left)
	cutpassed=true;


/*
 //8TeV case-LargeWidth
double ProbChi2_XY=TMath::Prob(atrk.ChiSquared(),(countercl1Hit*2-4));
if(ProbChi2_XY>0.01)
  if((fabs(atrk.Z_at_Rmin())<13500)||((atrk.Z_at_Rmin())*atrk.GetHitT2(0).GetHitZ()<0))   
    if(Z0impact<(ZImp_Right*1.3))
      if(Z0impact>(ZImp_Left*1.3))
	cutpassed=true;
*/

/*
 //8TeV case-LargeWidth-Nochi2
 if((fabs(atrk.Z_at_Rmin())<13500)||((atrk.Z_at_Rmin())*atrk.GetHitT2(0).GetHitZ()<0))   
    if(Z0impact<(ZImp_Right*1.3))
      if(Z0impact>(ZImp_Left*1.3))
	cutpassed=true;
*/



  /*
  if((countercl1Hit>=3)&&(atrk.GetHitEntries()>=4))
    if(T2CutsUtil.ChiCutCond(atrk, true, 0.01, 0.01))
       if(sqrt(atrk.bx_*atrk.bx_+atrk.by_*atrk.by_)<requiredVtxPositionDNDeta)
	 if(fabs(atrk.Eta())>4.5)
	   if((atrk.Eta()/fabs(atrk.Eta()))*atrk.Z_at_Rmin()<11000)
	     cutpassed=true;
  */	     
	     
  return cutpassed;
}







bool T2BackgroundAn::AreTrkFromIP_PyChParticle(int Hit_TrkId,const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, const HepMC::GenEvent* evt)
{
  
  bool out1=false;
  bool out2=false;
  SimTrack PythiaIPChParticleTrack;
  SimVertex simuvertex;

  if(getTrack(Hit_TrkId,theSimTracks,PythiaIPChParticleTrack) == true)
    {
      if(isTrackPrimary(Hit_TrkId, theSimTracks, theSimVertices))
	out1=true;
      else
	out1=false;
      
      out2=false;
      

      // Note that PythiaIPChParticleTrack.genpartIndex() is the barcode, not the index within GenParticleCollection, so I have to search the particle

      for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
	if((*p)->barcode()==PythiaIPChParticleTrack.genpartIndex())
	  {
	    (*p)->pdg_id();
	    if((*p)->status()==1)		 
	      out2=true;		
	  }
      }
	// std::cout<<"(Atrack.type()): "<<(PythiaIPChParticleTrack.type())<<std::endl; 
    }
  else
    {
     
      //	std::cout<<"AreTrkFromIP_PyChParticle: Hit_TrkId didn't associated to a Geant Trk: Geant fa cagare"<<std::endl;
    }

  bool outfinal=false;
  
  if(out2!=out1)
    std::cout<<"Incompatible results in AreTrkFromIP_PyChParticle"<<std::endl;
  else
    {
      outfinal=out1;
      
    }

  return outfinal;  
 
}






double T2BackgroundAn::MostProbableTrkEnergy(std::vector<double> &PrimaryTrkenergy_ForTheRecoHits){




  double theE2=0.;
  for (unsigned int h=0;h<PrimaryTrkenergy_ForTheRecoHits.size();h++)
    {
      theE2=theE2+PrimaryTrkenergy_ForTheRecoHits.at(h);
    }
  if(theE2>0.)
    theE2=theE2/PrimaryTrkenergy_ForTheRecoHits.size();
  else
    theE2=-1.;


  std::vector<std::vector<double> > SameEnergyHits;

  for(unsigned int y=0; y<PrimaryTrkenergy_ForTheRecoHits.size();y++){

    bool asameenergyHitfound=false;

    if(y==0){
      std::vector<double> EnergyV;
      EnergyV.push_back(PrimaryTrkenergy_ForTheRecoHits.at(y));
      SameEnergyHits.push_back(EnergyV);
    }
    else{
      
      for(unsigned int i=0;i<SameEnergyHits.size();i++)
	for(unsigned int j=0;j<SameEnergyHits.at(i).size();j++)
	  if(fabs(PrimaryTrkenergy_ForTheRecoHits.at(y)-SameEnergyHits.at(i).at(j))<1.0){
	    
	    SameEnergyHits.at(i).push_back(PrimaryTrkenergy_ForTheRecoHits.at(y));
	    asameenergyHitfound=true;
	    goto nexthit;
	  }
    }
    
   
   

  nexthit: // goto label
   
    if(asameenergyHitfound==false)
      {
	std::vector<double> EnergyV;
	EnergyV.push_back(PrimaryTrkenergy_ForTheRecoHits.at(y)); 
	SameEnergyHits.push_back(EnergyV);
      }
    
    
  }

  
for(unsigned int i=0;i<SameEnergyHits.size();i++)
  if(SameEnergyHits.at(i).size()>(0.7*PrimaryTrkenergy_ForTheRecoHits.size()))
    SameEnergyHits.at(i).at(0);


//return theE;
  return theE2;
}




bool T2BackgroundAn::CUTfordNdeta(T1T2Track Trk)
{
  bool toret=false;
  unsigned int counter=0;
  
  for(unsigned int u=0;u<Trk.GetHitEntries();u++)
   if(Trk.GetHitT2(u).GetHitClass()!=9)
     counter++;
    // if(Trk.GetHitT2(u).GetHitClass()==1)
    //counter++;

  //if(Trk.GetHitT2(u).GetHitClass()!=9)
  // counter++;
  
  if(counter>=4)///(*TCiter)
    if(T2CutsUtil.ChiCutCond(Trk))
      toret=true;
  
  return toret;
}



bool T2BackgroundAn::CUTfordNdetaVtx(T1T2Track trk, bool fitxy)
{
  bool toret=false;
  double requiredVtxPositionDNDeta = 60.;
  if((trk.bx_*trk.bx_+trk.by_*trk.by_)<requiredVtxPositionDNDeta)
    if(fabs(trk.Z_at_Rmin())<5000)    
      toret=true;  
  
  return toret;
}


 bool T2BackgroundAn::IsAcceptedForDNDeta(T1T2Track Trk)
{
  bool toret=false;
  bool cutfordNdeta= CUTfordNdeta(Trk);
  std::vector<T2Hit> hitvector;
  for(unsigned int mm=0;mm<Trk.GetHitEntries();mm++)
    hitvector.push_back(Trk.GetHitT2(mm));	  
	 
  T1T2Track trk2rz=T2CutsUtil.TrackFromHits(false,hitvector);//RZFIT
  std::vector<int>  QuarterWanted;QuarterWanted=T2_QuarterUsed;
  
  if(fabs(Trk.Eta())<maxetafordndetahisto)
    { 
      if(CUTfordNdetaVtx(trk2rz,false))
	if(cutfordNdeta)
	if(T2CutsUtil.TrkAlsoInQuarter(trk2rz,QuarterWanted)) //IMPORTANT: Here You need T2QuarterUsed=0
	  if(fabs(trk2rz.Eta())<maxetafordndetahisto)	
	    toret=true;
    }
  
  return toret;
}

bool T2BackgroundAn::VertexFromIonPump(double Zvtxposition)
{
  bool toret=false;
  //in cm
  if((fabs(Zvtxposition)>1300.)&&(fabs(Zvtxposition)<1370.))
    {
      toret=true;
    }
  
  return toret;
}


double T2BackgroundAn::GetEta2FromGeantHits(std::vector<std::pair<double,double> > SimHitRZ_ofGeantPrimaryTrk)
{
  
  double thetaHitAverage=0.;     
  double hemis=0.; 
  for(unsigned int u=0;u<SimHitRZ_ofGeantPrimaryTrk.size();u++)
    {
      std::pair<double, double> SimHitRZ_PrimTrk=SimHitRZ_ofGeantPrimaryTrk.at(u);
      thetaHitAverage += fabs(SimHitRZ_PrimTrk.first/SimHitRZ_PrimTrk.second);
      hemis=fabs(SimHitRZ_PrimTrk.second)/SimHitRZ_PrimTrk.second;
    }

  thetaHitAverage=thetaHitAverage/SimHitRZ_ofGeantPrimaryTrk.size();  
  double TrkEtaFromHitAverages_4Hits = (hemis*(-1.0)*log(thetaHitAverage/2.0)); //low angle approximation
  
  return TrkEtaFromHitAverages_4Hits;
    
}



bool T2BackgroundAn::ThisParticleHasStableDaughters(HepMC::GenEvent::particle_const_iterator p)
{

  bool daughtStable=false;
  for (HepMC::GenVertex::particle_iterator des =(*p)->end_vertex()->particles_begin(HepMC::children);
       des != (*p)->end_vertex()->particles_end(HepMC::children); ++des ) 
    {
      if((*des)->status()==1)
	daughtStable=true;
    }

  return daughtStable;
}



void T2BackgroundAn::PrintParents(HepMC::GenEvent::particle_const_iterator p)
{

  std::cout<<"Pythia Ancestors Particle-ID : ";
  for (HepMC::GenVertex::particle_iterator des =(*p)->end_vertex()->particles_begin(HepMC::ancestors);
       des != (*p)->end_vertex()->particles_end(HepMC::ancestors); ++des ) 
    {
     
      std::cout<<"  "<<(*des)->pdg_id();
    }
  std::cout<<std::endl;

  //return daughtStable;
}



int T2BackgroundAn::GetTrkOlderMotherPid(SimTrack aTrack,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
						const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,const HepMC::GenEvent* evt, int &motherbarcode,int &thePidGeneratorOfThisTrk)
{
  
  /////////////////////////////////
  int OlderMotherParticle_Id=0; //This will be returned.
  /////////////////////////////////

  unsigned int count=0;

  SimTrack DirectMotherTrk=aTrack;
  int DirectMotherTrk_Id=aTrack.trackId();
  
 
  for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
    if((*p)->barcode()==aTrack.genpartIndex())
      {
	// thePidGenerator=(*p)->pdg_id();
	((*p)->barcode());
	motherbarcode=((*p)->barcode());
	thePidGeneratorOfThisTrk=(*p)->pdg_id();	
	OlderMotherParticle_Id=(*p)->pdg_id();		
      }
  }
 //--  std::cout<<"GetTrkOlderMotherPid begin. Trk Id: "<<DirectMotherTrk_Id<<" Pythia pdg:"<<OlderMotherParticle_Id<<" Trk_PID:"<<(DirectMotherTrk.type())<<std::endl;

 //std::cout<<"Try to find mother for Particle_ID: "<<(aTrack.type())<<" with Trk_ID:"<<DirectMotherTrk_Id<<std::endl;

 while(DirectMotherTrk_Id!=(-1))//loop until you dont' find the stable mother or the source
   {
     count++;     

     SimVertex simuvertex;
 
     //getTrack: Returns track corresponding to DirectMotherTrk_Id. Returns true if found, false otherwise.
  
     if(getTrack(DirectMotherTrk_Id,theSimTracks,DirectMotherTrk) == true)
       {
	  //-- std::cout<<"Retrieved geant Trk with Particle_PID: "<<(DirectMotherTrk.type())<<" and Trk_ID:"<<DirectMotherTrk_Id<<std::endl;

	 OlderMotherParticle_Id=DirectMotherTrk.type();
	 //std::cout<<"Mother Trk id= "<<DirectMotherTrk_Id<<std::endl;
	 
	 //This look only if parentVtx is primary (Index<0) 
	 if(isTrackPrimary(DirectMotherTrk_Id, theSimTracks, theSimVertices)==false)
	   {	    
	     //-- std::cout<<"Not a primary track. Try to find originating vertex..."<<std::endl;
	     //Take the vtx generating the trk.
	     if(getVertex(DirectMotherTrk.vertIndex(), theSimVertices,simuvertex) == true) 
	       {	    
		 
		//--  std::cout<<"Originating vertex found. Try to found the ancestors of this track:"<<std::endl;
		 // Note that DirectMotherTrk.genpartIndex() is the barcode, 
		 // not the index within GenParticleCollection, so I have to search the particle
		  

		 //THIS IS NEVER WORKING ACTUALLY

		  for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
		    if((*p)->barcode()==DirectMotherTrk.genpartIndex())
		      {
			(*p)->pdg_id();
		   
			//-- std::cout<<"Pythia Ancestors Particle-IDs of the track: ";
			HepMC::GenVertex::particle_iterator ancestor;
			ancestor = (*p)->production_vertex()->particles_begin(HepMC::parents);
			for(;ancestor != (*p)->production_vertex()->particles_end(HepMC::parents); ++ancestor) {
			 
			  //std::cout<<"  "<<((*ancestor)->pdg_id());
			  //motherbarcode to return is saved.
			  motherbarcode=(*ancestor)->barcode();
			  
			  if((std::find(knowPart.begin(),knowPart.end(),fabs((*ancestor)->pdg_id()))!=knowPart.end()))
			    {
			      ((*ancestor)->barcode());
			      OlderMotherParticle_Id=(*ancestor)->pdg_id();
			      //--std::cout<<"Ancestor is a Known particle: Barcode:"<<OlderMotherId<<" PID:"<<OlderMotherParticle_Id<<std::endl;
			    }

			}
			//-- std::cout<<"OlderMotherParticle_Id is now: "<<	OlderMotherParticle_Id<<std::endl;		   		   		   
		      }
		    //std::cout<<"Begin-2"<<motherbarcode<<std::endl;
		  }
		 
		  //Find the mother track originating this vtx 
		  //This instruction allows to go back in the three to found the origin
		   //THIS IS THE ONLY IMPORTANT INSTRUCTION
		  DirectMotherTrk_Id=simuvertex.parentIndex(); 		  
		  //--std::cout<<"Mother Trk Id originating this vtx:"<<DirectMotherTrk_Id<<std::endl;
		
	       }
	    //--  else
	        //-- std::cout<<"Originating vertex NOT found:"<<std::endl;
	   }
	 else //Track is primary
	   { 
	    //--  std::cout<<"isTrackPrimary==true"<<std::endl;
	      //In this case DirectMotherTrk_Id should be <0 .
	     // std::cout<<"isTrackPrimary==true so Mother Trk is < 0, I want to store this info also. Track printing:"<<std::endl;
	     // Note that DirectMotherTrk.genpartIndex() is the barcode, not the index within GenParticleCollection, so I have to search the particle


	     for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
	       if((*p)->barcode()==DirectMotherTrk.genpartIndex())
		 {
		  (*p)->pdg_id();
		   
		   //   std::cout<<"Pythia Ancestors Particle-ID : ";
		   HepMC::GenVertex::particle_iterator ancestor;
		   ancestor = (*p)->production_vertex()->particles_begin(HepMC::parents);
		   for(;ancestor != (*p)->production_vertex()->particles_end(HepMC::parents); ++ancestor) {
		     
		     motherbarcode=(*ancestor)->barcode();
		     if((std::find(knowPart.begin(),knowPart.end(),fabs((*ancestor)->pdg_id()))!=knowPart.end()))
		       {
			 ((*ancestor)->barcode());
			 OlderMotherParticle_Id=(*ancestor)->pdg_id();
			 motherbarcode=(*ancestor)->barcode();
		       }
		     //otherwise keep what you already have, since pythia is not giving any true particle
		   }
		   		   		   		   
		 }
	       //std::cout<<"Begin-2"<<motherbarcode<<std::endl;
	     }
	     
	     //std::cout<<"Exit condition active since Mother Trk was primary"<<std::endl;
	     DirectMotherTrk_Id=-1;//Exit condition active
	   }
       }
     else
      DirectMotherTrk_Id=-1;//Exit condition active  but you have not recovered the trk  
     
   }//while end

 if(count==0) //The input trkId was already primary
   {


      for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
	if((*p)->barcode()==aTrack.genpartIndex())
	  {
	    (*p)->pdg_id();
	   ((*p)->barcode());
	    motherbarcode=((*p)->barcode());
	    OlderMotherParticle_Id=(*p)->pdg_id();		
	  }
      }

   }

 // std::cout<<"GetTrkOlderMotherPid Return Trk_Id:"<<OlderMotherId<<std::endl;
 return OlderMotherParticle_Id;

}





bool T2BackgroundAn::NeutralDecay(int pid){
 bool  toreturn=false;
  
  if(fabs(pid)==111)
    toreturn=true;
  
  
  return toreturn;
}


bool T2BackgroundAn::PrimaryTrackCondition(SimTrack atrk,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt,int &barcodeMother){


  bool primarycondition=false;
  
  int trackId=atrk.trackId();
  //SimTrack is stored and the generating vertex is at IP.
  
  bool isTrackPrimaryIP=false;
  if(fabs(atrk.charge())>0){
    isTrackPrimaryIP = isTrackPrimary(trackId, theSimTracks, theSimVertices);
  }
    
    //SimTrack PythiaIPChParticleTrack=(*ittrack);      
  int thePidGenerator=0; int thePidGeneratorOfThisTrk=0;
  
  thePidGenerator=GetTrkOlderMotherPid(atrk,theSimVertices,theSimTracks,evt,barcodeMother,thePidGeneratorOfThisTrk);

  if(thePidGeneratorOfThisTrk==0) //It means that Geant trk was not associated to a pythia particle
    thePidGeneratorOfThisTrk=atrk.type();

  //bool isTrackPrimaryIP = isTrackPrimary(trackId, theSimTracks, theSimVertices);     	        

  //if(verbosity>0)
  bool alloweddecayparticle=false;
  if(std::find(AllowedDecayToCount.begin(),AllowedDecayToCount.end(),fabs(thePidGenerator))!=AllowedDecayToCount.end())
    alloweddecayparticle=true;

  //I'm Excluding 211 also because I assume that a nonprimary older mother 211 have done Nuclear interaction
  if(fabs(atrk.charge())>0)
    {
      //&&(fabs(thePidGeneratorOfThisTrk)!=11)&&(fabs(thePidGeneratorOfThisTrk)!=22)
      if((isTrackPrimaryIP)||((isTrackPrimaryIP==false)&&(fabs(thePidGenerator)!=211)&&(thePidGeneratorOfThisTrk!=0)&&(alloweddecayparticle==true)))
	primarycondition=true;
    }

  //   std::cout<<"GeantTrkID:"<<trackId<<" PID_Mother:"<<thePidGenerator<<" Track PID: "<<thePidGeneratorOfThisTrk<<" Direct Primary_Condition: "<<isTrackPrimaryIP<<" Primary_Condition: "<<primarycondition<<std::endl;
  //Just a patch.
  primarycondition=false;
  if(isTrackPrimaryIP)
    primarycondition=true;
  return primarycondition;

}





//Take a RECO track and return an information about is nature: primary, secondary,none
//Secondary are also Associated to reco Hit which Matches with Simhit not pointing to any stored SimTrk.
//This is the returned vector for a RecoTrk given in input:
//trkinformation.push_back(numdifferentsimutrks);trkinformation.push_back(numprimary);
//trkinformation.push_back(numsecondary);trkinformation.push_back(nonmatchingHits);

std::vector<int>  T2BackgroundAn::RecotracksInfoFromSimu(const edm::Event& event,
							 const std::auto_ptr<edm::PSimHitContainer>& theSimHits,
							 const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
							 const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
							 const std::map<unsigned int, std::list<Local3DPoint> >& trackHitList,
							 const HepMC::GenEvent* evt,
							 T1T2Track recTrack,
							 std::vector<double> &PrimaryTrkenergy_ForTheRecoHits
							 )
{

  if (verbosity > 1){
    std::cout << "Matching simhits to recotracks.." << std::endl;
  }

  //Container for primary and secondary tracks
  auto_ptr<T1T2TrackCollection> theT2PrimaryTracks (new T1T2TrackCollection());
  auto_ptr<T1T2TrackCollection> theT2SecondaryTracks (new T1T2TrackCollection());

  //Container for simulated primary tracks found 
  std::map<unsigned int, bool> reconstructedTrackForSimulatedPrimaryTrackFound;
  std::map<unsigned int, bool> reconstructedTrackForSimulatedSecondaryTrackFound;

  //Container for all simulated positions
  //pair<row, column>, pair<trackId,isPrimary>
  std::vector<std::map<std::pair<int, int>, std::pair<unsigned int,bool> > > simPadHitPositionsInPlanes(40);
  
  //Convert simhits into information about plane, strip and pad.
  T2Geometry t2geom;
  T2GeometryUtil conv;     
  
  if (verbosity > 1){
    std::cout << "Analysing simhits.." << std::endl;
  }


  //Load simulated hits and associate for each corrisponding
  // T2-pad the  information about primary/secondary which generates it
  for(edm::PSimHitContainer::const_iterator simHitIt = theSimHits->begin(); simHitIt != theSimHits->end(); ++simHitIt){


    unsigned int detUnitId = simHitIt->detUnitId();    
    T2GeometryUtil::T2DetInfo detInfo = conv.GetT2Info(detUnitId);
    t2geom.setPlane(detUnitId);
    
    unsigned int trackIdFromHit = simHitIt->trackId();
    

    bool isHitPrimary=false; SimTrack Atrack;
   int motherbarcode=99;
    if(getTrack(trackIdFromHit,theSimTracks,Atrack) == true)
      {
	//ToDecomment
	
	isHitPrimary=PrimaryTrackCondition(Atrack,theSimVertices,theSimTracks,evt,motherbarcode);
      }

    // bool isHitPrimary = isTrackPrimary(trackIdFromHit, theSimTracks, theSimVertices);
    //SimTrack PythiaIPChParticleTrack=(*ittrack);
    //
    //thePidGenerator=GetTrkOlderMotherPid(PythiaIPChParticleTrack,theSimVertices,theSimTracks,evt,motherbarcode);
    //bool isHitPrimary = isTrackPrimary(trackIdFromHit, theSimTracks, theSimVertices);     	        

    //((isHitPrimary)||((isHitPrimary==false)&&(fabs(thePidGenerator)!=11)&&(fabs(thePidGenerator)!=22)&&(thePidGenerator!=0)))

    Local3DPoint hitPos = simHitIt->localPosition();    
    Local3DPoint corPos = CalculateCorrectXYZPos(Local3DPoint(hitPos.x(),hitPos.y(),hitPos.z()),detUnitId);
    
    int padRow = t2geom.getNearestPadRow(&corPos);
    int padCol = t2geom.getNearestPadCol(&corPos);
    
    
    reconstructedTrackForSimulatedSecondaryTrackFound[trackIdFromHit] = false;//inizializzazione?
    //Association DetId->Pad(Row,Col)->(IdTrack,IsfromPrimary)
    simPadHitPositionsInPlanes[detInfo.symb][std::pair<int, int>(padRow,padCol)] = std::pair<int,bool>(trackIdFromHit,isHitPrimary);    
    //For Each detector you have now a map that assign each SimHit-activated Pad row-col to a (TrkId,PrimaryNature)
  }


  
  if (verbosity > 1){    
    std::cout << "Analysing recotrack.." << std::endl;
  }

  
  
  int hitCount = recTrack.GetHitEntries();
  std::vector<int> allSimutrkId_foraRecoTrk;
  std::map<unsigned int, std::pair<int, bool> > recoTrackMatchInfo;
  std::list<int> notFoundCorrespondingSimulatedHit;
  std::map<unsigned int, std::pair<int, int> > mapTrkId_pairPrim_Second;
  


  for (int h = 0;h< hitCount;++h){
    T2Hit hit = recTrack.GetHitT2(h);
    unsigned int detUnitId = hit.GetHitDetRawId();
    T2GeometryUtil::T2DetInfo detInfo = conv.GetT2Info(detUnitId);
    t2geom.setPlane(detUnitId);
    Local3DPoint corPos(hit.GetHitX(),hit.GetHitY(),hit.GetHitZ());
    
    int padRow = t2geom.getNearestPadRow(&corPos);
    int padCol = t2geom.getNearestPadCol(&corPos);

    //Try to find simulated hit
    std::map<std::pair<int, int>, std::pair<unsigned int,bool> >::const_iterator it = simPadHitPositionsInPlanes[detInfo.symb].find(std::pair<int,int>(padRow,padCol));
    int hit_trkid=0;
    
    if (it != simPadHitPositionsInPlanes[detInfo.symb].end()){
      
      hit_trkid=it->second.first;
      allSimutrkId_foraRecoTrk.push_back(hit_trkid);
      
      if(it->second.second){ //This reco hit is from primary
	mapTrkId_pairPrim_Second[hit_trkid]=std::pair<int,int>(mapTrkId_pairPrim_Second[hit_trkid].first+1,mapTrkId_pairPrim_Second[hit_trkid].second);
	 SimTrack Atrack_;
		    
	 if(getTrack(hit_trkid,theSimTracks,Atrack_) == true)
	   {
	     PrimaryTrkenergy_ForTheRecoHits.push_back(Atrack_.momentum().e());
	   }

      }
      else
	mapTrkId_pairPrim_Second[hit_trkid]=std::pair<int,int>(mapTrkId_pairPrim_Second[hit_trkid].first,mapTrkId_pairPrim_Second[hit_trkid].second+1);
      
      //Found simulated hit                                                            /*[trkid]*/
      recoTrackMatchInfo[hit_trkid] = std::pair<int,bool>(recoTrackMatchInfo[hit_trkid].first+1,it->second.second);
      
    } else {
      //No simulated hit found for this Reco-Geometrical Hit
      
      //Assign the reco Hit to a fake hit_trkid=0.
      mapTrkId_pairPrim_Second[0]= std::pair<int,int>(recoTrackMatchInfo[0].first+1,recoTrackMatchInfo[0].second+1);
      recoTrackMatchInfo[0] = std::pair<int,bool>(recoTrackMatchInfo[0].first+1,false);
      notFoundCorrespondingSimulatedHit.push_back(h);
    }
  }



  //Now mapTrkId_pairPrim_Second contains a set of hit_trkid. Each  hit_trkid is 
  //associated to *the number* of primaryReco Hit and SecondaryReco Hit.
   
      unsigned int numdifferentsimutrks=0;unsigned int nonmatchingHits=0;
      unsigned int numprimary=0;unsigned int numsecondary=0;
      

      for(std::map<unsigned int, std::pair<int, int> >::const_iterator it =  mapTrkId_pairPrim_Second.begin(); it != mapTrkId_pairPrim_Second.end(); ++it){
	
	if(it->first!=0)
	  {
	    numprimary+=it->second.first;
	    numsecondary+=it->second.second;
	    numdifferentsimutrks++;
	  }
	else
	  {
	    nonmatchingHits+=it->second.first;
	  }
      }

      /*
      
      if (verbosity > -1){
	std::cout<<"---------------------------------------"<<std::endl;
	std::cout << "Reconstructed track is calculated from  " << hitCount << " hits." << std::endl;
	std::cout << "#Hit recognized as primary trkinformation[1]: " << numprimary <<std::endl; 
	std::cout << "#Hit recognized as secondary trkinformation[2]: " << numsecondary <<std::endl;
	std::cout << "# of different simu trks composing the reco trk trkinformation[0]:  " << numdifferentsimutrks <<std::endl;
	std::cout << "# Non-matching Hits trkinformation[3]:"<<nonmatchingHits<<std::endl;
	std::cout<<"---------------------------------------"<<std::endl<<std::endl;
      }
      */

     
      std::vector<int> trkinformation;
      trkinformation.push_back(numdifferentsimutrks);
      trkinformation.push_back(numprimary);
      trkinformation.push_back(numsecondary);
      trkinformation.push_back(nonmatchingHits);

      // std::cout<<"In generation: trkinfo[0]<<trkinfo[1]<<trkinfo[2]<<trkinfo[3]"<<trkinformation[0]<<trkinformation[1]<<trkinformation[2]<<trkinformation[3]<<std::endl;

      if(fastSimulation==false)
      if(numsecondary>numprimary)
	{
	  
	  SimVertex vertex;
	  //if(allSimutrkId_foraRecoTrk.size()==1) //Otherwise it is not clear what is happening
	    for(unsigned int h=0;h<allSimutrkId_foraRecoTrk.size();h++)
	      {
		if(h==0)
		  {
		    int ParenttrackId=allSimutrkId_foraRecoTrk.at(h);
		    SimTrack Atrack; SimTrack Mothertrack;  SimTrack OLDMothertrack;
		    
		    if(getTrack(ParenttrackId,theSimTracks,Atrack) == true)
		      {
		       PID_ofSecRecoTrack->Fill(Atrack.type());

		       
		       double thermin=recTrack.Rmin();
		       double thezmin=recTrack.Z_at_Rmin();
		    
		       if(recTrack.Eta()>0.)
			  if((thermin<600)&&(fabs(thezmin)<8000))
			    PID_ofSecRecoTrackAsaFakePrimary->Fill(Atrack.type());
		       

		       if(getVertex(Atrack.vertIndex(), theSimVertices,vertex) == true)   
			 {
			
			
			   if(recTrack.Eta()>0.)
			     if((thermin<600)&&(fabs(thezmin)<8000))
			       {
				 SimuVtxPositionForAFakeRecoPrimaryTrk->Fill(vertex.position().z(),vertex.position().y());//This is in cm
			      
				 //nt motherId = -1;
				 if(vertex.noParent()==false) { // there is a parent to this vertex
				   // geant id of the mother
				   // Now note that  motherGeantId is NOT an index, it's a track id
				   int motherGeantId =   vertex.parentIndex(); 
				   if(getTrack(motherGeantId,theSimTracks,Mothertrack) == true)
				     PID_DircetMotherForAFakeRecoPrimaryTrk->Fill(Mothertrack.type());

				   //GET BARCODE OF THE ORIGINATING PYTHIA TRACKS
				   int motherbarcode;int thePidGeneratorOfThisTrk=0;
				   int OldestmotherGeantPId = GetTrkOlderMotherPid(Atrack,theSimVertices,theSimTracks,evt,motherbarcode,thePidGeneratorOfThisTrk); 
				   PID_OldestMotherForAFakeRecoPrimaryTrk->Fill(OldestmotherGeantPId);
				   
				 }			     
			       }
			 }
		       
		      }



		  }
		
	      }
	}

      

  if (verbosity > 1){
    std::cout << "Done mathcing simhits to recotracks." << std::endl << std::endl;
  }
  
  return trkinformation;
}









T1T2Track T2BackgroundAn::RPhiFit(std::vector<T2Hit> &hitvec2)
{

  T1T2Track trk(2);
  trk.SetDet(2);

  //std::cout<<"T1T2Track BEGIN T2SelectionCutUtils::RPhiFit"<<std::endl;

  // std::vector<T2Hit> hitvec; //Only Real Hits hits
  std::vector<T2Hit> hitvecR;  //Real Hits and Vtx
  

  for(unsigned int jj =0; jj<hitvec2.size(); jj++)
   {         
	 hitvecR.push_back(hitvec2.at(jj));
      
   }
  //std::cout<<"Inside MyLinearfitCorr .. init"<<std::endl;

  
  //  unsigned int  sizeHitv=hitvec.size();
  unsigned int  sizeHitvR=hitvecR.size();
  if(sizeHitvR<2)
   {
     std::cout<<" T2SelectionCutUtils::RPhiFit problem: Track with less than 2 Cl1 hits!! Continue the fitting.."<<std::endl;
     
   }


  //std::vector<T2Hit> hitvecraw;
  
  //std::vector<T2Hit> doublehitz;
  //std::vector<T2Hit> hits0360;


 unsigned int num0=0;
 unsigned int num360=0;
 

 for(unsigned int jj =0; jj<sizeHitvR; jj++)
   { //std::cout<<"Hitphi "<<hitvec[jj].GetHitPhi()<<std::endl;
    // if((hitvecR[jj].GetHitPhi()>0.)&&(hitvecR[jj].GetHitPhi()<10.))
    if((hitvecR[jj].GetHitPhi()>0.)&&(hitvecR[jj].GetHitPhi()<80.))
      num0++;

    //    if((hitvecR[jj].GetHitPhi()>350.)&&(hitvecR[jj].GetHitPhi()<360.))
    if((hitvecR[jj].GetHitPhi()>280.)&&(hitvecR[jj].GetHitPhi()<360.))
      num360++;
  }


 if((num0>0) && (num360>0)) 
   {
    
     for(unsigned int jj =0; jj<sizeHitvR; jj++)
       {	
	 if((hitvecR[jj].GetHitPhi()>0.)&&(hitvecR[jj].GetHitPhi()<80.))
	   hitvecR[jj].SetHitPhi(hitvecR[jj].GetHitPhi()+360.);
       }
   }



 int hemisphere = (int)(hitvecR[0].GetHitZ()/fabs(hitvecR[0].GetHitZ()));

 std::vector<float> r;
 std::vector<float> phi;
 std::vector<float> z;
 std::vector<float> er;
 std::vector<float> ephi;
 std::vector<float> ez;


 float Sr=0;
 float Srz=0;
 float Szz_r=0;
 float Sz_r=0; 
 float S0_r=0; 
 float Swphi=0;
 float Sphiwphi=0;
 
 float cl1err=(0.12+0.05);
 
 //loop on cl1hits+vtx
  for(unsigned int jj =0; jj<sizeHitvR; jj++)
    {
      
      r.push_back(hitvecR[jj].GetHitR());
      
      z.push_back(hitvecR[jj].GetHitZ());    
    //  er.push_back(hitvecR[jj].GetHitDR());
      if(hitvecR[jj].GetHitClass()==1)
	 er.push_back(cl1err);
      else
	er.push_back(hitvecR[jj].GetHitDR());

      //	er.push_back(0.12+0.05);

      //ephi.push_back(hitvecR[jj].GetHitDPhi());
      
      ez.push_back(hitvecR[jj].GetHitDZ());
      
    
      
      Srz += r[jj]*z[jj]/er[jj]/er[jj];
      Szz_r += z[jj]*z[jj]/er[jj]/er[jj];
      Sz_r += z[jj]/er[jj]/er[jj];
      Sr += r[jj]/er[jj]/er[jj];
      S0_r += 1.0/er[jj]/er[jj];


      if(hitvecR.at(jj).GetHitClass()!=9){
	phi.push_back(hitvecR[jj].GetHitPhi());
	ephi.push_back(0.8);
	Swphi += 1.0/ephi.at(ephi.size()-1)/ephi.at(ephi.size()-1);
	Sphiwphi+=phi.at(phi.size()-1)/ephi.at(ephi.size()-1)/ephi.at(ephi.size()-1);
	//	std::cout<<phi.at(phi.size()-1)<<"  ->"<<Sphiwphi<<std::endl;
	  //Swphi += 1.0/ephi[jj]/ephi[jj];
	  //Sphiwphi+= phi[jj]/ephi[jj]/ephi[jj];
      }
    }


  double a_rz = (Srz*S0_r - Sz_r*Sr) / (Szz_r*S0_r - Sz_r*Sz_r);   // angular coefficient
  double b_rz = (Sr*Szz_r - Sz_r*Srz) / (Szz_r*S0_r - Sz_r*Sz_r);  // intercept   R=(a_rz)Z + b_rz
  double e_a_rz = sqrt( S0_r / (S0_r*Szz_r - Sz_r*Sz_r) );         
  double e_b_rz = sqrt( Szz_r / (S0_r*Szz_r - Sz_r*Sz_r) );
  double phim=Sphiwphi/Swphi;
  double e_phim= 1.0/sqrt(Swphi);
  
  //std::cout<<"->->---->"<<phim<<std::endl;
  double covab= - (Sz_r) / (Szz_r*S0_r - Sz_r*Sz_r);
  
  double chi2r = 0.0;
  double chi2p = 0.0;
  double chi2= 0.0;
  double normchi2red=0.0;


  unsigned int cl1count=0;
  
  for(unsigned int jj =0; jj<sizeHitvR; jj++)
    {
      if(hitvecR[jj].GetHitClass()==1)
	{
	  chi2r += (a_rz*hitvecR[jj].GetHitZ()+b_rz - hitvecR[jj].GetHitR())*(a_rz*hitvecR[jj].GetHitZ()+b_rz - hitvecR[jj].GetHitR())/cl1err/cl1err;
	  cl1count++;
	}
      else
	chi2r += (a_rz*hitvecR[jj].GetHitZ()+b_rz - hitvecR[jj].GetHitR())*(a_rz*hitvecR[jj].GetHitZ()+b_rz - hitvecR[jj].GetHitR())/hitvecR[jj].GetHitDR()/hitvecR[jj].GetHitDR();
      
    }
  

  
  for(unsigned int jjj =0; jjj<sizeHitvR; jjj++){          
    
    chi2p += (phim - phi[jjj])*(phim - phi[jjj])/ephi[jjj]/ephi[jjj];
  }
  
  
  
  if(cl1count>=2) 
    {      
      chi2=  (chi2p+chi2r)/2.0;
      normchi2red=(chi2p/((double)(cl1count)-1)+ chi2r/((double)(cl1count)-2))/2.0 ;
    }
  else
    {
      chi2=1.0;
      normchi2red=1.0;
      std::cout<<"Error: T2SelectionCutUtils RPhiFit without at least 2 hits"<<std::endl;
    }
  

  if(phim>360.0)                // This could happen when num0>0 AND num360>0
    phim=phim-360.0;
  
  
  phim=phim*3.14159265/180.0;
  e_phim=e_phim*3.14159265/180.0;
  
  //std::cout<<"rphifit: "<<phim*180/3.14159265<<std::endl;
  
  
  TVectorD vect(4);
  for(int oo=0; oo<4; oo++)
    vect[oo]=0.;
  
  
  vect[0]=phim;    // phi intercept
  vect[1]=0.;      // phi ang. coeff.
  vect[2]=b_rz;    // r intercept
  vect[3]=a_rz;    // r ang. coeff.
  

  
  TMatrixD mat(4,4);
  for(unsigned int oo=0; oo<4; oo++)
    for(unsigned int ooo=0;ooo<4; ooo++)
      mat[oo][ooo]=0.;
  

  mat[0][0]=e_phim*e_phim;
  mat[0][1]=0.;
  mat[0][2]=0.;
  mat[0][3]=0.;
  mat[1][0]=0.;
  mat[1][1]=0.;
  mat[1][2]=0.;
  mat[1][3]=0.;
  mat[2][0]=0.;
  mat[2][1]=0.;
  mat[2][2]=e_b_rz*e_b_rz;
  mat[2][3]=covab;
  mat[3][0]=0.;
  mat[3][1]=0.;
  mat[3][2]=covab;
  mat[3][3]=e_a_rz*e_a_rz;
  



//  double theeta = fabs(-log(tan(thetheta/2.))) * (double)hemisphere;

  //std::cout<<"calc phim "<<phim*180.0/3.14159265<<" "<<hitvec2[0].GetHitPhi()<<" "<<hitvec[0].GetHitPhi()<<std::endl;
  
  T1T2Track fittedtrack(vect,mat,a_rz, b_rz, phim, e_a_rz, e_b_rz, e_phim, chi2,chi2r,chi2p, normchi2red, hemisphere,2);   
  return fittedtrack; 
}




T1T2Track T2BackgroundAn::MyLinearfitCorrDEV(std::vector<T2Hit> hitvec2,TMatrixD &par_covariance_matrix,double &chi2_,bool Usestrip, int RoadID){

  std::vector<T2Hit> hitvec; hitvec.clear();
   
 
  double sigmaR=(0.12+0.05);
  double sigmaPhi=0.015;
  unsigned int numphimin20=0;
  unsigned int numphimag340=0;

  double phirad=0.;
  double r=0.;


  bool inserted=false;

  unsigned int numHit_t2=0;
  unsigned int numCl1HitHit_t2=0;
  unsigned int numStripOnlyHit_t2=0;
  unsigned int numPadOnlyHit_t2=0;
 
  if(verbosity){
    std::cout<<" MyLinearfitCorrDEV start with "<<hitvec2.size()<<std::endl;
    
  }

  //std::cout<<" DEBTRK MyLinearfitCorrDEV start with "<<hitvec2.size()<<std::endl;

  /////////////////////////////////////////////////////////////////////////////////
  ///           SAVE IN HITVEC ONLY WANTED HITS
  ///////////////////////////////////////////////////////////////////////////////

  for(unsigned int m=0;m<hitvec2.size();m++)
    { 

      if(hitvec2.at(m).GetHitR()<1)
	std::cout<<"-------> INPUT tracker Error: R hit too small:"<<std::endl;

      inserted=false;
      if(hitvec2.at(m).GetHitNumPad()>0)
	{
	  if((hitvec2.at(m).GetHitPhi()<85)&&(hitvec2.at(m).GetHitPhi()>0))//was <90
	    numphimin20++;
	
	  if((hitvec2.at(m).GetHitPhi()>280))//was>340
	    numphimag340++; 
	}
   
      
      if(hitvec2.at(m).GetHitNumPad()==0){
	if(Usestrip==true){
	  inserted=true ;
	  hitvec.push_back(hitvec2.at(m));
	}
      }
      else
	{
	  inserted=true ;
	  hitvec.push_back(hitvec2.at(m));
	}


      if(inserted){

	if(verbosity)
	  std::cout<<"Z-Phi: "<<hitvec2.at(m).GetHitZ()<<" "<<hitvec2.at(m).GetHitPhi()<<" Num Pad-Strip:"<<hitvec2.at(m).GetHitNumPad()<<" "<<hitvec2.at(m).GetHitNumStrip()<<std::endl;


	numHit_t2++;
	if((hitvec2.at(m).GetHitNumPad()>0)&&(hitvec2.at(m).GetHitNumStrip()==0))
	  numPadOnlyHit_t2++;

	if((hitvec2.at(m).GetHitNumPad()==0)&&(hitvec2.at(m).GetHitNumStrip()>0))
	  numStripOnlyHit_t2++;

	if((hitvec2.at(m).GetHitNumPad()>0)&&(hitvec2.at(m).GetHitNumStrip()>0))
	  {
	    numCl1HitHit_t2++;
	    hitvec2.at(m).GetHitPhi();
	    hitvec2.at(m).GetHitDPhi();
	  }

      }
      
    }

 
  /////////////////////////////////////////////////////////////////////////////////
  ///           Variable INITIALIZATION
  ///////////////////////////////////////////////////////////////////////////////


 int hemisphere = 0;

  unsigned int  sizeHitv=hitvec.size();
  if(sizeHitv<2)
    {
      std::cout<<" T2SelectionCutUtils::MyLinearfitCorr problem: Track with less than 2 Cl1 hits!! Continue the fitting.."<<std::endl;        }
  else
    hemisphere = (int)(hitvec[0].GetHitZ()/fabs(hitvec[0].GetHitZ()));

  
  TMatrixD ParCov_matr(4,4); //matrice di covarianza dei parametri trovati;
  unsigned int sizeArighe=sizeHitv*2;
  TMatrixD A(sizeArighe,4);
  TMatrixD At(4,sizeArighe);
  //A ?? la matrice per cui Mis=A(Param)
  TMatrixD Vy(sizeArighe,sizeArighe); //matrice di covarianza delle misure (una per ogni xy, quindi ?? diag a blocchi);
  TMatrixD Vym1(sizeArighe,sizeArighe); 

  TMatrixD Ls(4,4);//((A^T)(Vy^-1)A)
  // TMatrixD Cs(4,4);//(A^T)(Vy^-1)
  TMatrixD Cs(4,sizeArighe);//(A^T)(Vy^-1)

  TVectorD FittedParam(4);
 
  TVectorD Yvect(sizeArighe);


  Ls.Zero();
  TMatrixD Ls00(2,2);//un quarto della matrice LS.
  TMatrixD Ls01(2,2);
  TMatrixD Ls10(2,2);
  TMatrixD Ls11(2,2);

  TMatrixD Mia(2,2);
  TMatrixD Mib(2,2);

  TMatrixD MiaT(2,2);
  TMatrixD MibT(2,2);

  Ls00.Zero();
  Ls01.Zero();
  Ls10.Zero();
  Ls11.Zero();

  Vy.Zero();
  Vym1.Zero();

  if(verbosity)
    std::cout<<"MyLinearfitCorr Start computation  .. "<<std::endl;
  
  std::vector<std::vector<double> > All_YMeasures;
  std::vector<std::vector<double> > All_Vym1Measures;//11 12 21 22
  std::vector<std::vector<double> > All_VyMeasures;
  std::vector<double> All_ZMeasures; 





  for(unsigned int k =0; k<hitvec.size(); k++)
    {
      TMatrixD OneVy(2,2); 
      OneVy.Zero();
      
      phirad=hitvec[k].GetHitPhi()*3.14159/180.0;  
      sigmaPhi=hitvec[k].GetHitDPhi()*3.14159/180.0;
   
      r=hitvec[k].GetHitR();
      if(r<1)
	std::cout<<"Error: R hit too small:"<<r<<std::endl;

      sigmaR=hitvec[k].GetHitDR();

    //   
       //       std::cout<<"-DEBTRK sigmaR-R-phi-SigmaPhi:"<<sigmaR<<" "<<r<<" "<<phirad*180.0/3.14159<<" "<<sigmaPhi*180.0/3.14159<<" NumStrip: "<<hitvec[k].GetHitNumStrip()<<" NumPad:"<<hitvec[k].GetHitNumPad()<<std::endl;
 
   /////////////////////////////////////////////////////////////////////////////////
   ///         FITTING PART
   ///////////////////////////////////////////////////////////////////////////////



      OneVy(0,0)=sigmaR*sigmaR*cos(phirad)*cos(phirad)+r*r*sin(phirad)*sin(phirad)*sigmaPhi*sigmaPhi;//ex
      OneVy(0,1)=cos(phirad)*sin(phirad)*(sigmaR*sigmaR-r*r*sigmaPhi*sigmaPhi);  
      OneVy(1,0)=OneVy(0,1);
      OneVy(1,1)=sigmaR*sigmaR*sin(phirad)*sin(phirad)+ r*r*cos(phirad)*cos(phirad)*sigmaPhi*sigmaPhi;//ey
   
   //Invert Vy matrix
      TMatrixD OneVym1(2,2); 
      OneVym1=OneVy;
      Double_t deti;	
      OneVym1.Invert(&deti);
   
      //OneVy.Print();
      if(fabs(deti)<0.001)
	{
	  std::cout<<"Possible Vy Zero Determinant in one point error matrix!!!"<<std::endl;
	  std::cout<<"sigmaR-R-phi-SigmaPhi:"<<sigmaR<<" "<<r<<" "<<hitvec[k].GetHitR()<<" "<<phirad*180.0/3.14159<<" "<<sigmaPhi*180.0/3.14159<<" NumStrip: "<<hitvec[k].GetHitNumStrip()<<" NumPad:"<<hitvec[k].GetHitNumPad()<<std::endl;
	}

   Vym1(2*k,2*k)=OneVym1(0,0);
   Vym1(2*k,2*k+1)= OneVym1(0,1); 
   Vym1(2*k+1,2*k)=OneVym1(1,0);
   Vym1(2*k+1,2*k+1)=OneVym1(1,1);
   

   Vy(2*k,2*k)=OneVy(0,0);
   Vy(2*k,2*k+1)= OneVy(0,1); 
   Vy(2*k+1,2*k)=OneVy(1,0);
   Vy(2*k+1,2*k+1)=OneVy(1,1);

   
   std::vector<double>recYvect; 
   recYvect.push_back(hitvec[k].GetHitX());
   recYvect.push_back(hitvec[k].GetHitY());

   std::vector<double>  recVy;
   recVy.push_back(OneVy(0,0)); recVy.push_back(OneVy(0,1)); 
   recVy.push_back(OneVy(1,0)); recVy.push_back(OneVy(1,1));

   std::vector<double>  recVym1; 
   recVym1.push_back(OneVym1(0,0)); recVym1.push_back(OneVym1(0,1)); 
   recVym1.push_back(OneVym1(1,0)); recVym1.push_back(OneVym1(1,1));
   // std::vector<double> recZvect; recZvect.push_back(hitvec[k].GetHitZ());

   All_YMeasures.push_back(recYvect);
   All_Vym1Measures.push_back(recVym1);
   All_VyMeasures.push_back(recVy);
   All_ZMeasures.push_back(hitvec[k].GetHitZ());

   
   



   Yvect(2*k)=hitvec[k].GetHitX();
   Yvect(2*k+1)=hitvec[k].GetHitY();  

   A(2*k,0)=hitvec[k].GetHitZ();
   A(2*k,1)=1.0;
   A(2*k,2)=0.0;
   A(2*k,3)=0.0;
   A(2*k+1,0)=0.0;
   A(2*k+1,1)=0.0;
   A(2*k+1,2)=hitvec[k].GetHitZ();
   A(2*k+1,3)=1.0;

     

   Mia.Zero();
   Mib.Zero();
   
   
   Mia(0,0)=hitvec[k].GetHitZ();
   Mia(0,1)=0.;
   Mia(1,0)=1.;
   Mia(1,1)=0.;
   Mib(0,0)=0.;
   Mib(0,1)=hitvec[k].GetHitZ();
   Mib(1,0)=0.;
   Mib(1,1)=1.;
   
   MiaT=Mia;
   MibT=Mib;
   MiaT.Transpose(MiaT);
   MibT.Transpose(MibT);  
      

   Ls00=Ls00+Mia*OneVym1*MiaT;
   
   Ls01=Ls01+Mia*OneVym1*MibT;
   
   Ls10=Ls10+Mib*OneVym1*MiaT;
   
   Ls11=Ls11+Mib*OneVym1*MibT;  


 } //end loop on hits.

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

 if(verbosity)
   std::cout<<"MyLinearfitCorr Parameter calculation  .. "<<std::endl;

//Nota che qui Vy ?? in realt?? Vy^-1 

 Ls(0,0)=Ls00(0,0);  Ls(0,1)=Ls00(0,1);   Ls(0,2)=Ls01(0,0);   Ls(0,3)=Ls01(0,1);
 Ls(1,0)=Ls00(1,0);  Ls(1,1)=Ls00(1,1);   Ls(1,2)=Ls01(1,0);   Ls(1,3)=Ls01(1,1); 

 Ls(2,0)=Ls10(0,0);  Ls(2,1)=Ls10(0,1);   Ls(2,2)=Ls11(0,0);   Ls(2,3)=Ls11(0,1);
 Ls(3,0)=Ls10(1,0);  Ls(3,1)=Ls10(1,1);   Ls(3,2)=Ls11(1,0);   Ls(3,3)=Ls11(1,1);

 
 Double_t determ;	
 Ls.Invert(&determ);
 if(fabs(determ)<0.001)
     std::cout<<"WARNING: Possible Zero LS  Determinant!!!"<<std::endl;

 At.Zero();
 At.Transpose(A);
 Cs.Zero();
 Cs=At*Vym1;


 
 FittedParam= Ls*Cs*Yvect;//ax,bx.ay.by..
 
 //std::cout<<"MyLinearfitCorr FittedParam: "<<std::endl;
 //FittedParam.Print();

 ParCov_matr=At*Vym1*A;


 ParCov_matr.Invert(&determ);
 //std::cout<<"MyLinearfitCorr Error Matrix: "<<std::endl;
 //ParCov_matr.Print();

 if(fabs(determ)<0.001)
   std::cout<<"WARNING: Possible Zero ParCov_matr  Determinant!!!"<<std::endl;  

 //Calcolo chi2
 double chi2=0.;
 TVectorD Residui(sizeArighe);
 Residui=(Yvect-A*FittedParam);
 
 //TVectorD ResiduiT;
 //ResiduiT.Transpose(Residui);
 TVectorD VdotRes(sizeArighe);
 VdotRes=(Vym1*Residui);
 chi2=Residui*VdotRes;

 
 

 double chi2X=0.; double chi2Y=0.; 

 

   


 //////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////// SAVE THE OUTPUT
 /////////////////////////////////////////////////////////////////////////////////////

 //std::cout<<" DEBTRK SAVE THE Trk-1"<<std::endl;
 
 if(chi2<0)
   {
     std::cout<<"WARNING: MyLinearfitCorr Chi2: "<<chi2<<std::endl;
     std::cout<<"Residui"<<std::endl;
     Residui.Print();
     std::cout<<"Cov Matrix Vy"<<std::endl;
     Vy.Print();
     std::cout<<"Vym1*Residui"<<std::endl;
     VdotRes.Print();

   }
 par_covariance_matrix=ParCov_matr;
 
 //std::cout<<"MyLinearfitCorr Error Matrix: "<<std::endl;
 //par_covariance_matrix.Print();


 std::vector<double> vect;

 for (unsigned int i=0;i<10;i++)
   vect.push_back(0.);

 vect[0]=FittedParam(0);
 vect[1]=FittedParam(1);
 vect[2]=FittedParam(2);
 vect[3]=FittedParam(3);
 vect[4]=ParCov_matr(0,0);
 vect[5]=ParCov_matr(1,1);
 vect[6]=ParCov_matr(2,2);
 vect[7]=ParCov_matr(3,3);

 double correlabX=ParCov_matr(0,1);
 vect[8]=correlabX;
 double correlabY=ParCov_matr(2,3);
 vect[9]=correlabY;


 //Swap of positions in order to be compatible with the track object.
 TMatrixD ParCov_matrConv(4,4); 
 TVectorD FittedParamConv(4);
 

 for(unsigned int yy=0;yy<4;yy++)
   for(unsigned int ll=0;ll<4;ll++)
     ParCov_matrConv(yy,ll)=0.;


 ParCov_matrConv(0,0)=ParCov_matr(1,1);
 ParCov_matrConv(1,1)=ParCov_matr(3,3);
 ParCov_matrConv(2,2)=ParCov_matr(0,0);
 ParCov_matrConv(3,3)=ParCov_matr(2,2);

 FittedParamConv(0)=FittedParam(1);
 FittedParamConv(1)=FittedParam(3);
 FittedParamConv(2)=FittedParam(0);
 FittedParamConv(3)=FittedParam(2);

 //std::cout<<" DEBTRK SAVE THE Trk-2"<<std::endl;
  
  T1T2Track trk2rz=RPhiFit(hitvec);//RZFIT  
  
  double chi2Phi=trk2rz.ChiSquaredPhi();
  double chi2R=trk2rz.ChiSquaredR();
  
  double phiRZ=trk2rz._phiRZ;               //trk2rz.Phi(); 
  double e_phiRZ=trk2rz._e_phiRZ;
  
  double TanthetaRZ=trk2rz._TanthetaRZ;     //trk2rz.GetTy(); 
  double e_TanthetaRZ=trk2rz._e_TanthetaRZ; //trk2rz.GetTySigma();
  
  double bRZ=trk2rz._bRZ;          //trk2rz.GetTx(); 
  double e_bRZ=trk2rz._e_bRZ;      //trk2rz.GetTxSigma();
  
  if((fabs(vect[0])>1.8)||(fabs(vect[2])>1.8))
    {
      std::cout<<"WARNING in T2TrackProducer3 Ax:"<<(vect[0])<<"  Ay:"<<(vect[2])<<" too big!! Trk Hit (all) print :"<<std::endl;
      for(unsigned int jj =0; jj<hitvec.size(); jj++)
	{
	  std::cout<<"x-y-z-phi:  "<<hitvec[jj].GetHitX()<<"  "<<hitvec[jj].GetHitY()<<"  "<<hitvec[jj].GetHitZ()<<"  "<<hitvec[jj].GetHitPhi()<<"     DX:"<<hitvec[jj].GetHitDX()<<"     DY:"<<hitvec[jj].GetHitDY()<<" NumPad:"<<hitvec[jj].GetHitNumPad()<<std::endl;
	}
    }

  if(verbosity){
    std::cout<<"Before save: XY:"<<std::endl;
    // std::cout<<"ax.bx.ay.by:"<<vect[0]<<" "<<vect[1]<<" "<<vect[2]<<" "<<vect[3]<<std::endl;
  
    std::cout<<"Before save: RZ:"<<std::endl;
    std::cout<<"Phi "<<phiRZ*180.0/3.14159<<" TanTheta"<<TanthetaRZ<<" Brz:"<<bRZ<<" Eta:"<<trk2rz.Eta()<<std::endl;

    std::cout<<" TotHit:"<<numHit_t2<<" cl1Hit:"<<numCl1HitHit_t2<<" StripOnly:"<<numStripOnlyHit_t2<<" PadOnly:"<<numPadOnlyHit_t2<<std::endl;
  }


  T1T2Track fittedtrack(FittedParamConv,ParCov_matrConv,chi2,chi2X,chi2Y,hemisphere,2 , TanthetaRZ,  bRZ,   phiRZ,   e_TanthetaRZ,  e_bRZ,  e_phiRZ, chi2R, chi2Phi, RoadID, numHit_t2, numCl1HitHit_t2, numStripOnlyHit_t2, numPadOnlyHit_t2);
  

  if(verbosity){      
    std::cout<<"After save: :"<<std::endl;
    std::cout<<"PhiRZ:"<<fittedtrack._phiRZ*180.0/3.14159<<"  EtaRZ:"<< fittedtrack._etaRZ<<"  EtaXY:"<<fittedtrack.Eta()<<" ChirR/N:"<<(chi2R/(numCl1HitHit_t2+numPadOnlyHit_t2))<<" Chi2/N; "<<chi2/(numCl1HitHit_t2+numPadOnlyHit_t2)<<" NumHit:"<<fittedtrack._numHit_t2<<std::endl;
    
  }


  sizeHitv=hitvec.size();
  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      fittedtrack.AddHit(hitvec[jj]);
    }
  //std::cout<<" DEBTRK SAVE THE Trk size:"<<sizeHitv<<" hit saved"<<std::endl;
 // std::cout<<"Trk with "<<sizeHitv<<" hit saved"<<std::endl;
  return fittedtrack;
  
}







T1T2Track T2BackgroundAn::CalculateRecoTrkFromGeant(std::vector<Local3DPoint> PrimaryGeantHitPos, std::vector<double> &theZ,  std::vector<uint32_t> &thecmsswId){
  
  // std::cout<<"Inside CalculateRecoTrkFromGeantA"<<PrimaryGeantHitPos.size()<<std::endl;
  
    TMatrixD covmat(4,4);
    covmat.Zero();
    double chi2corr=0.;		
    
    T2Hit worsthit;
    unsigned int roadId=0;//theroad.RoadID;
    std::vector<T2Hit> hitvector;
    for(unsigned int i=0;i<PrimaryGeantHitPos.size();i++)
      {
	double X=PrimaryGeantHitPos.at(i).x();
	double Y=PrimaryGeantHitPos.at(i).y();
	double Z=theZ.at(i);
	double EX=2.;
	double EY=2.;
	double EZ=1.;
	uint32_t mycmsswid=thecmsswId.at(i);
	T2Hit arechit(X,  Y,  Z,  EX,  EY, EZ, mycmsswid);
	arechit.SetHitClass(1);
	arechit.SetHitNumStrip(1);
	arechit.SetHitNumPad(1);
	//set other properties of pad-strip composition
	hitvector.push_back(arechit);
      }
    // std::cout<<"Inside CalculateRecoTrkFromGeantB"<<hitvector.size()<<std::endl;
   
    T1T2Track firstTrkRaw=MyLinearfitCorrDEV(hitvector,covmat,chi2corr,false,roadId);  
    return firstTrkRaw;
    //std::cout<<"CDEBUG"<<std::endl;
}



//Take a Geant track and return the info about if it is not Reco (0) or is Reco from 1 Trk (1) or is Reco with more than one Trk (2)
//-1 means not enough geant hit in the geant trk or sw bug 
int  T2BackgroundAn::GeantTrkReconstructed(const edm::Event& event,
					   const std::auto_ptr<edm::PSimHitContainer>& theSimHits,
					   const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
					   const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
					   const std::map<unsigned int, std::list<Local3DPoint> >& trackHitList,
					   const std::auto_ptr<T1T2TrackCollection> &theT2Tracks2,
					   int GeantTrkIndex,
					   bool printActive,
					   unsigned int &NumGeantTrk_PassingCuts,
					   unsigned int &NumGeantTrk_PassingCuts_56Division,
					   unsigned int &NumGeantTrk_PassingCutsDiscrepancyZCutCond,
					   unsigned int &NumGeantTrk_PassingCutsDiscrepancy,
					   unsigned int &NumGeantTrk_PassingCuts_2,
					   /* std::vector<Local3DPoint> &PrimaryGeantHitPos,*/
					   T1T2Track &TheGeantTrkReconstr,bool &GeantTrkFound, double &effiasRecRec,int &BestTrkIdMatchedPosition, double  &lastPassingZmin,double  &lastPassingZImpact,
					   bool discrepancyStudies				   
					   )//discrepancyStudies is false as default ,int binmult	
{

  // std::cout <<"Geant Trk index "<< GeantTrkIndex <<" analysis"<< std::endl;
  if (verbosity > 1){ 
    std::cout << "Matching simhits to recotracks.." << std::endl;
  }


  bool bugflag=false;
  std::vector<Local3DPoint> PrimaryGeantHitPos;
  PrimaryGeantHitPos.clear();
  NumGeantTrk_PassingCutsDiscrepancyZCutCond=0;
  NumGeantTrk_PassingCutsDiscrepancy=0;
  NumGeantTrk_PassingCuts=0;
  NumGeantTrk_PassingCuts_2=0;
  //Container for primary and secondary tracks
  auto_ptr<T1T2TrackCollection> theT2PrimaryTracks (new T1T2TrackCollection());
  auto_ptr<T1T2TrackCollection> theT2SecondaryTracks (new T1T2TrackCollection());

  //Container for simulated primary tracks found 
  std::map<unsigned int, bool> reconstructedTrackForSimulatedPrimaryTrackFound;
  std::map<unsigned int, bool> reconstructedTrackForSimulatedSecondaryTrackFound;

  //Container for all simulated positions
  //pair<row, column>, pair<trackId,isPrimary>
  std::vector<std::map<std::pair<int, int>, std::pair<unsigned int,bool> > > simPadHitPositionsInPlanes(40);

  std::vector<int> PadRowHitsGeantTrk; 

  std::vector<int> PadColHitsGeantTrk;

  //std::vector<double> PadYHitsGeantTrk; std::vector<double> PadXHitsGeantTrk; std::vector<double> PadZHitsGeantTrk;
  std::vector<bool> IsPrimInfoHitsGeantTrk;
  std::vector<int> PadSymbdet_vect;
  //Convert simhits into information about plane, strip and pad.
  T2Geometry t2geom;
  T2GeometryUtil conv;     
  unsigned int geantHitCount=0;
  
  if (verbosity > 1){
    std::cout << "Analysing simhits.." << std::endl;
  }

 
  std::vector<unsigned int> GeantHiCmsswId;//Utilized for discrepancy studies.
  std::vector<double> GeantHitRefZ;//Utilized for discrepancy studies.
  
 std::vector<int> SymbplaneInGTrk;
  //Load simulated hits and associate for each corrisponding
  // T2-pad the  information about primary/secondary which generates it
  if(printActive)
    if (verbosity >= 1)
      std::cout<<"Info from SIMHIT in Wanted GTrkRecon:"<<std::endl;
  int padRowGeant=-111;

int padRow=-111;
 int padCol=-111; int trackId=0;
  for(edm::PSimHitContainer::const_iterator simHitIt = theSimHits->begin(); simHitIt != theSimHits->end(); ++simHitIt){

    trackId = simHitIt->trackId();

    if(GeantTrkIndex==trackId)
      {
	geantHitCount++;	
	int detUnitId = simHitIt->detUnitId();    
	if(detUnitId<0){
	  std::cout<<"ERROR detUnitId<0 !!! Forcing a crash.."<<std::endl;
	  std::vector<int> crashh;
	  crashh.at(100)=1;
	  std::cout<<crashh.at(1)<<std::endl;
	}
	//uint32_t
	simhitUnitID->Fill(detUnitId/10000);
	T2GeometryUtil::T2DetInfo detInfo = conv.GetT2Info(detUnitId);
	t2geom.setPlane(detUnitId); //getNearestPadRo and col does not use the det info. Are refered to one particular plane.
    
	Local3DPoint hitPos = simHitIt->localPosition();    
	Local3DPoint corPos = CalculateCorrectXYZPos(Local3DPoint(hitPos.x(),hitPos.y(),hitPos.z()),detUnitId);

	//padRow = t2geom.getNearestPadRo(&corPos);
	//padCol = t2geom.getNearestPadCo(&corPos); 
	double radiusGhit=sqrt(hitPos.x()*hitPos.x()+hitPos.y()*hitPos.y()); 
	double phiGHit=atan(fabs(hitPos.y())/fabs(hitPos.x()));//atan2(hitPos.y(),hitPos.x());

	
	phiGHit=phiGHit*180./3.14159;


	if((corPos.y()<0)&&(corPos.x()>0))
	  phiGHit=360.0-phiGHit;
	
	if((corPos.y()>0)&&(corPos.x()<0))
	  phiGHit=180.0-phiGHit;

	if((corPos.y()<0)&&(corPos.x()<0))
	  phiGHit=phiGHit+180.;

	padRow = t2geom.getNearestPadRow2_(&corPos,detUnitId);
	padCol =t2geom.getNearestPadCol_(radiusGhit,phiGHit,detUnitId);
	
	padRowGeant=padRow;


	//padRowReal = t2geom.getNearestPadRow(&hitPos);
	//padColReal = t2geom.getNearestPadCol(&hitPos);

	int padSymbDET =detInfo.symb;
	PadRowHitsGeantTrk.push_back(padRow); 
	PadColHitsGeantTrk.push_back(padCol); 
	PadSymbdet_vect.push_back(padSymbDET);
	
	if((padSymbDET/10)==0){
	  Q0PadRCPrimary->Fill(padRow,padCol);

	  if(padCol==32)
	    {
	      // RPhiRecHitinH0Col32->Fill();
	      RPhiGeantHitinH0Col32->Fill(radiusGhit,phiGHit);
	      XYGeantHitinH0Col32->Fill(corPos.x(),corPos.y());
	    }
	}
	if((padSymbDET/10)==1)
	  Q1PadRCPrimary->Fill(padRow,padCol);  
	

	if((padRow<0)||(padRow>=25))
	  {
	    std::cout<<"padRow out of range"<< padRow<<std::endl;
	    std::cout<<"Plane:"<<padSymbDET<<" Phi:"<<phiGHit<<" R:"<<radiusGhit<<std::endl;
	    //return -1;
	    bugflag=true;
	  }
	
	if((padCol<0)||(padCol>=65)){
	  std::cout<<"padCol out of range:"<< padCol<<"Plane:"<<padSymbDET<<" Phi:"<<phiGHit<<" R:"<<radiusGhit<<std::endl;
	  //return -1;
	  bugflag=true;
	}
	
	//return 
	
	
	  //It is assumed that the default track collection is the one you do the test
	  //The program only look at the first two and last two planes, this is assumed to be the reference.
	   int TherefSymbplane=RawtoSymb(detUnitId);
	   SymbplaneInGTrk.push_back(TherefSymbplane);
	   std::vector <int>::iterator iter;
	   //iter = std::find( refSymbplane.begin(), refSymbplane.end(), TherefSymbplane);
	   //Store the hitPosIn a vector.
	   //GeantHitRefX;
	   //GeantHitRefY;
	   //GeantHitRefZ;
	   double theZ=detInfo.Zdet;
	   GeantHitRefZ.push_back(theZ);
	   GeantHiCmsswId.push_back(detUnitId);
	   PrimaryGeantHitPos.push_back(corPos);
	 
	

      }// If GeantTrk==id
    
  }//End PSIMIHIT LOOP

  //std::cout<<"GeantTrkEffi Calculator: Simhit associated to original trk: "<<geantHitCount<<std::endl;
  //Iterating all reconstructed tracks to find the reco track matching with the input SimuTrk
  
  int numTrkRecontructingGeant=0; 
  if(verbosity>=1){
    int trksize=theT2Tracks2->begin()-theT2Tracks2->end();
    std::cout<<"Track size: "<<trksize<<std::endl;
  }


  double myeffiproj=0.;
  double myeffiprojRecCut=0.;

  if(PrimaryGeantHitPos.size()>=4){
     if(verbosity>1)
       std::cout<<"Calling From GeantTrkReconstructed the recon:"<<PrimaryGeantHitPos.size()<<std::endl;
    TheGeantTrkReconstr =CalculateRecoTrkFromGeant(PrimaryGeantHitPos,GeantHitRefZ,GeantHiCmsswId);
     if(verbosity>1)
       std::cout<<"Geant Trk recon with "<<TheGeantTrkReconstr.GetHitEntries()<<std::endl;
    GeantTrkFound=true;
  }
  else{
    GeantTrkFound=false;    
  }


  int myrefq=-1;//
  int dquarter=1;
  for(unsigned int i=0;i<SymbplaneInGTrk.size();i++)
    {
      if(i==0)
	myrefq=SymbplaneInGTrk.at(i)/10;
      else{
	//int quarter=SymbplaneInGTrk.at(i)/10;
	dquarter=(SymbplaneInGTrk.at(i)/10)-(SymbplaneInGTrk.at(i-1)/10);
      }
    }

 if(dquarter!=0)
   GeantTrkFound=false;//Means that are in 2 quarters

  std::vector<T1T2Track> allmatchTrk;
  int padrow2  =1;


  double lastbestchi2=10000.;


  for(T1T2TrackCollection::const_iterator recTrack = theT2Tracks2->begin(); recTrack!=theT2Tracks2->end();++recTrack){

    if(printActive)
      if(verbosity>=1)
	std::cout<<"Info from RECOTRK in GTrkRecon:"<<std::endl;
 
    unsigned int trkcounterPrim=0;

    int hitCount = recTrack->GetHitEntries();
    
    /*if(T2CutsUtil.TrkAlsoInQuarter((*recTrack),T2_QuarterUsed))*/
    
   
      for (int h = 0;h< hitCount;++h){
    
	T2Hit hit = recTrack->GetHitT2(h);
	int detUnitId = hit.GetHitDetRawId();

	rechitUnitID->Fill(hit.GetHitDetRawId()/10000);
	T2GeometryUtil::T2DetInfo detInfo2 = conv.GetT2Info(detUnitId);
      
	t2geom.setPlane(detUnitId);
	Local3DPoint corPos(hit.GetHitX(),hit.GetHitY(),hit.GetHitZ());
	//int padRow = t2geom.getNearestPadRo(&corPos);
	//int padCol = t2geom.getNearestPadCo(&corPos);

	
	int padRow = t2geom.getNearestPadRow2_(&corPos,detUnitId);
	int padCol = t2geom.getNearestPadCol_(/*&corPos,detInfo.cmsswid*/hit.GetHitR(),hit.GetHitPhi(),detUnitId);

  
	if((padRow<0)||(padRow>=25))
	  {
	    std::cout<<"padRow out of range2:"<< padRow<<std::endl;
	    bugflag=true;
	  }
	
	if((padCol<0)||(padCol>=65)){
	  std::cout<<"padCol out of range2:"<< padCol<<" Plane:"<<detInfo2.symb<<" Hit-Phi:"<<hit.GetHitPhi()<<" Hit-R:"<<hit.GetHitR()<<" Num Pad-Srip: "<<hit.GetHitNumPad()<<"-"<<hit.GetHitNumStrip()<<std::endl;
	  bugflag=true;
	}

	
	padrow2=padRow;

	//To be removed

	//if(printActive)
	//if(verbosity>=1)
	
	//	std::cout<<"PadR-PadCol-Plane in RECO Trk:"<<padRow<<" "<<padCol<<" "<<detInfo2.symb<<" || "<<corPos.x()<<" "<<corPos.y()<<" "<<corPos.z()<<std::endl;
	for(unsigned int u=0;u<PadRowHitsGeantTrk.size();u++)
	  {
	    //Very-Important-Parameter
	    if(PadSymbdet_vect.at(u) == (int) detInfo2.symb ){
	      
	      

	      if((padRow>-1)&&(padCol>-1))
	      if((abs(padRow-PadRowHitsGeantTrk.at(u))<=2)&&(abs(padCol-PadColHitsGeantTrk.at(u))<=2))
		{
		  //if(IsPrimInfoHitsGeantTrk.at(u)==true) 
		  //This above is not needed since the control about primary/sec shold be done at calling level. 
		  trkcounterPrim++;
		  //if((padRow==0)||(padRow==1))
		  //GeantTrkFound=false;

		  //	else
		  //trkcounterSec++;
		}
	    }
	  }
      }
    
      

      //This is very important. You require that the RecoTrk has at least 3 hit matching GEANT
      //Warning: the Normal RecoTrk are supposed to be the one of the testing plane.
      //Very-Important-Parameter
      if(((int)trkcounterPrim)>=numhitRequiredFormMatching)//We normally use 3. 1 is to be compatible with the discrepancyStudies  	
      {
	numTrkRecontructingGeant++;
	RecoTrkEta_MatchingPrimary->Fill(recTrack->Eta());
	T1T2Track atrk=(*recTrack);

	double thischi2=TMath::Prob((*recTrack).ChiSquared(),(*recTrack)._numCl1HitHit_t2*2-4);

	if(thischi2<lastbestchi2){
	  BestTrkIdMatchedPosition=recTrack-theT2Tracks2->begin();
	  lastbestchi2=thischi2;
	}

	if(RecoTrk_UnfoldingCut(atrk)){
	  NumGeantTrk_PassingCuts++;	
	  lastPassingZmin=atrk.Z_at_Rmin();
	  lastPassingZImpact=ZimpactFromRecTrack(atrk);
	}

	
	if(RecoTrk_UnfoldingCut_56Division(atrk))
	  NumGeantTrk_PassingCuts_56Division++;

	if(RecoTrk_UnfoldingCut_2(atrk))
	  NumGeantTrk_PassingCuts_2++;
      }    

   

  }//Reco Trk loop

   int trksize2=theT2Tracks2->begin()-theT2Tracks2->end();

   if(trksize2==1){
     if(myrefq==0){
       PadRowGeantRecoH0->Fill(padrow2-padRowGeant,padrow2);
       PadRowGeantRecoH0_bis->Fill(padrow2-padRowGeant,padRowGeant);	    
     }
     
     if(myrefq==1){
       PadRowGeantRecoH1->Fill(padrow2-padRowGeant,padrow2);
       PadRowGeantRecoH1_bis->Fill(padrow2-padRowGeant,padRowGeant);
     }
   }
  


  if (verbosity > 1){
    std::cout << "Done mathcing simhits to recotracks.Line3147" << std::endl << std::endl;
  }

  effiasRecRec=0.;
  if(myeffiproj>=1){
    effiasRecRec=1.;
  }else{
    effiasRecRec=0.;
  }

  if(myeffiprojRecCut>=1){
    effiasRecRec=3.;
  }
    
  //  std::cout<<"..Q"<<myrefq<<": T2-G4Hit associated to Beginning-TrkID: "<<GeantTrkIndex<<" #GeantHitInT2"<<geantHitCount<<" #RawTrk Associated:"<<numTrkRecontructingGeant<<" #CuttedRecoTrk:"<<NumGeantTrk_PassingCuts<<std::endl;
  
  if(geantHitCount<4)
    numTrkRecontructingGeant=-1;

  if(GeantTrkFound==false){
    std::cout<<"WARN Geant Trk not found in the efficalculator (maybe is an overlap) recGeantHitSize in Trk: "<<std::endl;
    numTrkRecontructingGeant=-1;
  }


if(PadSymbdet_vect.size()>3)
  {
    PadSymbdet_vect.at(0);
    /*
    for(unsigned int i=0;i<PadSymbdet_vect.size();i++){
      int actref=PadSymbdet_vect.at(i)/10;
      if(actref!=refquarterfirthistrk){
	//std::cout<<"Mixed quarter in Trk-GHits!"<<std::endl;
	numTrkRecontructingGeant=-1;
      }
    }
    */
 
  } 
 else{
   numTrkRecontructingGeant=-1;
   std::cout<<">>PadSymbdet_vect.size()<3"<<std::endl;
 }

 if(bugflag)
   numTrkRecontructingGeant=-1;

return numTrkRecontructingGeant;
}











double T2BackgroundAn::ZimpactFromRecTrack(T1T2Track &mytrk){

  double C0=(mytrk.GetHitT2(0).GetHitX()*mytrk.GetHitT2(0).GetHitX()+mytrk.GetHitT2(0).GetHitY()*mytrk.GetHitT2(0).GetHitY())/(mytrk.GetTx()*mytrk.GetHitT2(0).GetHitX()+mytrk.GetTy()*mytrk.GetHitT2(0).GetHitY());
  double Z0impact=mytrk.GetHitT2(0).GetHitZ()-C0;

  return Z0impact;
}



int T2BackgroundAn::getPythiaPID(const HepMC::GenEvent* evt,int &PyPid,int OldestmotherGeantId)
{

  for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
    if((*p)->barcode()==OldestmotherGeantId)
      {
	PyPid=(*p)->pdg_id();
      }
  }
  
  return PyPid;
}











unsigned int
T2BackgroundAn::loadTrackproducerRaw(const edm::Event &iEvent,
				     std::auto_ptr<T1T2TrackCollection> &theT2TracksRaw){


  theT2TracksRaw.get()->clear();
  t2trackVectorBothSideRAW.clear();

  edm::Handle<T1T2TrackCollection> myTrackCollection;
  
  // iEvent.getByLabel(TrackLabel2,"T2TrackColl2",myTrackCollection);
  iEvent.getByLabel(TrackLabel,"T2TrackColl",myTrackCollection);//TrackLabel=T2TrackColl3

  // unsigned int rawcounterchi01Point=0;
  //std::vector<int> dummy;
  //dummy.push_back(0);dummy.push_back(1);

  for(T1T2TrackCollection::const_iterator ittrack = myTrackCollection->begin();ittrack!=myTrackCollection->end();++ittrack){ 
    if(ittrack->GetHitEntries()>=4){

      t2trackVectorBothSideRAW.push_back((*ittrack));

      theT2TracksRaw->push_back(*ittrack);
    }
  }

return   theT2TracksRaw->size();
}


/*
  Reads all trackproducer 2 tracks into theT2Tracks
  Returns the number of tracks read.
 */
unsigned int
T2BackgroundAn::loadTrackproducer2Tracks(const edm::Event &iEvent,
					       std::auto_ptr<T1T2TrackCollection> &theT2Tracks2)
{ 
  std::cout << "loadTrackproducer2TracksIn"<<std::endl;
  
  int numtrkInit=0;

  T2SelectionCutUtils T2CutsUtile;
  theT2Tracks2.get()->clear();
  t2trackVectorBothSide.clear();

  edm::Handle<T1T2TrackCollection> myTrackCollection;
  
  // iEvent.getByLabel(TrackLabel2,"T2TrackColl2",myTrackCollection);
 
  iEvent.getByLabel(TrackLabel,"T2TrackColl",myTrackCollection);//TrackLabel=T2TrackColl3
if (verbosity > 1){
    std::cout << "Loading all calculated TRACK PRODUCER 2 tracks from root file...Load size" <<myTrackCollection->size()<< std::endl;
  }

 
  std::vector<int> OneHalf0;OneHalf0.push_back(0); 
  std::vector<int> OneHalf1;OneHalf1.push_back(1);
  std::vector<int> OneHalf2;OneHalf2.push_back(2);
  std::vector<int> OneHalf3;OneHalf3.push_back(3);

  for(T1T2TrackCollection::const_iterator ittrack = myTrackCollection->begin();ittrack!=myTrackCollection->end();++ittrack){
    if (verbosity > 1)
      std::cout<<"Trk with #Hit:"<<(*ittrack).GetHitEntries()<<std::endl;
    if(ittrack->GetHitEntries()>=4){
      numtrkInit++;
      T1T2Track debtrk=(*ittrack);double Zimp=ZimpactFromRecTrack(debtrk);      
      //std::cout << "muertoIn"<<std::endl;
      if((fabs(debtrk.Z_at_Rmin())<13500)||((debtrk.Z_at_Rmin())*debtrk.GetHitT2(0).GetHitZ()<0)){
	//std::cout << "muertoIn Zimp1:"<<Zimp<<" "<<debtrk.GetHitT2(0).GetHitZ()<<std::endl;
	if(T2CutsUtil.TrkAlsoInQuarter(debtrk,OneHalf0)){	
	  // std::cout << "muertoIn Zimp2 -0:"<<Zimp<<std::endl;		
	  CumulativeZimpH0->Fill(Zimp);
	}
	if(T2CutsUtil.TrkAlsoInQuarter(debtrk,OneHalf1)){ //std::cout << "muertoIn Zimp2 -1:"<<Zimp<<std::endl;			
	  CumulativeZimpH1->Fill(Zimp);
	}
	if(T2CutsUtil.TrkAlsoInQuarter(debtrk,OneHalf3)){ //std::cout << "muertoIn Zimp2 -2:"<<Zimp<<std::endl;			
	  CumulativeZimpH2->Fill(Zimp);
	}
	if(T2CutsUtil.TrkAlsoInQuarter(debtrk,OneHalf2)){ //std::cout << "muertoIn Zimp2 -3:"<<Zimp<<std::endl;			
	  CumulativeZimpH3->Fill(Zimp);
	}
      }
      //      std::cout << "muetoIn 2!!"<<std::endl;
    if((*ittrack).GetHitEntries()>=4)    
      t2trackVectorBothSide.push_back((*ittrack));
    //if(T2CutsUtile.ChiCutCond((*ittrack), true, 0.01, 0.01)) Removed on 15Ago
    theT2Tracks2->push_back(*ittrack);
    }    
  }  
  std::cout << "loadTrackproducer2TracksOut"<<std::endl;

  if (verbosity > 1){
    std::cout << "Done reading TRACK PRODUCER 2 tracks. # Trks:" << theT2Tracks2->size() << " while at beginning:"<<numtrkInit<< std::endl;    
  }

  return theT2Tracks2->size();
}



/*
  Reads all simulated tracks into theSimTracks.
  Returns number of tracks read.
 */ 
unsigned int
T2BackgroundAn::loadSimulatedTracks(const edm::Event& iEvent,
					  std::auto_ptr<edm::SimTrackContainer>& theSimTracks)
{
  if (verbosity > 1){
    std::cout << "Loading all simtracks from root file..." << std::endl;
  }

  edm::Handle<edm::SimTrackContainer> simTrackCollection;
  iEvent.getByLabel(SimTrackContainerLabel, simTrackCollection);
  theSimTracks.get()->clear();

  for(edm::SimTrackContainer::const_iterator ittrack = simTrackCollection->begin();ittrack != simTrackCollection->end(); ++ittrack){ 
    SimTrack thetrack=(*ittrack);
    theSimTracks->push_back(thetrack);
  }

  if (verbosity > 1){
    std::cout << "Done reading simtracks." << std::endl << std::endl;
  }

  return theSimTracks->size();
}

/*
  Reads all simulated vertices into theSimVertices.
  Returns number of vertices read.
 */
unsigned int
T2BackgroundAn::loadSimulatedVertices(const edm::Event& iEvent,
					    std::auto_ptr<edm::SimVertexContainer>& theSimVertices)
{
  if (verbosity > 1){
    std::cout << "Loading all simvertices from root file..." << std::endl;
  }

  edm::Handle<edm::SimVertexContainer> simVertexCollection;
  iEvent.getByLabel(SimVertexContainerLabel, simVertexCollection);
  theSimVertices.get()->clear();

  for(edm::SimVertexContainer::const_iterator itvertex = simVertexCollection->begin();itvertex != simVertexCollection->end(); ++itvertex){ 
    SimVertex thevertex=(*itvertex);   
    theSimVertices->push_back(thevertex);
  }

  if (verbosity > 1){
    std::cout << "Done reading simvertices." << std::endl << std::endl;    
  }

  return theSimVertices->size(); //bbb

}

/*
  Reads all simulated hits into theSimHits.
  Returns number of simhits read.
 */
unsigned int 
T2BackgroundAn::loadSimulatedHits(const edm::Event& iEvent, 
					const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
					const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
					std::auto_ptr<edm::PSimHitContainer>& theSimHits,
					std::map<unsigned int, std::list<Local3DPoint> >& trackHitList,
					std::map<unsigned int, int >& unknownTrackList)
{
  if (verbosity > 1){
    std::cout << "Loading simulated hits..." << std::endl;    
  }

  edm::Handle<edm::PSimHitContainer> psimHitCollection;
  iEvent.getByLabel("g4SimHits","TotemHitsT2Gem",psimHitCollection);
  
  theSimHits.get()->clear();
  trackHitList.clear();
  
  for(edm::PSimHitContainer::const_iterator ithit = psimHitCollection->begin();ithit != psimHitCollection->end(); ++ithit){ 
    theSimHits->push_back(*ithit);  
    PSimHitT2_EnergyLoss->Fill(((*ithit).energyLoss())*1000.0*1000.0*1000.0); //EnergyLoss is in GeV, I want eV
    Local3DPoint localPosition = ithit->localPosition();

    unsigned int detUnitId = ithit->detUnitId();
    int particleType = ithit->particleType();
    unsigned int trackId = ithit->trackId();

    //If not found track with TrackId add it to unknown tracks with particle code
    bool trackIdFound = false;
    for(edm::SimTrackContainer::const_iterator ittrack = theSimTracks->begin();ittrack != theSimTracks->end(); ++ittrack){     
      if (ittrack->trackId() == trackId){
	trackIdFound = true;
	break;
      }
    }

    if (trackIdFound == false){
      //Check whether unknown track id has been inserted before
      std::map<unsigned int, int>::const_iterator ittrack = unknownTrackList.find(trackId);
      if (ittrack == unknownTrackList.end()){
	//Not inserted before
	unknownTrackList[trackId] = particleType;
	//Insert in PDG histograms
      } else {
	//Inserted before
	if (ittrack->second != particleType){
	  if (verbosity > 0){
	    std::cout << "!!! Known unknown track with different PDG !!!" << " Track ID: " << trackId << "  Old PDG: " << ittrack->second << " New PDG: " << particleType << std::endl
		      << "Skipped inserting unknown track." << std::endl;
	  }
	}
      }       
    }
    
    T2GeometryUtil conv;     
    T2GeometryUtil::T2DetInfo detInfo = conv.GetT2Info(detUnitId);    
    Local3DPoint corr_pos = CalculateCorrectXYZPos(Local3DPoint(localPosition.x(),localPosition.y(),detInfo.Zdet),detUnitId);
    
    trackHitList[trackId].push_back(corr_pos);    
  }

  if (verbosity > 1){
    std::cout << "Done loading simulated hits." << std::endl << std::endl;    
  }

  return theSimHits.get()->size();
 
}




double T2BackgroundAn::PartCharge(int pdgcode){
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




/*
T2BackgroundAn::GetPadColR_PHI_Z(double &R,double &Phi, double &Z, int col, int row, uint32 rawdet_id)
{
  T2ROGeometry t2rogeo(rawdet_id);
  double myphimin=t2rogeo.GetPadPhiMin(row,col);
  double myphimax=t2rogeo.GetPadPhiMin(row,col);
  double myphiAvg=(myphimin+myphimax)/2.0;
  int pCol=-1;
  int pColRiprova=col;
  if(mydet.planeSide()==1)//Onother way to refer col numbering respect to the plane closer to the ip.
    pColRiprova=(64-col);

  std::cout<<col<<" "<<pColRiprova<<" | "<<std::endl;
  
  unsigned int ht=0;
  if(detnumb<10)
    ht=0;
  if(detnumb>=10)
    ht=1;

  T2DetId DetInRefPlane(0,ht,0,0);    //I refer the numbering respect to the first plane.
  T2ROGeometry t2roRefPlNearHS0(DetInRefPlane.rawId());
				 
  bool foundcol=false;
  int j=0;
  while((foundcol==false)&&(j<65))
    {
      double refphiminj=t2roRefPlNearHS0.GetPadPhiMin(row,j);//riga-colonna
      double refphimaxj=t2roRefPlNearHS0.GetPadPhiMax(row,j);//riga-colonna
      double refphiavg=(refphiminj+refphimaxj)/2.0;  //riga-colonna
      //	std::cout<<"??? row: "<<row<<"  Phi: "<<phi<<" Phimin: "<<phiminj<<" Phimax: "<<phimaxj<<std::endl;
      //There cpuld be a gap of about 0.1 between pad. Try to look for this->enlarge pad azimuthal area.
      if((myphiAvg>=(refphiminj-0.1))&&(myphiAvg<=(refphimaxj+0.1)))
	{
	  pCol=j;
	  foundcol=true;
	}
      else
	{
	  if(fabs(refphiminj-refphimaxj)>350.) // This can happen: refphiminj: 358.534 refphimaxj: 1.41393
	    {
	      if((fabs(myphiAvg-refphimaxj)<2.0)||(fabs(myphiAvg-refphiminj)<2.0))
		{
		  pCol=j;
		  foundcol=true;
		}
	    }
	}
      j++;
    }
  

  if(pCol==-1)
    {
      std::cout<<"Error in getNearestPadCol_ : col not found !!!!. Report of the problem:";
      bool foundcol=false;
      unsigned int j=0;
      while((foundcol==false)&&(j<65))
	{
	  double phiminj=t2rogeo.GetPadPhiMin(row,j);//riga-colonna
	  double phimaxj=t2rogeo.GetPadPhiMax(row,j);//riga-colonna
	  
	  double refphiminj=t2roRefPlNearHS0.GetPadPhiMin(row,j);//riga-colonna
	  double refphimaxj=t2roRefPlNearHS0.GetPadPhiMax(row,j);//riga-colonna
	  double refphiavg=(refphiminj+refphimaxj)/2.0;  //riga-colonna
	  //	std::cout<<"??? row: "<<row<<"  Phi: "<<phi<<" Phimin: "<<phiminj<<" Phimax: "<<phimaxj<<std::endl;
	  //There cpuld be a gap of about 0.1 between pad. Try to look for this->enlarge pad azimuthal area.
	  if((myphiAvg>=(refphiminj-0.1))&&(myphiAvg<=(refphimaxj+0.1)))
	    {
	      pCol=j;
	      foundcol=true;
	    }
	  else
	    {
	      if(fabs(refphiminj-refphimaxj)>350.) // This can happen: refphiminj: 358.534 refphimaxj: 1.41393
		{
		  
		  if((fabs(myphiAvg-refphimaxj)<2.0)||(fabs(myphiAvg-refphiminj)<2.0))
		    {
		      pCol=j;
		      foundcol=true;
		    }
		}
	    }
	  
	  std::cout<<"DEB: row: "<<row<<"  MyPhiAvg: "<< myphiAvg<<" Phimin: "<<refphiminj-0.1<<" Phimax: "<<refphimaxj+0.1<<std::endl;
	  
	  j++;
	}
      
    }
  else //You can Fill your Histograms
    {
      std::cout<<detnumb<<" "<<pCol<<" "<<row<<std::endl;
      OneEventPlanePadRowColConv[detnumb]->Fill(pCol,row);
      OneEventPlanePadRowCol[detnumb]->Fill(pColRiprova,row);
      if(detnumb<10)
	OneEventPlanePadRowColGlob->Fill(pColRiprova,row);
    }
}
*/

unsigned int T2BackgroundAn::RawtoSymb(uint32_t thedet)
{
  T2DetId converter;
  unsigned int pl=converter.plane(thedet);
  unsigned int pls=converter.planeSide(thedet);
  unsigned int ht=converter.halfTelescope(thedet);
  unsigned int arm=converter.arm(thedet);
  unsigned int symbolic=pl*2+pls+ht*10+20*arm;	  
 
  return symbolic;
}


double T2BackgroundAn::GeneratorEtaFromTrkBarcode(int barcode,const HepMC::GenEvent* evt)
{

  double geneta=0.; 

  for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
    if((*p)->barcode()==barcode)
      {
	geneta=((*p)->momentum().eta());
      }
  }

  return geneta;
  
}


/*
  Calculates correct XYZ pos from simhit pos
 */
Local3DPoint
T2BackgroundAn::CalculateCorrectXYZPos(const Local3DPoint &pos, const unsigned int detUnitId)
{
  T2GeometryUtil conv;     
  T2GeometryUtil::T2DetInfo detInfo = conv.GetT2Info(detUnitId);  

  double corr_x_pos = pos.x();
  double corr_y_pos = pos.y();
  double corr_z_pos = pos.z();

  if (detInfo.ht == 1 && detInfo.plside == 1){
    corr_x_pos = -corr_x_pos;      
     if (detInfo.arm == 0){
      corr_y_pos = -corr_y_pos;       
       }
  } else if (detInfo.ht == 1 && detInfo.plside == 0){
    corr_x_pos = -corr_x_pos;
    if (detInfo.arm == 1){
      corr_y_pos = -corr_y_pos;       
        }
  } else if (detInfo.ht == 0 && detInfo.plside == 1){  //REFERENCE
     if (detInfo.arm == 1){
     corr_y_pos = -corr_y_pos;       
      }
  } else if (detInfo.ht == 0 && detInfo.plside == 0){
    if (detInfo.arm == 0){
      corr_y_pos = -corr_y_pos;       
     }
  }    
  
  return Local3DPoint(corr_x_pos,corr_y_pos,corr_z_pos);
}

/*
  Returns track corresponding to trackId 
  Returns true if found, false otherwise.
 */
bool
T2BackgroundAn::getTrack(const unsigned int trackId,  const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, SimTrack &track)
{
  if(trackId < 1){
    return false;
  }

  //Guesss
  if (trackId <= theSimTracks.get()->size()){
    track = (*theSimTracks.get())[trackId-1];
    if (track.trackId() == trackId){
      return true;
    }
  }

  //Iterate all
  for(edm::SimTrackContainer::const_iterator ittrack = theSimTracks->begin();ittrack != theSimTracks->end(); ++ittrack){     
    if (ittrack->trackId() == trackId){
      track  = (*ittrack);
      return true;
    }
  }
  
  return false;  
}

/*
  Returns vertex corresponding to vertexId 
  Returns true if found, false otherwise.
*/
bool
T2BackgroundAn::getVertex(const  int vertexId,  const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, SimVertex &vertex)
{
  if (vertexId < (int) theSimVertices.get()->size() && vertexId >=0){
    vertex = (*theSimVertices.get())[vertexId];
    return true;
  }
  return false;  
}
/*
bool
T2BackgroundAn::getVertex(const unsigned int vertexId,  const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, SimVertex &vertex)
{
  if (vertexId < theSimVertices.get()->size() && vertexId >=0){
    vertex = (*theSimVertices.get())[vertexId];
    return true;
  }
  return false;  
}
*/





bool 
T2BackgroundAn::isTrackPrimary(const int trackId,
					 const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
					 const std::auto_ptr<edm::SimVertexContainer>& theSimVertices)
{
  SimTrack track;
  if (getTrack(trackId,theSimTracks,track) == true){
    SimVertex vertex;    
    if (getVertex(track.vertIndex(), theSimVertices,vertex) == true){
      if (vertex.parentIndex() < 0){
	return true;
      }
    }
  }
  return false;
}





double T2BackgroundAn::CalculateEfficiencyAsRecRec(T1T2Track &TheGeantTrkReconstr, T1T2Track &atrk)
{
  
  T2SelectionCutUtils T2CutUtils;





  


  


  




  


  unsigned int sizeref=TheGeantTrkReconstr.GetT2TrackHits().size();
  unsigned int sizetest=atrk.GetT2TrackHits().size();
  bool TrkMatch=false;
  double trkX1=0.;
  double trkY1=0.;
  double trkR1=0.;
  for(unsigned int jj=0;jj<sizeref;jj++)
    {
      trkR1+=TheGeantTrkReconstr.GetHitT2(jj).GetHitR();
      trkX1+=TheGeantTrkReconstr.GetHitT2(jj).GetHitX();
      trkY1+=TheGeantTrkReconstr.GetHitT2(jj).GetHitY();  
    }
  
  trkR1=trkR1/((double)sizeref);  
  trkY1=trkY1/((double)sizeref);  
  trkX1=trkX1/((double)sizeref);  
  
  for(unsigned kk=0;kk<sizetest;kk++){
    double z=atrk.GetHitT2(kk).GetHitZ();	    	    
    double xp=z*TheGeantTrkReconstr.GetTx()+TheGeantTrkReconstr.X0();
    double yp=z*TheGeantTrkReconstr.GetTy()+TheGeantTrkReconstr.Y0();
    trkR1=sqrt(xp*xp+yp*yp);
    
    
   T2CutUtils.PhiFitAverage(atrk.GetT2TrackHits());
    double trkR2=atrk.GetHitT2(kk).GetHitR();//TestArmXTrk->at(j).GetHitT2(0).GetHitR();
    double TrkDY=yp-atrk.GetHitT2(kk).GetHitY();//TestArmXTrk->at(j).GetHitT2(0).GetHitY();//CleanRefHalfXTrk->at(i).GetHitT2(0).GetHitY();
    double TrkDX=xp-atrk.GetHitT2(kk).GetHitX();//TestArmXTrk->at(j).GetHitT2(0).GetHitX();
	      
    if((fabs(trkR2-trkR1)<3*atrk.GetHitT2(kk).GetHitDR())&&(fabs(TrkDY)<3*atrk.GetHitT2(kk).GetHitDY())&&(fabs(TrkDX)<3*atrk.GetHitT2(kk).GetHitDX())){
      TrkMatch=true;
    }
  }
  /*
  double trkR1=0.;
    for(unsigned int i=0;i<CleanRefHalfXTrk->size();i++)
      {
	double trkphi1=T2CutUtils.PhiFitAverage(CleanRefHalfXTrk->at(i).GetT2TrackHits());
	double trketa1=T2CutUtils.EtaFromAverages(CleanRefHalfXTrk->at(i));
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
	trkY1=trkY1/((double)sizeref);  trkX1=trkX1/((double)sizeref);  

	double trk1Rmin=0.;
	double trk1Rmax=0.;
	

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

	


	bool cutRefZImpact=false;
	cutRefZImpact=RecoTrk_UnfoldingCut(CleanRefHalfXTrk->at(i));
	if(cutRefZImpact)
	  numreftrksZImpactCut++;



	for(unsigned int j=0;j<TestArmXTrk->size();j++)//TestArmXTrk
	  {

	    unsigned int sizetest=TestArmXTrk->at(j).GetT2TrackHits().size();
	    bool TrkMatch=false;

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
	      
	      
	      if((fabs(trkR2-trkR1)<3*TestArmXTrk->at(j).GetHitT2(kk).GetHitDR())&&(fabs(TrkDY)<3*TestArmXTrk->at(j).GetHitT2(kk).GetHitDY())&&(fabs(TrkDX)<3*TestArmXTrk->at(j).GetHitT2(kk).GetHitDX())){
		TrkMatch=true;
		
	      }
	    }

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

		if(cutRefZImpact)
		  internalcounterZImpactCut++;
		
	      } //If Trk_Match
	    
   
	  }//End-For TestTrack
	
	if(internalcounter>0)
	  trackMatch++;

	if(internalcounterEtaCut>0)
	  trackMatchEtaCut++;

	
	if(internalcounterZImpactCut>0)
	  trackMatchZImpactCut++;
      }
    
    
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
	if(refQuarter==0)
	  H0TrkEffi_AxisAvgMult_ZImpactCut->Fill(avgMult,(trackMatchZImpactCut/numreftrksZImpactCut));
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
	  H0TrkEffi_AxisAvgMult_EtaCut->Fill(avgMult,(trackMatch/numreftrks));	 
	}

	if(refQuarter==1){
	  H1TrkEffi_AxisRawTrk_EtaCut->Fill(numreftrks_raw,(trackMatchEtaCut/numreftrksEtaCut));
	  H1TrkEffi_AxisAvgMult_EtaCut->Fill(avgMult,(trackMatch/numreftrks));

	}
    
	if(refQuarter==2){
	  H2TrkEffi_AxisRawTrk_EtaCut->Fill(numreftrks_raw,(trackMatchEtaCut/numreftrksEtaCut));
	  H2TrkEffi_AxisAvgMult_EtaCut->Fill(avgMult,(trackMatch/numreftrks));
	}

	if(refQuarter==3){
	  H3TrkEffi_AxisRawTrk_EtaCut->Fill(numreftrks_raw,(trackMatchEtaCut/numreftrksEtaCut));
	  H3TrkEffi_AxisAvgMult_EtaCut->Fill(avgMult,(trackMatch/numreftrks));	
	}


      }


      
    }
  */


  double effiret=0.;
  if(TrkMatch)
    effiret=1.;

  return effiret;
  
    //   std::cout<<"CalculateEfficiency-END"<<std::endl;
}


void T2BackgroundAn::Trigger(const T2PadDigiCollection* PadDigiptr, bool & H0trigger, bool & H1trigger, bool & H2trigger, bool & H3trigger)
{
  int Super_Pad_matrix_supr[8][13][40];               //definisco la matrice di piano [row][col][detID] 


  //std::vector<std::vector< std::vector<int> > > Super_Pad_matrix_supr;
  //std::vector<std::vector< std::vector<int> > > Trk_revealed


  
  unsigned long count=0;
  int aux1=0,aux2=0, aux3=0; 
  int Srow=0; 
  int Scol=0; 
  int hit_mf=0, hit_mn=0, hit_pf=0, hit_pn=0; 



  //variable event-type dependent 

	 






  int trg_track_cnt=0; 
  int trg_track_cnt_pf=0; 
  int trg_track_cnt_pn=0; 
  int trg_track_cnt_mf=0; 
  int trg_track_cnt_mn=0; 

  
  
   T2PadDigiCollection::DigiRangeIterator detUnitItP;
  //std::vector<int> howmanypads;



  std::vector<int> vectpadRow;
  std::vector<int> vectpadCol;
  std::vector<int> vectpadId;
  std::vector<int>  Super_Pad_row;    //Creo il vettore che conterr?? le SuperPad row accese 
  std::vector<int>  Super_Pad_column; //Creo il vettore che conterr?? le SuperPad colonne accese 
   
   DigiContainerIterator<T2DetId, T2PadDigi> itp;

   for(itp= PadDigiptr->begin(); itp!=PadDigiptr->end(); ++itp)
       {
	 T2DetId mydet=(*itp).first; 
	 int detnumb = mydet.arm()*20 + mydet.halfTelescope()*10 + mydet.plane()*2 + mydet.planeSide();

	 //Maschero a mano H2
	 if((detnumb!=21)&&(detnumb!=23)&&(detnumb!=25)&&(detnumb!=27)&&(detnumb!=29))
	 for(std::vector<T2PadDigi>::const_iterator itpad =(*itp).second.first; itpad !=(*itp).second.second; itpad++)
	   {
	      vectpadRow.push_back((*itpad).getRow());
	      vectpadCol.push_back((*itpad).getCol());
	      vectpadId.push_back(detnumb);
	      Super_Pad_row.push_back((*itpad).getRow()); 
	      Super_Pad_column.push_back((*itpad).getCol());
	      
	   }

       }

 

  
  { 

    //////////////////////////////////////// 
    // initializations 
    /////////////////////////////////////// 



     							//possono contenere informazioni ridondanti!! 
    //Multy_histo_cut->Fill(the_t2ev->TrkAx.size());//BE		 

    trg_track_cnt=0; 
    trg_track_cnt_pf=0; 
    trg_track_cnt_pn=0; 
    trg_track_cnt_mn=0; 
    trg_track_cnt_mf=0; 


		    
    for (aux1=0; aux1<8; aux1++) for (aux2=0; aux2<13; aux2++) for (aux3=0; aux3<40; aux3++) Super_Pad_matrix_supr[aux1][aux2][aux3]=0; 
    

    //////////////////////////////////////// 
    // Star trigger recontruction 
    /////////////////////////////////////// 
    ////////////////////////////////////// 
    //Remapping degli hit in SuperPad 
    ////////////////////////////////////// 
    
    
    for(count=0; count<vectpadRow.size(); count++)  //Pad_row.size()=Pad_col.size() 
      { 
	//if (the_t2ev->Pad_det.at(count)!=2) continue;  //seleziono solo gli hit sul piano 2 
	
	Srow= (int) vectpadRow.at(count)/3; 
	Scol= (int) vectpadCol.at(count)/5; 
	
	//if (the_t2ev->Pad_row.at(count)==0) printf("found 0"); 
	Super_Pad_row[count]= Srow; 
	Super_Pad_column[count]= Scol; 
	Super_Pad_matrix_supr[Srow][Scol][vectpadId.at(count)]=1; 
      } 
    
    ///////////////////////////////////////////////////////////// 
    //Ricostruzione dell'evento dal punto di vista del trigger. 
    // criterio di coincidenza 
    ///////////////////////////////////////////////////////////// 

    
    
    //printf("Super_Pad_matrix\n"); 
    for (aux1=0; aux1<8; aux1++) 
      { 
	//printf("\n"); 
	for (aux2=0; aux2<13; aux2++) 
	  { 
	    hit_mf=0; 
	    hit_mn=0; 
	    hit_pf=0; 
	    hit_pn=0; 
	    //for (aux3=0; aux3<40; aux3++) 
	    for (aux3=0; aux3<40; aux3+=2) 
	      { 
					  
		if (Super_Pad_matrix_supr[aux1][aux2][aux3]==1) 
		  { 
		    if (aux3<10) hit_pn++; 
		    else if (aux3<20) hit_pf++; 
		    else if (aux3<30) hit_mn++; 
		    else if (aux3<40) hit_mf++; 
		    else printf("Ma che succede?\n"); 
		  } 
	      } 
	    for (aux3=1; aux3<40; aux3+=2) 
	      { 
		if (Super_Pad_matrix_supr[aux1][12-aux2][aux3]==1) //odd plane have reverse column id 
		  { 
		    if (aux3<10) hit_pn++; 
		    else if (aux3<20) hit_pf++; 
		    else if (aux3<30) hit_mn++; 
		    else if (aux3<40) hit_mf++; 
		    else printf("Ma che succede?\n"); 
		  } 
	      } 
	    
	    
	    if (hit_mn>10 | hit_mf>10 | hit_pn>10 | hit_pf>10 ) printf("WARNING:TO MANY HIT!"); 
	    if (hit_mn>4)    //check if the trigger road have been detected by 4 or more planes 
	      { 
		

		trg_track_cnt++; 
		trg_track_cnt_mn++; 
		
	      } 
	    if (hit_mf>4) 
	      { 
		

		trg_track_cnt++; 
		trg_track_cnt_mf++; 
		
	      } 
	    if (hit_pn>4) 
	      { 
		

		trg_track_cnt++; 
		trg_track_cnt_pn++; 
	      } 
	    if (hit_pf>4) 
	      { 
		

		trg_track_cnt++; 
		trg_track_cnt_pf++; 
	      }	 
	  } 
	
      }  
  }
  
  if(trg_track_cnt_pf>0)
    H1trigger=true;
  else
    H1trigger=false;
  
  if(trg_track_cnt_pn>0)
    H0trigger=true;
  else
    H0trigger=false;

  if(trg_track_cnt_mf>0)
    H3trigger=true;
  else
    H3trigger=false;
  
  if(trg_track_cnt_mn>0)
    H2trigger=true;
  else
    H2trigger=false;

}









void T2BackgroundAn::ZImpactRangeFromFit(int quarter, double eta2, double & ZImpLeft, double & ZImpRight){

  double eta3=fabs(eta2);
  unsigned int index=etaToPos(eta3);

  if(quarter==0){    
    if(index>=ZimpCutat95DataVH0.size())
      std::cout<<"Error inside ZImpactRangeFromFit"<<std::endl;
    else{
      ZImpLeft=ZimpCutat95DataLVH0.at(index);
      ZImpRight=ZimpCutat95DataRVH0.at(index);
    }
  }

  if(quarter==1){
     if(index>=ZimpCutat95DataVH1.size())
      std::cout<<"Error inside ZImpactRangeFromFit"<<std::endl;
    else{
      ZImpLeft=ZimpCutat95DataLVH1.at(index);
      ZImpRight=ZimpCutat95DataRVH1.at(index);
    }
  }
  if(quarter==2){
     if(index>=ZimpCutat95DataVH2.size())
      std::cout<<"Error inside ZImpactRangeFromFit"<<std::endl;
    else{
      ZImpLeft=ZimpCutat95DataLVH2.at(index);
      ZImpRight=ZimpCutat95DataRVH2.at(index);
    }
  }
  if(quarter==3){
     if(index>=ZimpCutat95DataVH3.size())
      std::cout<<"Error inside ZImpactRangeFromFit"<<std::endl;
    else{
      ZImpLeft=ZimpCutat95DataLVH3.at(index);
      ZImpRight=ZimpCutat95DataRVH3.at(index);
    }
  }
  
}



unsigned int  T2BackgroundAn::etaToPos(double centre){
 
  unsigned int pos=0;
  centre=fabs(centre); double deltaMezzi=0.025;
  
  if(centre>5.375-deltaMezzi) 
    if(centre<=5.375+deltaMezzi) 
      pos=0;  

  if(centre>5.425-deltaMezzi) 
    if(centre<=5.425+deltaMezzi) 
      pos=1;


  if(centre>5.625-deltaMezzi) 
    if(centre<=5.625+deltaMezzi) 
      pos=2;
 
  if(centre>5.675-deltaMezzi) 
    if(centre<=5.675+deltaMezzi) 
      pos=3;
  
  if(centre>5.725-deltaMezzi) 
    if(centre<=5.725+deltaMezzi) 
      pos=4;
  
  if(centre>5.775-deltaMezzi) 
    if(centre<=5.775+deltaMezzi) 
      pos=5;
  
  if(centre>5.825-deltaMezzi) 
    if(centre<=5.825+deltaMezzi) 
      pos=6;
  
  if(centre>5.875-deltaMezzi) 
    if(centre<=5.875+deltaMezzi) 
      pos=7;
  
  if(centre>5.925-deltaMezzi) 
    if(centre<=5.925+deltaMezzi) 
      pos=8;
  
  if(centre>5.975-deltaMezzi) 
    if(centre<=5.975+deltaMezzi) 
      pos=9;
  
  if(centre>6.025-deltaMezzi) 
    if(centre<=6.025+deltaMezzi) 
      pos=10;
  
  if(centre>6.075-deltaMezzi) 
    if(centre<=6.075+deltaMezzi) 
      pos=11;
 
  if(centre>6.125-deltaMezzi) 
    if(centre<=6.125+deltaMezzi) 
      pos=12;
 
  if(centre>6.175-deltaMezzi) 
    if(centre<=6.175+deltaMezzi) 
      pos=13;
 
  if(centre>6.225-deltaMezzi) 
    if(centre<=6.225+deltaMezzi) 
      pos=14;
 
  if(centre>6.275-deltaMezzi) 
    if(centre<=6.275+deltaMezzi) 
      pos=15;
 
  if(centre>6.325-deltaMezzi) 
    if(centre<=6.325+deltaMezzi) 
      pos=16;
 
  if(centre>6.375-deltaMezzi) 
    if(centre<=6.375+deltaMezzi) 
      pos=17;
 
  if(centre>6.375+deltaMezzi)
    pos=17;

  return pos;
}

// ------------ method called once each job just before starting event loop  ------------
void 
T2BackgroundAn::beginJob()//const edm::EventSetup&
{
  if(verbosity>0)
    std::cout<<"BeginJon"<<std::endl;
  numevent=0;

  numshowevent=0;
  numtotevent=0;
  numcleanevent=0;


  K0EventCounter=0;
  HadrEventCounter=0;
  knowPart.push_back(11); knowPart.push_back(12);  knowPart.push_back(13);  knowPart.push_back(14); knowPart.push_back(15); 
  knowPart.push_back(16); knowPart.push_back(17); knowPart.push_back(18); knowPart.push_back(22); knowPart.push_back(23);
  knowPart.push_back(24); 

  knowPart.push_back(130); knowPart.push_back(310); knowPart.push_back(311); 
  
  knowPart.push_back(111); knowPart.push_back(211); knowPart.push_back(113); knowPart.push_back(213); knowPart.push_back(215); 
  knowPart.push_back(117); knowPart.push_back(217); knowPart.push_back(119); knowPart.push_back(219); knowPart.push_back(221); 
  knowPart.push_back(331); knowPart.push_back(223); knowPart.push_back(333); knowPart.push_back(225); knowPart.push_back(335); 
  knowPart.push_back(227); 
  knowPart.push_back(337); knowPart.push_back(229);  knowPart.push_back(411); knowPart.push_back(421); knowPart.push_back(413); 
  knowPart.push_back(423); knowPart.push_back(415); knowPart.push_back(425); knowPart.push_back(431); knowPart.push_back(433); 
  knowPart.push_back(435); knowPart.push_back(511); knowPart.push_back(521); knowPart.push_back(513); knowPart.push_back(523); 
  knowPart.push_back(515); knowPart.push_back(525); knowPart.push_back(531); knowPart.push_back(533); knowPart.push_back(535); 
  knowPart.push_back(541); knowPart.push_back(543); knowPart.push_back(545); knowPart.push_back(441); knowPart.push_back(443); 
  knowPart.push_back(445); knowPart.push_back(551); knowPart.push_back(553); knowPart.push_back(555); knowPart.push_back(557); 
  knowPart.push_back(2212); knowPart.push_back(2112); knowPart.push_back(2224); knowPart.push_back(2214); knowPart.push_back(2114); 
  knowPart.push_back(1114); knowPart.push_back(3122); knowPart.push_back(3222); knowPart.push_back(3212); 
  knowPart.push_back(3112); knowPart.push_back(3224); knowPart.push_back(3214); knowPart.push_back(3114); knowPart.push_back(3322); 
  knowPart.push_back(3312); knowPart.push_back(3324); knowPart.push_back(3314); knowPart.push_back(3334); knowPart.push_back(4122); 
  knowPart.push_back(4222); knowPart.push_back(4212); knowPart.push_back(4112); knowPart.push_back(4224); knowPart.push_back(4214); 
  knowPart.push_back(4114); knowPart.push_back(4232); knowPart.push_back(4132); knowPart.push_back(4322); knowPart.push_back(4312); 
  knowPart.push_back(4324); knowPart.push_back(4314); knowPart.push_back(4332); knowPart.push_back(4334); knowPart.push_back(4412); 
  knowPart.push_back(4422); knowPart.push_back(4414); knowPart.push_back(4424); knowPart.push_back(4432); knowPart.push_back(4434); 
  knowPart.push_back(4444); knowPart.push_back(5122); knowPart.push_back(5112); knowPart.push_back(5212); knowPart.push_back(5222); 
  knowPart.push_back(5114); knowPart.push_back(5214); knowPart.push_back(5224); knowPart.push_back(5132); knowPart.push_back(5232); 
  knowPart.push_back(5312); knowPart.push_back(5322); knowPart.push_back(5314); knowPart.push_back(5324); knowPart.push_back(5332); 
  knowPart.push_back(5334); knowPart.push_back(5142); knowPart.push_back(5242); knowPart.push_back(5412); knowPart.push_back(5422); 
  knowPart.push_back(5414); knowPart.push_back(5424); knowPart.push_back(5342); knowPart.push_back(5432); knowPart.push_back(5434); 
  knowPart.push_back(5442); knowPart.push_back(5444); knowPart.push_back(5512); knowPart.push_back(5522); knowPart.push_back(5514); 
  knowPart.push_back(5524); knowPart.push_back(5532); knowPart.push_back(5534); knowPart.push_back(5542); knowPart.push_back(5544); 
  knowPart.push_back(5554);



  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  
  
  //AllowedDecayToCount.push_back(11); AllowedDecayToCount.push_back(12); electr and it neutrino Excluded (not decayng particle)
  //AllowedDecayToCount.push_back(13);  AllowedDecayToCount.push_back(14); muon ans its neutrino Excluded (lifetime too big ~10-6 to decay before T2)

  //WARNING: 66% of Tau decay involve hadrons in the final states.
  AllowedDecayToCount.push_back(15); //10-12s 10-13s (Tau), will decay before T2.
  //AllowedDecayToCount.push_back(16); AllowedDecayToCount.push_back(17); AllowedDecayToCount.push_back(18); Tau and its neutrino excluded

  AllowedDecayToCount.push_back(23); //W 10-25s
  AllowedDecayToCount.push_back(24); //Z 

  AllowedDecayToCount.push_back(130); AllowedDecayToCount.push_back(310); AllowedDecayToCount.push_back(311); 
  //K0L, K0S,K0. 10-8 and 10-11 sec can be seen in T2. Have pi0 and charged decay mode: To review
   //K+ lifetime 10-8
  
  //AllowedDecayToCount.push_back(111); Pi0 Excluded
  
  
  AllowedDecayToCount.push_back(211); //Pi+- lifetime is 10-8 but gamma~10 for E= 1GeV Can be seen in T2.
  
  AllowedDecayToCount.push_back(113); //pi+ + pi- decay
  
  //WARNING: to rewievw: pi+ + pi- and 2pi0 decay
  AllowedDecayToCount.push_back(213);  //rho+
  
  AllowedDecayToCount.push_back(215); //a2 meson ?


 
  //Other quite rare I=1 mesons
  AllowedDecayToCount.push_back(117); AllowedDecayToCount.push_back(217); AllowedDecayToCount.push_back(119); AllowedDecayToCount.push_back(219); 

  //WARNING: Neutral decay mode is 72%. At the moment I drop it. Very short lifetime
  //AllowedDecayToCount.push_back(221); AllowedDecayToCount.push_back(331);//eta, eta prime
   

  AllowedDecayToCount.push_back(223); AllowedDecayToCount.push_back(333); //w phi
  AllowedDecayToCount.push_back(225); AllowedDecayToCount.push_back(335); //f2 f2'
  AllowedDecayToCount.push_back(227); //w3 phi3
  AllowedDecayToCount.push_back(337); AllowedDecayToCount.push_back(229);  
  AllowedDecayToCount.push_back(411); AllowedDecayToCount.push_back(421); AllowedDecayToCount.push_back(413); 
  AllowedDecayToCount.push_back(423); AllowedDecayToCount.push_back(415); AllowedDecayToCount.push_back(425); AllowedDecayToCount.push_back(431); AllowedDecayToCount.push_back(433); 
  AllowedDecayToCount.push_back(435); AllowedDecayToCount.push_back(511); AllowedDecayToCount.push_back(521); AllowedDecayToCount.push_back(513); AllowedDecayToCount.push_back(523); 
  AllowedDecayToCount.push_back(515); AllowedDecayToCount.push_back(525); AllowedDecayToCount.push_back(531); AllowedDecayToCount.push_back(533); AllowedDecayToCount.push_back(535); 
  AllowedDecayToCount.push_back(541); AllowedDecayToCount.push_back(543); AllowedDecayToCount.push_back(545); AllowedDecayToCount.push_back(441); AllowedDecayToCount.push_back(443); 
  AllowedDecayToCount.push_back(445); AllowedDecayToCount.push_back(551); AllowedDecayToCount.push_back(553); AllowedDecayToCount.push_back(555); AllowedDecayToCount.push_back(557); 
  
  //AllowedDecayToCount.push_back(2212); AllowedDecayToCount.push_back(2112); //p and n excluded

  AllowedDecayToCount.push_back(2224); AllowedDecayToCount.push_back(2214); AllowedDecayToCount.push_back(2114); 
  AllowedDecayToCount.push_back(1114); AllowedDecayToCount.push_back(3122); AllowedDecayToCount.push_back(3222); AllowedDecayToCount.push_back(3212); 
  AllowedDecayToCount.push_back(3112); AllowedDecayToCount.push_back(3224); AllowedDecayToCount.push_back(3214); AllowedDecayToCount.push_back(3114); AllowedDecayToCount.push_back(3322); 
  AllowedDecayToCount.push_back(3312); AllowedDecayToCount.push_back(3324); AllowedDecayToCount.push_back(3314); AllowedDecayToCount.push_back(3334); AllowedDecayToCount.push_back(4122); 
  AllowedDecayToCount.push_back(4222); AllowedDecayToCount.push_back(4212); AllowedDecayToCount.push_back(4112); AllowedDecayToCount.push_back(4224); AllowedDecayToCount.push_back(4214); 
  AllowedDecayToCount.push_back(4114); AllowedDecayToCount.push_back(4232); AllowedDecayToCount.push_back(4132); AllowedDecayToCount.push_back(4322); AllowedDecayToCount.push_back(4312); 
  AllowedDecayToCount.push_back(4324); AllowedDecayToCount.push_back(4314); AllowedDecayToCount.push_back(4332); AllowedDecayToCount.push_back(4334); AllowedDecayToCount.push_back(4412); 
  AllowedDecayToCount.push_back(4422); AllowedDecayToCount.push_back(4414); AllowedDecayToCount.push_back(4424); AllowedDecayToCount.push_back(4432); AllowedDecayToCount.push_back(4434); 
  AllowedDecayToCount.push_back(4444); AllowedDecayToCount.push_back(5122); AllowedDecayToCount.push_back(5112); AllowedDecayToCount.push_back(5212); AllowedDecayToCount.push_back(5222); 
  AllowedDecayToCount.push_back(5114); AllowedDecayToCount.push_back(5214); AllowedDecayToCount.push_back(5224); AllowedDecayToCount.push_back(5132); AllowedDecayToCount.push_back(5232); 
  AllowedDecayToCount.push_back(5312); AllowedDecayToCount.push_back(5322); AllowedDecayToCount.push_back(5314); AllowedDecayToCount.push_back(5324); AllowedDecayToCount.push_back(5332); 
  AllowedDecayToCount.push_back(5334); AllowedDecayToCount.push_back(5142); AllowedDecayToCount.push_back(5242); AllowedDecayToCount.push_back(5412); AllowedDecayToCount.push_back(5422); 
  AllowedDecayToCount.push_back(5414); AllowedDecayToCount.push_back(5424); AllowedDecayToCount.push_back(5342); AllowedDecayToCount.push_back(5432); AllowedDecayToCount.push_back(5434); 
  AllowedDecayToCount.push_back(5442); AllowedDecayToCount.push_back(5444); AllowedDecayToCount.push_back(5512); AllowedDecayToCount.push_back(5522); AllowedDecayToCount.push_back(5514); 
  AllowedDecayToCount.push_back(5524); AllowedDecayToCount.push_back(5532); AllowedDecayToCount.push_back(5534); AllowedDecayToCount.push_back(5542); AllowedDecayToCount.push_back(5544); 
  AllowedDecayToCount.push_back(5554);





  std::cout<<"Loading All ZImpact ranges"<<std::endl;


  TFile *ffunct6325H0 = 0;
  TFile *ffunct6375H0 = 0;
  TFile *ffunct6225H0 = 0;
  TFile *ffunct6275H0 = 0;
  TFile *ffunct6125H0 = 0;
  TFile *ffunct6175H0 = 0;
  TFile *ffunct6025H0 = 0;
  TFile *ffunct6075H0 = 0;
  TFile *ffunct5925H0 = 0;
  TFile *ffunct5975H0 = 0;
  TFile *ffunct5825H0 = 0;
  TFile *ffunct5875H0 = 0;
  TFile *ffunct5775H0 = 0;
  TFile *ffunct5725H0 = 0;
  TFile *ffunct5675H0 = 0;
  TFile *ffunct5625H0 = 0;
  TFile *ffunct5375H0 = 0;
  TFile *ffunct5425H0 = 0;
 
  TFile *ffunct6325H1 = 0;
  TFile *ffunct6375H1 = 0;
  TFile *ffunct6225H1 = 0;
  TFile *ffunct6275H1 = 0;
  TFile *ffunct6125H1 = 0;
  TFile *ffunct6175H1 = 0;
  TFile *ffunct6025H1 = 0;
  TFile *ffunct6075H1 = 0;
  TFile *ffunct5925H1 = 0;
  TFile *ffunct5975H1 = 0;
  TFile *ffunct5825H1 = 0;
  TFile *ffunct5875H1 = 0;
  TFile *ffunct5775H1 = 0;
  TFile *ffunct5725H1 = 0;
  TFile *ffunct5675H1 = 0;
  TFile *ffunct5625H1 = 0;
  TFile *ffunct5375H1 = 0;
  TFile *ffunct5425H1 = 0;
 

  TFile *ffunct6325H2 = 0;
  TFile *ffunct6375H2 = 0;
  TFile *ffunct6225H2 = 0;
  TFile *ffunct6275H2 = 0;
  TFile *ffunct6125H2 = 0;
  TFile *ffunct6175H2 = 0;
  TFile *ffunct6025H2 = 0;
  TFile *ffunct6075H2 = 0;
  TFile *ffunct5925H2 = 0;
  TFile *ffunct5975H2 = 0;
  TFile *ffunct5825H2 = 0;
  TFile *ffunct5875H2 = 0;
  TFile *ffunct5775H2 = 0;
  TFile *ffunct5725H2 = 0;
  TFile *ffunct5675H2 = 0;
  TFile *ffunct5625H2 = 0;
  TFile *ffunct5375H2 = 0;
  TFile *ffunct5425H2 = 0;
 
  TFile *ffunct6325H3 = 0;
  TFile *ffunct6375H3 = 0;
  TFile *ffunct6225H3 = 0;
  TFile *ffunct6275H3 = 0;
  TFile *ffunct6125H3 = 0;
  TFile *ffunct6175H3 = 0;
  TFile *ffunct6025H3 = 0;
  TFile *ffunct6075H3 = 0;
  TFile *ffunct5925H3 = 0;
  TFile *ffunct5975H3 = 0;
  TFile *ffunct5825H3 = 0;
  TFile *ffunct5875H3 = 0;
  TFile *ffunct5775H3 = 0;
  TFile *ffunct5725H3 = 0;
  TFile *ffunct5675H3 = 0;
  TFile *ffunct5625H3 = 0;
  TFile *ffunct5375H3 = 0;
  TFile *ffunct5425H3 = 0;
 

  bool initialized = false;
  //OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6175_P8Simulation8TeV.root _MaxClu45
  if((MaxPadCluOfFittingCurves==65)||(MaxPadCluOfFittingCurves==65)){

    if(NameOfGenerator=="Py8"){
    initialized = true;
      std::cout<<"Loading of H0 Fit Curves 65..."<<std::endl;

      ffunct6325H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6325_P8Simulation8TeV.root");
      ffunct6375H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6375_P8Simulation8TeV.root");

      ffunct6225H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6225_P8Simulation8TeV.root");
      ffunct6275H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6275_P8Simulation8TeV.root");

      ffunct6125H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6125_P8Simulation8TeV.root");
      ffunct6175H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6175_P8Simulation8TeV.root");
      
      std::cout<<" ..."<<std::endl;
      
      ffunct6025H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6025_P8Simulation8TeV.root");
      ffunct6075H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6075_P8Simulation8TeV.root");

      ffunct5925H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5925_P8Simulation8TeV.root");
      ffunct5975H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5975_P8Simulation8TeV.root");
      
      std::cout<<" ..."<<std::endl;

      ffunct5825H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5825_P8Simulation8TeV.root");
      ffunct5875H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5875_P8Simulation8TeV.root");
      
      std::cout<<" ..."<<std::endl;
      
      ffunct5775H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5775_P8Simulation8TeV.root");
      ffunct5725H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5725_P8Simulation8TeV.root");

      ffunct5675H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5675_P8Simulation8TeV.root");
      ffunct5625H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5625_P8Simulation8TeV.root");

      ffunct5375H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5375_P8Simulation8TeV.root");
      ffunct5425H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5425_P8Simulation8TeV.root");
  
      std::cout<<"Loading of H1 Fit Curves..."<<std::endl;
      
 
      
      ffunct6325H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6325_P8Simulation8TeV.root");
      ffunct6375H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6375_P8Simulation8TeV.root");

      ffunct6225H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6225_P8Simulation8TeV.root");
      ffunct6275H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6275_P8Simulation8TeV.root");

      ffunct6125H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6125_P8Simulation8TeV.root");
      ffunct6175H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6175_P8Simulation8TeV.root");

      ffunct6025H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6025_P8Simulation8TeV.root");
      ffunct6075H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6075_P8Simulation8TeV.root");

      ffunct5925H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5925_P8Simulation8TeV.root");
      ffunct5975H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5975_P8Simulation8TeV.root");

      ffunct5825H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5825_P8Simulation8TeV.root");
      ffunct5875H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5875_P8Simulation8TeV.root");

      ffunct5775H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5775_P8Simulation8TeV.root");
      ffunct5725H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5725_P8Simulation8TeV.root");

      ffunct5675H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5675_P8Simulation8TeV.root");
      ffunct5625H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5625_P8Simulation8TeV.root");

      ffunct5375H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5375_P8Simulation8TeV.root");
      ffunct5425H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5425_P8Simulation8TeV.root");
  
      



      std::cout<<"Loading of H2 Fit Curves..."<<std::endl;

    
      ffunct6325H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6325_P8Simulation8TeV.root");
      ffunct6375H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6375_P8Simulation8TeV.root");

      ffunct6225H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6225_P8Simulation8TeV.root");
      ffunct6275H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6275_P8Simulation8TeV.root");

      ffunct6125H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6125_P8Simulation8TeV.root");
      ffunct6175H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6175_P8Simulation8TeV.root");

      ffunct6025H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6025_P8Simulation8TeV.root");
      ffunct6075H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6075_P8Simulation8TeV.root");

      ffunct5925H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5925_P8Simulation8TeV.root");
      ffunct5975H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5975_P8Simulation8TeV.root");

      ffunct5825H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5825_P8Simulation8TeV.root");
      ffunct5875H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5875_P8Simulation8TeV.root");

      ffunct5775H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5775_P8Simulation8TeV.root");
      ffunct5725H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5725_P8Simulation8TeV.root");

      ffunct5675H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5675_P8Simulation8TeV.root");
      ffunct5625H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5625_P8Simulation8TeV.root");

      ffunct5375H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5375_P8Simulation8TeV.root");
      ffunct5425H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5425_P8Simulation8TeV.root");
    


      std::cout<<"Loading of H3 Fit Curves..."<<std::endl;      
    
      ffunct6325H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6325_P8Simulation8TeV.root");
      ffunct6375H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6375_P8Simulation8TeV.root");

      ffunct6225H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6225_P8Simulation8TeV.root");
      ffunct6275H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6275_P8Simulation8TeV.root");

      ffunct6125H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6125_P8Simulation8TeV.root");
      ffunct6175H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6175_P8Simulation8TeV.root");

      ffunct6025H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6025_P8Simulation8TeV.root");
      ffunct6075H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6075_P8Simulation8TeV.root");

      ffunct5925H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5925_P8Simulation8TeV.root");
      ffunct5975H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5975_P8Simulation8TeV.root");

      ffunct5825H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5825_P8Simulation8TeV.root");
      ffunct5875H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5875_P8Simulation8TeV.root");

      ffunct5775H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5775_P8Simulation8TeV.root");
      ffunct5725H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5725_P8Simulation8TeV.root");

      ffunct5675H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5675_P8Simulation8TeV.root");
      ffunct5625H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5625_P8Simulation8TeV.root");

      ffunct5375H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5375_P8Simulation8TeV.root");
      ffunct5425H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5425_P8Simulation8TeV.root");
    }
  }



  

  if(MaxPadCluOfFittingCurves==45){
	  initialized = true;
      std::cout<<"Loading of H0 Fit Curves 45..."<<std::endl;
      if(NameOfGenerator=="Py8"){

      ffunct6325H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6325_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6375H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6375_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6225H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6225_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6275H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6275_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      std::cout<<"..."<<std::endl;

      ffunct6125H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6125_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6175H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6175_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6025H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6025_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6075H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6075_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      std::cout<<"..."<<std::endl;
      
      ffunct5925H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5925_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5975H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5975_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5825H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5825_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5875H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5875_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      std::cout<<"..."<<std::endl;

      ffunct5775H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5775_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5725H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5725_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5675H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5675_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5625H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5625_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5375H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5375_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5425H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5425_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
  
      std::cout<<"Loading of H1 Fit Curves..."<<std::endl;
      
 
      
      ffunct6325H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6325_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6375H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6375_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6225H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6225_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6275H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6275_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6125H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6125_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6175H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6175_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6025H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6025_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6075H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6075_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5925H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5925_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5975H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5975_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5825H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5825_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5875H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5875_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5775H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5775_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5725H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5725_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5675H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5675_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5625H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5625_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5375H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5375_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5425H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5425_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
  
      



      std::cout<<"Loading of H2 Fit Curves..."<<std::endl;

    
      ffunct6325H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6325_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6375H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6375_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6225H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6225_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6275H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6275_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6125H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6125_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6175H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6175_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6025H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6025_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6075H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6075_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5925H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5925_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5975H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5975_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5825H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5825_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5875H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5875_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5775H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5775_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5725H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5725_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5675H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5675_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5625H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5625_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5375H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5375_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5425H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5425_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
    


      std::cout<<"Loading of H3 Fit Curves..."<<std::endl;      
    
      ffunct6325H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6325_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6375H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6375_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6225H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6225_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6275H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6275_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6125H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6125_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6175H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6175_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6025H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6025_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6075H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6075_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5925H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5925_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5975H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5975_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5825H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5825_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5875H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5875_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5775H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5775_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5725H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5725_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5675H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5675_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5625H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5625_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5375H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5375_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5425H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5425_MCPy8_424_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      }

      
      if(NameOfGenerator=="EPOS"){
	//       OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5875_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root
      std::cout<<"Loading of H0 Fit Curves EPOS-45..."<<std::endl;      
    
      ffunct6325H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6325_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6375H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6375_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6225H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6225_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6275H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6275_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      std::cout<<"..."<<std::endl;

      ffunct6125H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6125_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6175H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6175_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6025H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6025_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6075H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin6075_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      std::cout<<"..."<<std::endl;
      
      ffunct5925H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5925_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5975H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5975_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5825H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5825_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5875H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5875_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      std::cout<<"..."<<std::endl;

      ffunct5775H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5775_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5725H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5725_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5675H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5675_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5625H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5625_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5375H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5375_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5425H0=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H0_ZminCut13500_All_Bin5425_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
  
      std::cout<<"Loading of H1 Fit Curves..."<<std::endl;
      
 
      
      ffunct6325H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6325_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6375H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6375_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6225H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6225_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6275H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6275_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6125H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6125_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6175H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6175_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6025H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6025_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6075H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin6075_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5925H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5925_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5975H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5975_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5825H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5825_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5875H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5875_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5775H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5775_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5725H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5725_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5675H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5675_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5625H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5625_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5375H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5375_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5425H1=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H1_ZminCut13500_All_Bin5425_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
  
      



      std::cout<<"Loading of H2 Fit Curves..."<<std::endl;

    
      ffunct6325H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6325_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6375H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6375_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6225H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6225_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6275H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6275_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6125H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6125_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6175H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6175_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6025H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6025_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6075H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin6075_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5925H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5925_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5975H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5975_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5825H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5825_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5875H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5875_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5775H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5775_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5725H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5725_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5675H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5675_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5625H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5625_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5375H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5375_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5425H2=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H2_ZminCut13500_All_Bin5425_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
    


      std::cout<<"Loading of H3 Fit Curves..."<<std::endl;      
    
      ffunct6325H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6325_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6375H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6375_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6225H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6225_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6275H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6275_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6125H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6125_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6175H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6175_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct6025H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6025_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct6075H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin6075_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5925H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5925_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5975H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5975_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5825H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5825_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5875H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5875_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5775H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5775_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5725H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5725_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5675H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5675_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5625H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5625_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");

      ffunct5375H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5375_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      ffunct5425H3=new TFile("/afs/cern.ch/exp/totem/scratch/berretti/Test_424_Realese/CMSSW_4_2_4/src/TotemT1T2Validation/T2BackgroundAnalysis/test/FitOutput_MaxClu45_EPOS/OutpFitsTrkZImpact_H3_ZminCut13500_All_Bin5425_MCEPOS_445_Tun8372__Plot_ZImpCurvemerged_MaxClu45.root");
      }
  }

	if(!initialized){
		 std::cout<<"Fit Curves NOT LOADED!"<<std::endl;
		 return;
	}else{
		  std::cout<<"Fit Curves LOADED!"<<std::endl;
	}

   
         
      //   etaToPos
      etacentreBinV.clear();
      ZimpCutat95DataVH0.clear();ZimpCutat95DataLVH0.clear();ZimpCutat95DataRVH0.clear();
      ZimpCutat95DataVH1.clear();ZimpCutat95DataLVH1.clear();ZimpCutat95DataRVH1.clear();
      ZimpCutat95DataVH2.clear();ZimpCutat95DataLVH2.clear();ZimpCutat95DataRVH2.clear();
      ZimpCutat95DataVH3.clear();ZimpCutat95DataLVH3.clear();ZimpCutat95DataRVH3.clear();
      double ZimpCutat95Data=0.;double ZimpCutat95DataR=0.;double ZimpCutat95DataL=0.;

      
      etacentreBinV.push_back(5.375);etacentreBinV.push_back(5.425);etacentreBinV.push_back(5.625);etacentreBinV.push_back(5.675);
      etacentreBinV.push_back(5.725);etacentreBinV.push_back(5.775);etacentreBinV.push_back(5.825);etacentreBinV.push_back(5.875);
      etacentreBinV.push_back(5.925);etacentreBinV.push_back(5.975);etacentreBinV.push_back(6.025);etacentreBinV.push_back(6.075);
      etacentreBinV.push_back(6.125);etacentreBinV.push_back(6.175);etacentreBinV.push_back(6.225);etacentreBinV.push_back(6.275);
      etacentreBinV.push_back(6.325);etacentreBinV.push_back(6.375);


      ///////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////
      ZimpCutat95Data=((TH1F*)ffunct5375H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5375H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5375H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct5425H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5425H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5425H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct5625H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5625H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5625H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);
       
      ZimpCutat95Data=((TH1F*)ffunct5675H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5675H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5675H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);
      
      ZimpCutat95Data=((TH1F*)ffunct5725H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5725H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5725H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct5775H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5775H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5775H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);
     
      ZimpCutat95Data=((TH1F*)ffunct5825H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5825H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5825H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;  
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);
      
      ZimpCutat95Data=((TH1F*)ffunct5875H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5875H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5875H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);
       
      ZimpCutat95Data=((TH1F*)ffunct5925H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5925H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5925H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);
      
      ZimpCutat95Data=((TH1F*)ffunct5975H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5975H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5975H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6025H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6025H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6025H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6075H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6075H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6075H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6125H0->Get("ZimpCutat95H"))->GetMean();    
      ZimpCutat95DataL=((TH1F*)ffunct6125H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6125H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6175H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6175H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6175H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6225H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6225H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6225H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;  
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6275H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6275H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6275H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6325H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6325H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6325H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6375H0->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6375H0->Get("MuofDG"))->GetMean()-ZimpCutat95Data; 
      ZimpCutat95DataR=((TH1F*)ffunct6375H0->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH0.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH0.push_back(ZimpCutat95DataR);ZimpCutat95DataVH0.push_back(ZimpCutat95Data);
      
       ///////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////
      ZimpCutat95Data=((TH1F*)ffunct5375H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5375H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5375H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct5425H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5425H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5425H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct5625H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5625H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5625H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);
       
      ZimpCutat95Data=((TH1F*)ffunct5675H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5675H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5675H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);
      
      ZimpCutat95Data=((TH1F*)ffunct5725H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5725H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5725H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct5775H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5775H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5775H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);
     
      ZimpCutat95Data=((TH1F*)ffunct5825H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5825H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5825H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;  
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);
      
      ZimpCutat95Data=((TH1F*)ffunct5875H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5875H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5875H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);
       
      ZimpCutat95Data=((TH1F*)ffunct5925H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5925H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5925H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);
      
      ZimpCutat95Data=((TH1F*)ffunct5975H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5975H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5975H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6025H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6025H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6025H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6075H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6075H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6075H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6125H1->Get("ZimpCutat95H"))->GetMean();    
      ZimpCutat95DataL=((TH1F*)ffunct6125H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6125H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6175H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6175H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6175H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6225H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6225H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6225H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;  
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6275H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6275H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6275H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6325H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6325H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6325H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6375H1->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6375H1->Get("MuofDG"))->GetMean()-ZimpCutat95Data; 
      ZimpCutat95DataR=((TH1F*)ffunct6375H1->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH1.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH1.push_back(ZimpCutat95DataR);ZimpCutat95DataVH1.push_back(ZimpCutat95Data);
       
     
       ///////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////
      ZimpCutat95Data=((TH1F*)ffunct5375H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5375H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5375H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct5425H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5425H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5425H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct5625H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5625H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5625H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);
       
      ZimpCutat95Data=((TH1F*)ffunct5675H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5675H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5675H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);
      
      ZimpCutat95Data=((TH1F*)ffunct5725H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5725H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5725H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct5775H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5775H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5775H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);
     
      ZimpCutat95Data=((TH1F*)ffunct5825H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5825H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5825H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;  
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);
      
      ZimpCutat95Data=((TH1F*)ffunct5875H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5875H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5875H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);
       
      ZimpCutat95Data=((TH1F*)ffunct5925H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5925H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5925H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);
      
      ZimpCutat95Data=((TH1F*)ffunct5975H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5975H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5975H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6025H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6025H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6025H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6075H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6075H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6075H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6125H2->Get("ZimpCutat95H"))->GetMean();    
      ZimpCutat95DataL=((TH1F*)ffunct6125H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6125H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6175H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6175H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6175H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6225H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6225H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6225H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;  
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6275H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6275H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6275H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6325H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6325H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6325H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6375H2->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6375H2->Get("MuofDG"))->GetMean()-ZimpCutat95Data; 
      ZimpCutat95DataR=((TH1F*)ffunct6375H2->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH2.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH2.push_back(ZimpCutat95DataR);ZimpCutat95DataVH2.push_back(ZimpCutat95Data);



       ///////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////
      ZimpCutat95Data=((TH1F*)ffunct5375H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5375H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5375H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct5425H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5425H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5425H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct5625H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5625H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5625H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);
       
      ZimpCutat95Data=((TH1F*)ffunct5675H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5675H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5675H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);
      
      ZimpCutat95Data=((TH1F*)ffunct5725H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5725H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5725H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct5775H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5775H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5775H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);
     
      ZimpCutat95Data=((TH1F*)ffunct5825H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5825H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5825H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;  
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);
      
      ZimpCutat95Data=((TH1F*)ffunct5875H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5875H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5875H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);
       
      ZimpCutat95Data=((TH1F*)ffunct5925H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5925H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5925H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);
      
      ZimpCutat95Data=((TH1F*)ffunct5975H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct5975H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct5975H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6025H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6025H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6025H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6075H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6075H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6075H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6125H3->Get("ZimpCutat95H"))->GetMean();    
      ZimpCutat95DataL=((TH1F*)ffunct6125H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6125H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6175H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6175H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6175H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6225H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6225H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6225H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;  
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6275H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6275H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6275H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6325H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6325H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data;
      ZimpCutat95DataR=((TH1F*)ffunct6325H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);

      ZimpCutat95Data=((TH1F*)ffunct6375H3->Get("ZimpCutat95H"))->GetMean();
      ZimpCutat95DataL=((TH1F*)ffunct6375H3->Get("MuofDG"))->GetMean()-ZimpCutat95Data; 
      ZimpCutat95DataR=((TH1F*)ffunct6375H3->Get("MuofDG"))->GetMean()+ZimpCutat95Data;
      ZimpCutat95DataLVH3.push_back(ZimpCutat95DataL);ZimpCutat95DataRVH3.push_back(ZimpCutat95DataR);ZimpCutat95DataVH3.push_back(ZimpCutat95Data);
      





















 etamax=7.5;
 etamin=4.5;
 totnumberbinfordndeta=480;
 maxetafordndetahisto=12.0; //range assumed simmetric around eta=0
 etabinsize=(maxetafordndetahisto*2)/totnumberbinfordndeta;


 
 char sZname2[1024];
 char sZnamehist[1024];
 

 T2GeometryUtil conver;  

 



 
 std::ostringstream ssE;
 ssE<<PtCutinPrimaryEfficiency;//EnergyCutinPrimaryEfficiency;

 string globZc="Zcut=";
 std::ostringstream ss;
 for(unsigned int ii=0;ii<ZEffiCutImpact.size();ii++){
   ss<<ZEffiCutImpact.at(ii);
   ss<<"|";   
 }
 globZc="|"+globZc+ss.str()+"_ECut="+ssE.str();   

 for(unsigned int m=0;m<50;m++)//50=Number of Geant "Primary Particle"
    {

      sprintf(sZname2, "ArmPlus_TrkEtaEfficiency_FixedGeantMult %d", m); 
      sprintf(sZnamehist, "ArmPlus_TrkEtaEfficiency_FixedGeantMult-%d", m); 
     
      ArmPlus_TrkRecoCutMult_FixedGeantMult[m]=std::auto_ptr<TH1D>(new TH1D(sZname2,sZnamehist,100,-0.5,99.5));
      ArmPlus_TrkRecoCutMult_FixedGeantMult[m]->SetDirectory(0);

      sprintf(sZname2, "ArmPlus_TrkEtaEfficiency_FixedGeantMult2 %d", m); 
      sprintf(sZnamehist, "ArmPlus_TrkEtaEfficiency_FixedGeantMult2-%d", m); 
     
      ArmPlus_TrkRecoCutMult_FixedGeantMult2[m]=std::auto_ptr<TH1D>(new TH1D(sZname2,sZnamehist,100,-0.5,99.5));
      ArmPlus_TrkRecoCutMult_FixedGeantMult2[m]->SetDirectory(0);
      
      sprintf(sZname2, "ArmPlus_TrkRecoCutMult_FixedGeantMult %d", m); 
      sprintf(sZnamehist, "ArmPlus_TrkRecoCutMult_FixedGeantMult-%d", m);
      ArmPlus_TrkEtaEfficiency_FixedGeantMult[m]= std::auto_ptr<TProfile>(new TProfile(sZname2,sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
      ArmPlus_TrkEtaEfficiency_FixedGeantMult[m]->SetDirectory(0);

      //PrimaryRecoDnDeta2NoCut= std::auto_ptr<TProfile> (new TProfile("PrimaryRecoDnDeta2NoCut"," DN/D#eta of Reconstructed Primary Trk (No cuts) vs #eta",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

    }


for(unsigned int m=0;m<10; m++)
   {
     sprintf(sZname2, "H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult %d", m); 
     sprintf(sZnamehist, "H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult-AvgPadClu= %d", m); 
     string histostr(sZnamehist);
     histostr=histostr+globZc;

     H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,histostr.c_str(),totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

    
     sprintf(sZname2, "H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult %d", m); 
     sprintf(sZnamehist, "H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult-AvgPadClu= %d", m); 
     H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,histostr.c_str(),totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

     sprintf(sZname2, "H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult %d", m); 
     sprintf(sZnamehist, "H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult-AvgPadClu= %d", m); 
     H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,histostr.c_str(),totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

     sprintf(sZname2, "H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_PhiCut %d", m); 
     sprintf(sZnamehist, "H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_PhiCut-AvgPadClu= %d", m); 
     H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_PhiCut[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,histostr.c_str(),totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));


     
     sprintf(sZname2, "H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult %d", m); 
     sprintf(sZnamehist, "H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult-AvgPadClu= %d", m); 
     H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,histostr.c_str(),totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
     
     sprintf(sZname2, "H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_ZcutOnGeant %d", m); 
     sprintf(sZnamehist, "H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_ZcutOnGeant-AvgPadClu= %d", m); 
     H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_ZcutOnGeant[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));


     sprintf(sZname2, "H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_56Division %d", m); 
     sprintf(sZnamehist, "H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_56Division %d", m);
     H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_56Division[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));


     sprintf(sZname2, "H0_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult %d", m); 
     sprintf(sZnamehist, "H0_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult %d", m);
     H0_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,sZnamehist,15,0.,360.));

     sprintf(sZname2, "H1_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult %d", m); 
     sprintf(sZnamehist, "H1_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult %d", m);
     H1_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult[m]= std::auto_ptr<TProfile> (new TProfile(sZname2,sZnamehist,15,0.,360.));

     
   }
 
     H0_TrkPhiEfficiencyCumulative_UnfoldCutsCumul= std::auto_ptr<TProfile> (new TProfile("H0_TrkPhiEfficiencyCumulative_UnfoldCutsCumul","H0_TrkPhiEfficiencyCumulative_UnfoldCutsCumul",60,0.,360.));
     
     H1_TrkPhiEfficiencyCumulative_UnfoldCutsCumul= std::auto_ptr<TProfile> (new TProfile("H1_TrkPhiEfficiencyCumulative_UnfoldCutsCumul","H1_TrkPhiEfficiencyCumulative_UnfoldCutsCumul",60,0.,360.));

     NumEventGeneratorTrigger= std::auto_ptr<TH1D> (new TH1D("NumEventGeneratorTrigger","NumEventGeneratorTrigger",3,-0.5,2.5)); 
     NumEventTrackTrigger= std::auto_ptr<TH1D> (new TH1D("NumEventTrackTrigger","NumEventTrackTrigger",3,-0.5,2.5));


for(unsigned int m=0;m<50; m++)
   {
     sprintf(sZname2, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadColvsPadRow %d", m); 
     sprintf(sZnamehist, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadColvsPadRow-%d", m); 
     
     SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadColvsPadRow[m]=std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,66,-0.5,65.5,25,-0.5,24.5));
     SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadColvsPadRow[m]->SetDirectory(0);

     sprintf(sZname2, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOXY %d", m); 
     sprintf(sZnamehist, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOXY-%d", m); 
     SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOXY[m]=std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,280,-140.,140.,280,-140.,140.));
     SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOXY[m]->SetDirectory(0);
   
     sprintf(sZname2, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEORPhi %d", m); 
     sprintf(sZnamehist, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEORPhi-%d", m); 
     SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEORPhi[m]=std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,120,-0.5,359.5,70,40.,150.));
     SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEORPhi[m]->SetDirectory(0);
   
     sprintf(sZname2, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOXZ %d", m); 
     sprintf(sZnamehist, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOXZ-%d", m); 
     SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOXZ[m]=std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,280,-140.,140.));

     sprintf(sZname2, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOYZ %d", m); 
     sprintf(sZnamehist, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOYZ-%d", m); 
     SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOYZ[m]=std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,280,-140.,140.));

     sprintf(sZname2, "SelectedEventGeantPadTrkWithMoreThan3TracingGEOYZ_NoelectronPrimary %d", m); 
     sprintf(sZnamehist, "SelectedEventGeantPadTrkWithMoreThan3TracingGEOYZ_NoelectronPrimary-%d", m); 
     SelectedEventGeantPadTrkWithMoreThan3TracingGEOYZ_NoelectronPrimary[m]=std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,280,-140.,140.));

     sprintf(sZname2, "SelectedEventGeantPadTrkWithMoreThan3TracingGEOXZ_NoelectronPrimary %d", m); 
     sprintf(sZnamehist, "SelectedEventGeantPadTrkWithMoreThan3TracingGEOXZ_NoelectronPrimary-%d", m); 
     SelectedEventGeantPadTrkWithMoreThan3TracingGEOXZ_NoelectronPrimary[m]=std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,280,-140.,140.));

     sprintf(sZname2, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadRowvsZ %d", m); 
     sprintf(sZnamehist, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadRowvsZ-%d", m); 
     SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadRowvsZ[m]=std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,25,-0.5,24.5));
     sprintf(sZname2, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadColvsZ %d", m); 
     sprintf(sZnamehist, "SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadColvsZ-%d", m); 
     SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadColvsZ[m]=std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,66,-0.5,65.5));  

     sprintf(sZname2, "SelectedEventGeantPadTrkWithMoreThan3TracingPadRowvsZ_NoelectronPrimary %d", m); 
     sprintf(sZnamehist, "SelectedEventGeantPadTrkWithMoreThan3TracingPadRowvsZ_NoelectronPrimary-%d", m); 
     SelectedEventGeantPadTrkWithMoreThan3TracingPadRowvsZ_NoelectronPrimary[m]=std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,25,-0.5,24.5));  
       
     sprintf(sZname2, "SelectedEventGeantPadTrkWithMoreThan3TracingPadRowvsZ_NoelectronPrimary %d", m); 
     sprintf(sZnamehist, "SelectedEventGeantPadTrkWithMoreThan3TracingPadRowvsZ_NoelectronPrimary-%d", m); 
     SelectedEventGeantPadTrkWithMoreThan3TracingPadColvsZ_NoelectronPrimary[m]=std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,66,-0.5,65.5));  


     sprintf(sZname2, "SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY_NoelectronPrimary %d", m); 
     sprintf(sZnamehist, "SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY_NoelectronPrimary-%d", m); 
     SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY_NoelectronPrimary[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,140,-140.,140.,140,-140.,140.));



     
     sprintf(sZname2, "SelectedEventDIGIPad_GTrkWithMoreThan3_GEOYZ %d", m);
     sprintf(sZnamehist, "SelectedEventDIGIPad_GTrkWithMoreThan3_GEOYZ-%d", m); 
     SelectedEventDIGIPad_GTrkWithMoreThan3_GEOYZ[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,280,-140.,140.));

     
     sprintf(sZname2, "SelectedEventDIGIPad_GTrkWithMoreThan3_GEOXZ %d", m);
     sprintf(sZnamehist, "SelectedEventDIGIPad_GTrkWithMoreThan3_GEOXZ-%d", m); 
     SelectedEventDIGIPad_GTrkWithMoreThan3_GEOXZ[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,280,-140.,140.));

     sprintf(sZname2, "SelectedEventDIGIPad_GTrkWithMoreThan3_PadColvsPadRow %d", m);
     sprintf(sZnamehist, "SelectedEventDIGIPad_GTrkWithMoreThan3_PadColvsPadRow-%d", m); 
     SelectedEventDIGIPad_GTrkWithMoreThan3_PadColvsPadRow[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,66,-0.5,65.5,25,-0.5,24.5));

     sprintf(sZname2, "SelectedEventDIGIPad_GTrkWithMoreThan3_PadColvsZ %d", m);
     sprintf(sZnamehist, "SelectedEventDIGIPad_GTrkWithMoreThan3_PadColvsZ-%d", m); 
     SelectedEventDIGIPad_GTrkWithMoreThan3_PadColvsZ[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,66,-0.5,65.5));

     sprintf(sZname2, "SelectedEventDIGIPad_GTrkWithMoreThan3_PadRowvsZ %d", m);
     sprintf(sZnamehist, "SelectedEventDIGIPad_GTrkWithMoreThan3_PadRowvsZ-%d", m); 
     SelectedEventDIGIPad_GTrkWithMoreThan3_PadRowvsZ[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,25,-0.5,24.5));

     sprintf(sZname2, "SelectedEventDIGIPad_GTrkWithMoreThan3_GEOXY %d", m);
     sprintf(sZnamehist, "SelectedEventDIGIPad_GTrkWithMoreThan3_GEOXY-%d", m); 
     SelectedEventDIGIPad_GTrkWithMoreThan3_GEOXY[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,280,-140.,140.,280,-140.,140.));






     sprintf(sZname2, "SelectedEventRoadPadFinderClu_PadColvsZ %d", m);
     sprintf(sZnamehist, "SelectedEventRoadPadFinderClu_PadColvsZ-%d", m); 
     SelectedEventRoadPadFinderClu_PadColvsZ[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,66,-0.5,65.5));     


     sprintf(sZname2, "SelectedEventRoadPadFinderClu_PadRowvsZ %d", m);
     sprintf(sZnamehist, "SelectedEventRoadPadFinderClu_PadRowvsZ-%d", m); 
     SelectedEventRoadPadFinderClu_PadRowvsZ[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,25,-0.5,24.5));  
 
     sprintf(sZname2, "SelectedEventRoadPadFinderClu_GEOXZ %d", m);
     sprintf(sZnamehist, "SelectedEventRoadPadFinderClu_GEOXZ-%d", m); 
     SelectedEventRoadPadFinderClu_GEOXZ[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,280,-140.,140.));
      
     sprintf(sZname2, "SelectedEventRoadPadFinderClu_GEOYZ %d", m);
     sprintf(sZnamehist, "SelectedEventRoadPadFinderClu_GEOYZ-%d", m); 
     SelectedEventRoadPadFinderClu_GEOYZ[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,280,-140.,140.));
      
      
     sprintf(sZname2, "SelectedEventRoadPadFinder_FullR_Clu_GEOYZ %d", m);
     sprintf(sZnamehist, "SelectedEventRoadPadFinder_FullR_Clu_GEOYZ-%d", m); 
     SelectedEventRoadPadFinder_FullR_Clu_GEOYZ[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,280,-140.,140.));
     

     
     sprintf(sZname2, "SelectedEventRoadPadFinder_FullR_Clu_PadColvsZ %d", m);
     sprintf(sZnamehist, "SelectedEventRoadPadFinder_FullR_Clu_PadColvsZ-%d", m); 
     SelectedEventRoadPadFinder_FullR_Clu_PadColvsZ[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,66,-0.5,65.5));     

     sprintf(sZname2, "SelectedEventRoadPadFinder_FullR_Clu_PadRowvsZ %d", m);
     sprintf(sZnamehist, "SelectedEventRoadPadFinder_FullR_Clu_PadRowvsZ-%d", m); 
     SelectedEventRoadPadFinder_FullR_Clu_PadRowvsZ[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,25,-0.5,24.5));  
 
     sprintf(sZname2, "SelectedEventRoadPadFinder_FullR_Clu_GEOXZ %d", m);
     sprintf(sZnamehist, "SelectedEventRoadPadFinder_FullR_Clu_GEOXZ-%d", m); 
     SelectedEventRoadPadFinder_FullR_Clu_GEOXZ[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,280,-140.,140.));
      
     sprintf(sZname2, "SelectedEventRoadPadFinder_FullR_Clu_GEOYZ %d", m);
     sprintf(sZnamehist, "SelectedEventRoadPadFinder_FullR_Clu_GEOYZ-%d", m); 
     SelectedEventRoadPadFinder_FullR_Clu_GEOYZ[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,80,13800,14300,280,-140.,140.));
      

    
    
   }




 
 DNDetaT2ALL= std::auto_ptr<TH1D> (new TH1D("DNDetaT2ALL","DNDetaT2ALL",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 
 DNDetaT2Show= std::auto_ptr<TH1D> (new TH1D("DNDetaT2Show","DNDetaT2Show",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 
 DNDetaT2Clean= std::auto_ptr<TH1D> (new TH1D("DNDetaT2Clean","DNDetaT2Clean",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 
 
 Histonumtoeventgen= std::auto_ptr<TH1D> (new TH1D("Histonumtoeventgen","Histonumtoeventgen",2,-0.5,1.5));
 Histonumshowevent= std::auto_ptr<TH1D> (new TH1D("Histonumshowevent","Histonumshowevent",2,-0.5,1.5));  
 Histonumcleanevent= std::auto_ptr<TH1D> (new TH1D("Histonumcleanevent","Histonumcleanevent",2,-0.5,1.5));
 
 


 NumChPartInT2PerEvt= std::auto_ptr<TH1D>(new TH1D("NumChPartInT2PerEvt","NumChPartInT2PerEvt",300,-0.5, 299.5)); 
 NumTrkInT2PerEvt= std::auto_ptr<TH1D>(new TH1D("NumTrkInT2PerEvt","NumTrkInT2PerEvt",300,-0.5, 299.5)); 
 NumTrkInT2PerEvt7mcut= std::auto_ptr<TH1D>(new TH1D("NumTrkInT2PerEvt7mcut","NumTrkInT2PerEvt7mcut",300,-0.5, 299.5));

 PSimHitT2_EnergyLoss= std::auto_ptr<TH1D>(new TH1D("PSimHitT2_EnergyLoss","PSimHitT2_EnergyLoss",1000,0., 10000));//Interesting: 1ev-10KeV
 PSimHitT2_EnergyLoss->SetXTitle("E-loss (eV)");

 GeantThetaXofPrimary= std::auto_ptr<TH1D>(new TH1D("GeantThetaXofPrimary","GeantThetaXofPrimary",100,-0.05, 0.05));
 GeantThetaYofPrimary= std::auto_ptr<TH1D>(new TH1D("GeantThetaYofPrimary","GeantThetaYofPrimary",100,-0.05, 0.05));

 RecoTrkEta_MatchingPrimary= std::auto_ptr<TH1D>(new TH1D("RecoTrkEta_MatchingPrimary","RecoTrkEta_MatchingPrimary",100,-9.05, 9.05));
 
 PID_GeantTrk_atleast4HitInT2= std::auto_ptr<TH1D>(new TH1D("PID_GeantTrk_atleast4HitInT2","Geant Trk Particle ID, (at least 4 Hits)",11001,-5500.5, 5500.5));

  PID_GeantTrk_atleast4HitInT2_StableAtIP= std::auto_ptr<TH1D>(new TH1D("PID_GeantTrk_atleast4HitInT2_StableAtIP","Geant Trk Particle ID, (at least 4 Hits) Stable at IP",11001,-5500.5, 5500.5));

  PID_GeantTrk_atleast4HitInT2_ConsideredAsPrimary= std::auto_ptr<TH1D>(new TH1D("PID_GeantTrk_atleast4HitInT2_ConsideredAsPrimary","Geant Trk Particle ID, (at least 4 Hits) Considered As Primary",11001,-5500.5, 5500.5));

 Pythia_PDG_HavingStableDesc_Theta003= std::auto_ptr<TH1D>(new TH1D("Pythia_PDG_HavingStableDesc_Theta003","Pythia Unstable Particle ID, |#theta|<0.03",11001,-5500.5, 5500.5));  

 Pythia_PDG_Status2_Theta003= std::auto_ptr<TH1D>(new TH1D("Pythia_PDG_Status2_Theta003","Pythia Particle ID Status 2 |#theta|<0.03",11001,-5500.5, 5500.5)); 
  
 Pythia_PDG_HavingStableDesc_Theta001= std::auto_ptr<TH1D>(new TH1D("Pythia_PDG_HavingStableDesc_Theta001","Pythia Unstable Particle ID, |#theta|<0.01",11001,-5500.5, 5500.5));

 Pythia_UnstablePDG= std::auto_ptr<TH1D>(new TH1D("Pythia_UnstablePDG","Pythia Unstable Particle ID",11001,-5500.5, 5500.5)); 
 Pythia_PDG_HavingStableDesc= std::auto_ptr<TH1D>(new TH1D("Pythia_PDG_HavingStableDesc","Particle ID of mother particle having stable daughters",11001,-5500.5, 5500.5)); 
 
 PID_ofSecRecoTrack=std::auto_ptr<TH1D>(new TH1D("PID_ofSecRecoTrack","Geant PID of secondary tracks",11001,-5500.5, 5500.5)); 
 
 
 H0_sameTrk_StripNumb= std::auto_ptr<TH1D>(new TH1D("H0_sameTrk_StripNumb","Strip number active in the same trk",512,-0.5, 511.5));
 
 
 H0_sameTrk_PadRow= std::auto_ptr<TH1D>(new TH1D("H0_sameTrk_PadRow","Pad Row active",26,-0.5, 25.5));
  
  H0_sameTrk_PadCol= std::auto_ptr<TH1D>(new TH1D("H0_sameTrk_PadCol","Pad Col active",66,-0.5, 65.5));

  H0_sameTrk_PadColMaxDist= std::auto_ptr<TH1D>(new TH1D("H0_sameTrk_PadColMaxDist","Pad Col max distance in the same trk",66,-0.5, 65.5));

  H0_sameTrk_PadRowMaxDist= std::auto_ptr<TH1D>(new TH1D("H0_sameTrk_PadRowMaxDist","Pad Row max distance in the same trk",26,-0.5, 25.5));

  
  H0_sameTrk_StripRowMaxDist= std::auto_ptr<TH1D>(new TH1D("H0_sameTrk_StripRowMaxDist","Strip Row max distance in the same trk",26,-0.5, 25.5));

  H0_sameTrk_StripColMaxDist= std::auto_ptr<TH1D>(new TH1D("H0_sameTrk_StripColMaxDist","Strip Col max distance in the same trk",26,-0.5, 25.5));
  
  DigiPadOccupancy= std::auto_ptr<TProfile>(new TProfile("DigiPadOccupancy","Plane Pad occupancy",40,-0.5,39.5,""));
  DigiPadOccupancy->GetYaxis()->SetTitle("Occupancy ");
  DigiPadOccupancy->GetXaxis()->SetTitle("Detector plane");
  DigiPadOccupancy->SetDirectory(0);

  CumulativeGeantTrkHadronR=std::auto_ptr<TH1D>(new TH1D("CumulativeGeantTrkHadronR","CumulativeGeantTrkHadronR",150,0,150)); 
  CumulativeGeantTrkElectronR=std::auto_ptr<TH1D>(new TH1D("CumulativeGeantTrkElectronR","CumulativeGeantTrkElectronR",150,0,150));
  CumulativeGeantTrkElectronXY=std::auto_ptr<TH2D>(new TH2D("CumulativeGeantTrkElectronXY","CumulativeGeantTrkElectronXY",150,-150,150,150,-150,150));

  PID_DircetMotherForAFakeRecoPrimaryTrk=std::auto_ptr<TH1D>(new TH1D("PID_DircetMotherForAFakeRecoPrimaryTrk","Geant PID of secondary tracks Mother Reco as primary",11001,-5500.5, 5500.5)); 

  PID_OldestMotherForAFakeRecoPrimaryTrk=std::auto_ptr<TH1D>(new TH1D("PID_OldestMotherForAFakeRecoPrimaryTrk","Geant PID of secondary tracks older Mother Reco as primary",11001,-5500.5, 5500.5)); 
  
  PID_ofSecRecoTrackAsaFakePrimary=std::auto_ptr<TH1D>(new TH1D("PID_ofSecRecoTrackAsaFakePrimary","Geant PID of secondary tracks  Reco as primary",11001,-5500.5, 5500.5)); 


  SecondaryTrksEtaFromIonPump= std::auto_ptr<TH1D>(new TH1D("SecondaryTrksEtaFromIonPump","SecondaryTrksEtaFromIonPump",240, -8, 8)); 
    
  SecondaryTrksZvsPId= std::auto_ptr<TH2D>(new TH2D("SecondaryTrksZvsPId","SecondaryTrksZvsPId",700, -3500, 3500,2500,0,2500)); 
    
  SecondaryTrksZvsEta= std::auto_ptr<TH2D>(new TH2D("SecondaryTrksZvsEta","SecondaryTrksZvsEta",700, 0, 3500,240, -8, 8));

  SimuVtxPositionForAFakeRecoPrimaryTrk= std::auto_ptr<TH2D>(new TH2D("SimuVtxPositionForAFakeRecoPrimaryTrk","Geant Vertex Position for fake primary reconstructed tracks",400, -2000, 2000,500,-500,500)) ;

  PdgSecondaryInT2EtaAcceptance= std::auto_ptr<TH1D>(new TH1D("PdgSecondaryInT2EtaAcceptance","PdgSecondaryInT2EtaAcceptance",7000, -3500, 3500));
  ParticlePdgvsEinT2 = std::auto_ptr<TH2F>(new TH2F("ParticlePdgvsEinT2","ParticlePdgvsEinT2",700, -3500, 3500,700,0,7000));
  ParticlePdgvsEinT2->SetDirectory(0); 
  ParticlePdgvsENotinT2 = std::auto_ptr<TH2F>(new TH2F("ParticlePdgvsENotinT2","ParticlePdgvsENotinT2",700, -3500, 3500,700,0,7000));
  ParticlePdgvsENotinT2->SetDirectory(0);
  PrimaryTrksEtavsZatRmin = std::auto_ptr<TH2F>(new TH2F("PrimaryTrksEtavsZatRmin","PrimaryTrksEtavsZatRmin",120,-8,8,300,-1500,1500));
  PrimaryTrksEtavsZatRmin->SetDirectory(0);
  PrimaryTrksEta = std::auto_ptr<TH1D>(new TH1D("PrimaryTrksEta","PrimaryTrksEta",120, -8, 8));
  PrimaryTrksEta->SetDirectory(0);
  AllTrksEta= std::auto_ptr<TH1D>(new TH1D("AllTrksEta","AllTrkEta",120, -8, 8));
  AllTrksEta->SetDirectory(0);
  SecondaryTrksEta= std::auto_ptr<TH1D>(new TH1D("SecondaryTrksEta","SecondaryTrksEta",240, -8, 8));
  SecondaryTrksEta->SetDirectory(0);

  
  
  Reco_Trk_DeltaThetaX_vs_DeltaThetaY_FromPrimary = std::auto_ptr<TH2D>(new TH2D("Reco_Trk_DeltaThetaX_vs_DeltaThetaY_FromPrimary","Reco_Trk_DeltaThetaX_vs_DeltaThetaY_FromPrimary",200,-0.1,0.1,200,-0.1,0.1));

  
  Reco_Trk_DeltaThetaX_vs_DeltaThetaY_FromSecondary = std::auto_ptr<TH2D>(new TH2D("Reco_Trk_DeltaThetaX_vs_DeltaThetaY_FromSecondary","Reco_Trk_DeltaThetaX_vs_DeltaThetaY_FromSecondary",200,-0.1,0.1,200,-0.1,0.1));

 

   H0_SecondaryRecoTrksZ0_Vs_Z0OrtogImpact_BisMatch=  std::auto_ptr<TH2D>(new TH2D("H0_SecondaryRecoTrksZ0_Vs_Z0OrtogImpact_BisMatch", "H0_SecondaryRecoTrksZ0_Vs_Z0OrtogImpact_BisMatch", 300, -2000, 2000,300, -2000, 2000)) ;
 
   H0_PrimaryRecoTrksZ0_Vs_Z0OrtogImpact_BisMatch=  std::auto_ptr<TH2D>(new TH2D("H0_PrimaryRecoTrksZ0_Vs_Z0OrtogImpact_BisMatch", "H0_PrimaryRecoTrksZ0_Vs_Z0OrtogImpact_BisMatch", 300, -2000, 2000,300, -2000, 2000)) ;

   H0_SecondaryRecoTrksZ0Ortog_Vs_Y0Ortog_BisMatch = std::auto_ptr<TH2D>(new TH2D("H0_SecondaryRecoTrksZ0Ortog_Vs_Y0Ortog_BisMatch", "H0_SecondaryRecoTrksZ0Ortog_Vs_Y0Ortog_BisMatch", 300, -2000, 2000,300, -50, 50)); 
   H0_PrimaryRecoTrksZ0Ortog_Vs_Y0Ortog_BisMatch = std::auto_ptr<TH2D>(new TH2D("H0_PrimaryRecoTrksZ0Ortog_Vs_Y0Ortog_BisMatch", "H0_PrimaryRecoTrksZ0Ortog_Vs_Y0Ortog_BisMatch", 300, -2000, 2000,300, -50, 50));

   H0_SecondaryRecoTrksZ0Ortog_Vs_X0Ortog_BisMatch  = std::auto_ptr<TH2D>(new TH2D("H0_SecondaryRecoTrksZ0Ortog_Vs_X0Ortog_BisMatch", "H0_SecondaryRecoTrksZ0Ortog_Vs_X0Ortog_BisMatch", 300, -2000, 2000,300, -50, 50));   
   H0_PrimaryRecoTrksZ0Ortog_Vs_X0Ortog_BisMatch = std::auto_ptr<TH2D>(new TH2D("H0_PrimaryRecoTrksZ0Ortog_Vs_X0Ortog_BisMatch", "H0_PrimaryRecoTrksZ0Ortog_Vs_X0Ortog_BisMatch", 300, -2000, 2000,300, -50, 50));

   H0_SecondaryRecoTrksX0Ortog_Vs_Y0Ortog_BisMatch  = std::auto_ptr<TH2D>(new TH2D("H0_SecondaryRecoTrksX0Ortog_Vs_Y0Ortog_BisMatch", "H0_SecondaryRecoTrksX0Ortog_Vs_Y0Ortog_BisMatch", 300, -2000, 2000,300, -50, 50));
   H0_PrimaryRecoTrksX0Ortog_Vs_Y0Ortog_BisMatch = std::auto_ptr<TH2D>(new TH2D("H0_PrimaryRecoTrksX0Ortog_Vs_Y0Ortog_BisMatch", "H0_PrimaryRecoTrksX0Ortog_Vs_Y0Ortog_BisMatch", 300, -2000, 2000,300, -50, 50));

   H0_SecondaryRecoTrksZ0Ortog_Vs_R0Ortog_BisMatch  = std::auto_ptr<TH2D>(new TH2D("H0_SecondaryRecoTrksZ0Ortog_Vs_R0Ortog_BisMatch", "H0_SecondaryRecoTrksZ0Ortog_Vs_R0Ortog_BisMatch", 300, -2000, 2000,300, 0, 50));
   H0_PrimaryRecoTrksZ0Ortog_Vs_R0Ortog_BisMatch = std::auto_ptr<TH2D>(new TH2D("H0_PrimaryRecoTrksZ0Ortog_Vs_R0Ortog_BisMatch", "H0_PrimaryRecoTrksZ0Ortog_Vs_R0Ortog_BisMatch", 300, -2000, 2000,300, 0., 50));
  
  

  Reco_TrkZ0_vs_R0_AllH0_FromPrimary = std::auto_ptr<TH2D>(new TH2D("Reco_TrkZ0_vs_R0_AllH0_FromPrimary","Track Z at Rmin (x-y) vs R@IP",100,-25000,25000,100,0,100));
  Reco_TrkZ0_vs_R0_AllH0_FromPrimary->SetXTitle("Z (mm)");
  Reco_TrkZ0_vs_R0_AllH0_FromPrimary->SetYTitle("R @ IP (mm)");

  Reco_TrkZ0_vs_R0_AllH0_FromSecondary = std::auto_ptr<TH2D>(new TH2D("Reco_TrkZ0_vs_R0_AllH0_FromSecondary","Track Z at Rmin (x-y) vs R@IP",100,-25000,25000,100,0,100));
  Reco_TrkZ0_vs_R0_AllH0_FromSecondary->SetXTitle("Z (mm)");
  Reco_TrkZ0_vs_R0_AllH0_FromSecondary->SetYTitle("R @ IP (mm)");

  Reco_TrkZ0_vs_Eta_AllH0_FromPrimary = std::auto_ptr<TH2D>(new TH2D("Reco_TrkZ0_vs_Eta_AllH0_FromPrimary","Track Z at Rmin (x-y) vs #eta",100,-25000,25000,100,-10,10));
  Reco_TrkZ0_vs_Eta_AllH0_FromPrimary->SetXTitle("Z (mm)");


  Reco_TrkZ0_vs_Eta_AllH0_FromSecondary = std::auto_ptr<TH2D>(new TH2D("Reco_TrkZ0_vs_Eta_AllH0_FromSecondary","Track Z at Rmin (x-y) vs #eta",100,-25000,25000,100,-10,10));
  Reco_TrkZ0_vs_Eta_AllH0_FromSecondary->SetXTitle("Z (mm)");
  

  Reco_TrkZ0_vs_Rmin_AllH0_FromSecondary = std::auto_ptr<TH2D>(new TH2D("Reco_TrkZ0_vs_Rmin_AllH0_FromSecondary","Track Z at Rmin (x-y) vs R_min",100,-25000,25000,100,0,100));
  Reco_TrkZ0_vs_Rmin_AllH0_FromSecondary->SetXTitle("Z (mm)");
  Reco_TrkZ0_vs_Rmin_AllH0_FromSecondary->SetYTitle("R_min (cm)");
  
  
  Reco_TrkZ0_vs_Rmin_AllH0_FromPrimary = std::auto_ptr<TH2D>(new TH2D("Reco_TrkZ0_vs_Rmin_AllH0_FromPrimary","Track Z at Rmin (x-y) vs R_min",100,-25000,25000,100,0,100));
  Reco_TrkZ0_vs_Rmin_AllH0_FromPrimary->SetXTitle("Z (mm)");
  Reco_TrkZ0_vs_Rmin_AllH0_FromPrimary->SetYTitle("R_min (cm)");
    

  SelectedEventGeantEntryPad= std::auto_ptr<TH2F>(new TH2F("SelectedEventGeantEntryPad","SelectedEventGeantEntryPad",66,-0.5,65.5,25,-0.5,24.5)); 
  SelectedEventGeantExitPad= std::auto_ptr<TH2F>(new TH2F("SelectedEventGeantExitPad","SelectedEventGeantExitPad",66,-0.5,65.5,25,-0.5,24.5));

  SelectedEventGeantPad= std::auto_ptr<TH2F>(new TH2F("SelectedEventGeantPad","SelectedEventGeantPad",66,-0.5,65.5,25,-0.5,24.5));
  SelectedEventGeantPad->SetXTitle("Pad Col (#Phi direction)");
  SelectedEventGeantPad->SetYTitle("Pad Row (Radial direction)");

  
  SelectedEventGeantPadTrkWithMoreThan3_PythiaChinIP= std::auto_ptr<TH2F>(new TH2F("SelectedEventGeantPadTrkWithMoreThan3_PythiaChinIP","SelectedEventGeantPadTrkWithMoreThan3_PythiaChinIP",66,-0.5,65.5,25,-0.5,24.5));
  SelectedEventGeantPadTrkWithMoreThan3_PythiaChinIP->SetXTitle("Pad Col (#Phi direction)");
  SelectedEventGeantPad->SetYTitle("Pad Row (Radial direction)");


  SelectedEventGeantPadTrkWithMoreThan3Hit= std::auto_ptr<TH2F>(new TH2F("SelectedEventGeantPadTrkWithMoreThan3Hit","SelectedEventGeantPadTrkWithMoreThan3Hit",66,-0.5,65.5,25,-0.5,24.5));
  SelectedEventGeantPadTrkWithMoreThan3Hit->SetXTitle("Pad Col (#Phi direction)");
  SelectedEventGeantPadTrkWithMoreThan3Hit->SetYTitle("Pad Row (Radial direction)");

  SelectedEventGeantPadTrkWithMoreThan3HitEntry_Weight= std::auto_ptr<TH2F>(new TH2F("SelectedEventGeantPadTrkWithMoreThan3HitEntry_Weight","SelectedEventGeantPadTrkWithMoreThan3HitEntry_Weight",66,-0.5,65.5,25,-0.5,24.5));

  SelectedEventGeantPadTrkWithMoreThan3HitExit_Weight= std::auto_ptr<TH2F>(new TH2F("SelectedEventGeantPadTrkWithMoreThan3HitExit_Weight","SelectedEventGeantPadTrkWithMoreThan3HitExit_Weight",66,-0.5,65.5,25,-0.5,24.5));


SelectedEventGeantPadTrkWithMoreThan3TracingMultiCount= std::auto_ptr<TH2F>(new TH2F("SelectedEventGeantPadTrkWithMoreThan3TracingMultiCount","SelectedEventGeantPadTrkWithMoreThan3TracingMultiCount",66,-0.5,65.5,25,-0.5,24.5)); 

 SelectedEventGeantPadTrkWithMoreThan3Tracing= std::auto_ptr<TH2F>(new TH2F("SelectedEventGeantPadTrkWithMoreThan3Tracing","SelectedEventGeantPadTrkWithMoreThan3Tracing",66,-0.5,65.5,25,-0.5,24.5)); 
 
 SelectedEventGeantPadTrkWithMoreThan3Tracing->SetXTitle("Pad Col (#Phi direction)");
 SelectedEventGeantPadTrkWithMoreThan3Tracing->SetYTitle("Pad Row (Radial direction)");
 
 SelectedEventGeantPadTrkWithMoreThan3TracingGEO= std::auto_ptr<TH2F>(new TH2F("SelectedEventGeantPadTrkWithMoreThan3TracingGEO","SelectedEventGeantPadTrkWithMoreThan3TracingGEO",120,-0.5,359.5,70,40.,150.));
 SelectedEventGeantPadTrkWithMoreThan3TracingGEO->SetXTitle("Pad Col #Phi (deg)");
 SelectedEventGeantPadTrkWithMoreThan3TracingGEO->SetYTitle("Pad Col R (mm)");
  



 SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY= std::auto_ptr<TH2F>(new TH2F("SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY","SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY",140,-140.,140.,140,-140.,140.));
 SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY->SetXTitle("X mm");
 SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY->SetYTitle("Y mm");
 
 SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY_Noelectron= std::auto_ptr<TH2F>(new TH2F("SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY_Noelectron","SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY_Noelectron",140,-140.,140.,140,-140.,140.));
 SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY_Noelectron->SetXTitle("X mm");
 SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY_Noelectron->SetYTitle("Y mm");

 
 SelectedEventGeantPadTrkWithMoreThan3_PythiaChinIPGEOXY= std::auto_ptr<TH2F>(new TH2F("SelectedEventGeantPadTrkWithMoreThan3_PythiaChinIPGEOXY","PythiaChinIPGEOXY",140,-140.,140.,140,-140.,140.));
 SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY->SetXTitle("X mm");
 SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY->SetYTitle("Y mm");
  

 RecoTrksEtaALL_EtaCut4_R0Cut20= std::auto_ptr<TH1D>(new TH1D("RecoTrksEtaALL_EtaCut4_R0Cut20","RecoTrksEtaALL_EtaCut4_R0Cut20",240, -8, 8));
 SecondaryRecoTrksEta_2Arms_EtaCut4_R0Cut20= std::auto_ptr<TH1D>(new TH1D("SecondaryRecoTrksEta_2Arms_EtaCut4_R0Cut20","SecondaryRecoTrksEta_2Arms_EtaCut4_R0Cut20",240, -8, 8));
 PrimaryRecoTrksEta_2Arms_EtaCut4_R0Cut20= std::auto_ptr<TH1D>(new TH1D("PrimaryRecoTrksEta_2Arms_EtaCut4_R0Cut20","PrimaryRecoTrksEta_2Arms_EtaCut4_R0Cut20",240, -8, 8));
 
 NumHitInRecTrk= std::auto_ptr<TH1D>(new TH1D("NumHitInRecTrk","NumHitInRecTrk", 21,-0.5, 20.5));
 PrimaryVertexPositionYZ= std::auto_ptr<TH2D>(new TH2D("PrimaryVertexPositionYZ","PrimaryVertexPositionYZ",400, -2000, 2000,500,-500,500)) ;
 PrimaryVertexPositionYZ->SetDirectory(0);
 SecondaryVertexPositionYZ= std::auto_ptr<TH2D>(new TH2D("SecondaryVertexPositionYZ","SecondaryVertexPositionYZ",400, -2000, 2000,500,-500,500)) ;
  SecondaryVertexPositionYZ->SetDirectory(0);
  ALLVertexPositionYZ= std::auto_ptr<TH2D>(new TH2D("ALLVertexPositionYZ","ALLVertexPositionYZ",400, -2000, 2000,500,-50,50)) ;
  ALLVertexPositionYZ->SetDirectory(0);

  PrimaryVtxPositionContributingInT2= std::auto_ptr<TH2D>(new TH2D("PrimaryVtxPositionContributingInT2","PrimaryVtxPositionContributingInT2",400, -2000, 2000,500,-200,200)) ;
  SecondaryVtxPositionContributingInT2 = std::auto_ptr<TH2D>(new TH2D("SecondaryVtxPositionContributingInT2","SecondaryVtxPositionContributingInT2",400, -2000, 2000,100,-50,50)) ;
 
  SecondaryVtxPositionContributingInT2MultiCount= std::auto_ptr<TH2D>(new TH2D("SecondaryVtxPositionContributingInT2MultiCount","SecondaryVtxPositionContributingInT2MultiCount",400, -2000, 2000,100,-50,50)) ;


 SecondaryVtx_Z_ContributingInT2 = std::auto_ptr<TH1D>(new TH1D("SecondaryVtx_Z_ContributingInT2","SecondaryVtx_Z_ContributingInT2",300, -2000, 2000)) ; 
 SecondaryVtx_R_ContributingInT2 = std::auto_ptr<TH1D>(new TH1D("SecondaryVtx_R_ContributingInT2","SecondaryVtx_R_ContributingInT2",100,-50,50)) ;
  
 SecondaryVtx_ZvsR_ContributingInT2= std::auto_ptr<TH2D>(new TH2D("SecondaryVtx_ZvsR_ContributingInT2","SecondaryVtxPositionContributingInT2",4000, -2000, 2000,1000,-50,50)) ;

 RecoVtx_ZvsR= std::auto_ptr<TH2D>(new TH2D("RecoVtx_ZvsR","RecoVtx_ZvsR",400, -2000, 2000,500,-50,50)) ;
 CumulativeZimpH2= std::auto_ptr<TH1D>(new TH1D("CumulativeZimpH2","CumulativeZimpH2 ",500, -15000, 15000));
 CumulativeZimpH3= std::auto_ptr<TH1D>(new TH1D("CumulativeZimpH3","CumulativeZimpH3",500, -15000, 15000));

 CumulativeZimpH0= std::auto_ptr<TH1D>(new TH1D("CumulativeZimpH0","CumulativeZimpH0 ",500, -15000, 15000));
 CumulativeZimpH1= std::auto_ptr<TH1D>(new TH1D("CumulativeZimpH1","CumulativeZimpH1",500, -15000, 15000));
 AllRecoVtxYZ= std::auto_ptr<TH2D>(new TH2D("AllRecoVtxYZ","AllRecoVtxYZ",400, -2000, 2000,500,-50,50)) ;

 AllRecoVtxZ = std::auto_ptr<TH1D>(new TH1D("AllRecoVtxZ","AllRecoVtxZ",300, -2000, 2000)) ; 

 eta2geant_forEffHisto1= std::auto_ptr<TH1D>(new TH1D("eta2geant_forEffHisto1","eta2geant_forEffHisto1",240, -8, 8));
 eta2geant_forEffHisto2= std::auto_ptr<TH1D>(new TH1D("eta2geant_forEffHisto2","eta2geant_forEffHisto2",240, -8, 8));
 eta2geant_forEffHisto3= std::auto_ptr<TH1D>(new TH1D("eta2geant_forEffHisto3","eta2geant_forEffHisto3",240, -8, 8));
 
 RecoTrksEtaALL= std::auto_ptr<TH1D>(new TH1D("RecoTrksEtaALL","RecoTrksEtaALL ",240, -8, 8));
 
 SecondaryRecoTrksR02PointingToVtx = std::auto_ptr<TH1D>(new TH1D("SecondaryRecoTrksR02PointingToVtx","SecondaryRecoTrksR02PointingToVtx",400,0,40.));
 SecondaryRecoTrksZ0PointingToVtx= std::auto_ptr<TH1D>(new TH1D("SecondaryRecoTrksZ0PointingToVtx","SecondaryRecoTrksZ0PointingToVtx",300, -2000, 2000)); 
 SecondaryRecoTrksR0PointingToVtx= std::auto_ptr<TH1D>(new TH1D("SecondaryRecoTrksR0PointingToVtx","SecondaryRecoTrksR0PointingToVtx",100, 0, 20)); 
 SecondaryRecoTrksEtaPointingToVtx= std::auto_ptr<TH1D>(new TH1D("SecondaryRecoTrksEtaPointingToVtx","SecondaryRecoTrksEtaPointingToVtx ",240, -8, 8));
 SecondaryRecoTrksEta= std::auto_ptr<TH1D>(new TH1D("SecondaryRecoTrksEta","SecondaryRecoTrksEta ",240, -8, 8));

 PrimaryRecoTrksZ0PointingToVtx = std::auto_ptr<TH1D>(new TH1D("PrimaryRecoTrksZ0PointingToVtx","PrimaryRecoTrksZ0PointingToVtx",300, -2000, 2000));
 PrimaryRecoTrksrR02PointingToVtx  = std::auto_ptr<TH1D>(new TH1D("PrimaryRecoTrksrR02PointingToVtx","PrimaryRecoTrksrR02PointingToVtx",400,0,40.));
 PrimaryRecoTrksrR0PointingToVtx= std::auto_ptr<TH1D>(new TH1D("PrimaryRecoTrksrR0PointingToVtx","PrimaryRecoTrksrR0PointingToVtx",100,0,20.));

 SecondaryRecoTrksR02_2Arms  = std::auto_ptr<TH1D>(new TH1D("SecondaryRecoTrksR02_2Arms","SecondaryRecoTrksR02_2Arms",200,0,50.));
 SecondaryRecoTrksZ0_2Arms  = std::auto_ptr<TH1D>(new TH1D("SecondaryRecoTrksZ0_2Arms","SecondaryRecoTrksZ0_2Arms",300, -2000, 2000));

 PrimaryRecoTrksR02_2Arms  = std::auto_ptr<TH1D>(new TH1D("PrimaryRecoTrksR02_2Arms","PrimaryRecoTrksR02_2Arms",200,0,50.));
 PrimaryRecoTrksZ0_2Arms = std::auto_ptr<TH1D>(new TH1D("PrimaryRecoTrksZ0_2Arms","PrimaryRecoTrksZ0_2Arms",300, -2000, 2000));

PrimaryRecoTrksZ0_2ArmsMixedEnergy = std::auto_ptr<TH1D>(new TH1D("PrimaryRecoTrksZ0_2ArmsMixedEnergy","PrimaryRecoTrksZ0_2ArmsMixedEnergy",300, -2000, 2000));
PrimaryRecoTrksZ0_2ArmsEbig5 = std::auto_ptr<TH1D>(new TH1D("PrimaryRecoTrksZ0_2ArmsEbig5","PrimaryRecoTrksZ0_2ArmsEbig5",300, -2000, 2000));
PrimaryRecoTrksZ0_2ArmsEsmall5 = std::auto_ptr<TH1D>(new TH1D("PrimaryRecoTrksZ0_2ArmsEsmall5","PrimaryRecoTrksZ0_2ArmsEsmall5",300, -2000, 2000));




 H0_PrimaryRecoTrksZ0_2ArmsEbig5R0Less20= std::auto_ptr<TH1D>(new TH1D("H0_PrimaryRecoTrksZ0_2ArmsEbig5R0Less20","H0_PrimaryRecoTrksZ0_2ArmsEbig5R0Less20",300, -2000, 2000));
 H0_PrimaryRecoTrksEta_2ArmsEbig5 = std::auto_ptr<TH1D>(new TH1D("H0_PrimaryRecoTrksEta_2ArmsEbig5","H0_PrimaryRecoTrksEta_2ArmsEbig5",240, -8, 8));
 H0_PrimaryRecoTrksEta_2ArmsEbig5R0Less20 = std::auto_ptr<TH1D>(new TH1D("H0_PrimaryRecoTrksEta_2ArmsEbig5R0Less20","H0_PrimaryRecoTrksEta_2ArmsEbig5R0Less20",240, -8, 8));
 H0_PrimaryRecoTrksR0_2ArmsEbig5= std::auto_ptr<TH1D>(new TH1D("H0_PrimaryRecoTrksR0_2ArmsEbig5","H0_PrimaryRecoTrksR0_2ArmsEbig5",400,0,40.));


   H0_RecoTrksEta= std::auto_ptr<TH1D>(new TH1D("H0_RecoTrksEta","H0_RecoTrksEta",240, -8, 8));
   H0_SecondaryRecoTrksEta= std::auto_ptr<TH1D>(new TH1D("H0_SecondaryRecoTrksEta","H0_SecondaryRecoTrksEta",240, -8, 8));
   H0_SecondaryRecoTrksEtaR0Small20= std::auto_ptr<TH1D>(new TH1D("H0_SecondaryRecoTrksEtaR0Small20","H0_SecondaryRecoTrksEtaR0Small20",240, -8, 8));
   
   
   H0_RecoTrksEtaR0Small20= std::auto_ptr<TH1D>(new TH1D("H0_SecondaryRecoTrksEtaR0Small20","H0_SecondaryRecoTrksEtaR0Small20",240, -8, 8));

   H0_PrimaryRecoTrksBY_2Arms= std::auto_ptr<TH1D>(new TH1D("H0_PrimaryRecoTrksBY_2Arms","H0_PrimaryRecoTrksBY_2Arms",500, -250, 250));
   H0_PrimaryRecoTrksBX_2Arms= std::auto_ptr<TH1D>(new TH1D("H0_PrimaryRecoTrksBX_2Arms","H0_PrimaryRecoTrksBX_2Arms",500, -250, 250));  

   H0_SecondaryRecoTrksBX_2Arms= std::auto_ptr<TH1D>(new TH1D("H0_SecondaryRecoTrksBX_2Arms","H0_SecondaryRecoTrksBX_2Arms",500, -250, 250));
   H0_SecondaryRecoTrksBY_2Arms= std::auto_ptr<TH1D>(new TH1D("H0_SecondaryRecoTrksBY_2Arms","H0_SecondaryRecoTrksBY_2Arms",500, -250, 250));

   H0_PrimaryRecoTrksZ0_2ArmsEbig5 = std::auto_ptr<TH1D>(new TH1D("H0_PrimaryRecoTrksZ0_2ArmsEbig5","H0_PrimaryRecoTrksZ0_2ArmsEbig5",300, -2000, 2000));
   H0_PrimaryRecoTrksZ0_2Arms = std::auto_ptr<TH1D>(new TH1D("H0_PrimaryRecoTrksZ0_2Arms","H0_PrimaryRecoTrksZ0_2Arms",300, -2000, 2000));

   RecoTrkPrimaryNumHitEnergy = std::auto_ptr<TH1D>(new TH1D("RecoTrkPrimaryNumHitEnergy","RecoTrkPrimaryNumHitEnergy",21, -0.5, 20.5));
   RecoTrkPrimaryMostProbEnergy = std::auto_ptr<TH1D>(new TH1D("RecoTrkPrimaryMostProbEnergy","RecoTrkPrimaryMostProbEnergy",1000, -1.5, 998.5));

   H0_PrimaryXY_OrtogImpactEbig5=  std::auto_ptr<TH2D>(new TH2D("TrkXYImpact_H0", "TrkXYImpact_H0", 500, -200, 200, 500, -200, 200)) ;
   H0_PrimaryXY_OrtogImpact=  std::auto_ptr<TH2D>(new TH2D("TrkXYImpact_H0", "TrkXYImpact_H0", 500, -200, 200, 500, -200, 200)); 
   H0_XY_OrtogImpact=  std::auto_ptr<TH2D>(new TH2D("TrkXYImpact_H0", "TrkXYImpact_H0", 500, -200, 200, 500, -200, 200));


   H0_PrimaryZ0_OrtogImpact= std::auto_ptr<TH1D>(new TH1D("H0_PrimaryZ0_OrtogImpact","H0_PrimaryZ0_OrtogImpact",300, -2000, 2000));
   H0_Z0_OrtogImpact= std::auto_ptr<TH1D>(new TH1D("H0_Z0_OrtogImpact","H0_Z0_OrtogImpact",300, -2000, 2000));
   H0_RecoTrksZ0= std::auto_ptr<TH1D>(new TH1D("RecoTrksZ0","RecoTrksZ0",300, -2000, 2000));

   H0_PrimaryZ0_OrtogImpactEbig5  = std::auto_ptr<TH1D>(new TH1D("H0_PrimaryZ0_OrtogImpactEbig5","H0_PrimaryZ0_OrtogImpactEbig5",300, -2000, 2000));

  
 

   H0_PrimaryXY_OrtogImpact_BisMatch=  std::auto_ptr<TH2D>(new TH2D("H0_PrimaryXY_OrtogImpact_BisMatch", "H0_PrimaryXY_OrtogImpact_BisMatch", 500, -200, 200, 500, -200, 200)) ;
   H0_PrimaryZ0_OrtogImpact_BisMatch= std::auto_ptr<TH1D>(new TH1D("H0_PrimaryZ0_OrtogImpact_BisMatch","H0_PrimaryZ0_OrtogImpact_BisMatch",300, -2000, 2000));

   H0_PrimaryRecoTrksZ0_2Arms_BisMatch= std::auto_ptr<TH1D>(new TH1D("H0_PrimaryRecoTrksZ0_2Arms_BisMatch","H0_PrimaryRecoTrksZ0_2Arms_BisMatch",300, -2000, 2000));
   //   H0_PrimaryZ0_BisMatch= std::auto_ptr<TH1D>(new TH1D("H0_PrimaryZ0_BisMatch","H0_PrimaryZ0_BisMatch",300, -2000, 2000));


   PrimaryRecoTrksEta = std::auto_ptr<TH1D>(new TH1D("PrimaryRecoTrksEta","PrimaryRecoTrksEta",240, -8, 8));
 
   PrimaryRecoTrksR0_2Arms = std::auto_ptr<TH1D>(new TH1D("PrimaryRecoTrksR0_2Arms","PrimaryRecoTrksR0_2Arms",400,0,40.));

   SecondaryRecoTrksR0_2Arms= std::auto_ptr<TH1D>(new TH1D("SecondaryRecoTrksR0_2Arms","SecondaryRecoTrksR0_2Arms",400,0,40.));

 RecoTrksR0ALL= std::auto_ptr<TH1D>(new TH1D("RecoTrksR0ALL","RecoTrksR0ALL",400,0,40.));


 RecoTrksZ0ALL= std::auto_ptr<TH1D>(new TH1D("RecoTrksZ0ALL","RecoTrksZ0ALL",300, -2000, 2000));
  
 SecondaryRecoTrksZ0_2Arms_R0Cut20= std::auto_ptr<TH1D>(new TH1D("SecondaryRecoTrksZ0_2Arms_R0Cut20","SecondaryRecoTrksZ0_2Arms_R0Cut20",300, -2000, 2000));
PrimaryRecoTrksZ0_2Arms_EtaCut4_R0Cut20= std::auto_ptr<TH1D>(new TH1D("PrimaryRecoTrksZ0_2Arms_EtaCut4_R0Cut20","PrimaryRecoTrksZ0_2Arms_EtaCut4_R0Cut20",300, -2000, 2000));

 RecoTrksZ0ALL_R0Cut20= std::auto_ptr<TH1D>(new TH1D("RecoTrksZ0ALL_R0Cut20","RecoTrksZ0ALL_R0Cut20",300, -2000, 2000));
  
 RecoTrksZ0ALL_R0Cut10= std::auto_ptr<TH1D>(new TH1D("RecoTrksZ0ALL_R0Cut10","RecoTrksZ0ALL_R0Cut10",300, -2000, 2000));

 RecoTrksZ0ALL_EtaCut4 = std::auto_ptr<TH1D>(new TH1D("RecoTrksZ0ALL_EtaCut4","RecoTrksZ0ALL_EtaCut4",300, -2000, 2000));

 SecondarySmallZ0_EtavsR0= std::auto_ptr<TH2D>(new TH2D("SecondarySmallZ0_EtavsR0","SecondarySmallZ0_EtavsR0",240, -8, 8,400,0,40));


 SecondaryRecoTrksZ0_2Arms_EtaCut4_R0Cut20= std::auto_ptr<TH1D>(new TH1D("SecondaryRecoTrksZ0_2Arms_EtaCut4_R0Cut20","SecondaryRecoTrksZ0_2Arms_EtaCut4_R0Cut20",300, -2000, 2000));

 RecoTrksZ0ALL_EtaCut4_R0Cut20= std::auto_ptr<TH1D>(new TH1D("RecoTrksZ0ALL_EtaCut4_R0Cut20","RecoTrksZ0ALL_EtaCut4_R0Cut20",300, -2000, 2000));

 SecondaryBigZ0_EtavsR0= std::auto_ptr<TH2D>(new TH2D("SecondaryBigZ0_EtavsR0","SecondaryBigZ0_EtavsR0",240, -8, 8,400,0,40));

 
 RecoTrksR0ALL_EtaCut4_R0Cut20= std::auto_ptr<TH1D>(new TH1D("RecoTrksR0ALL_EtaCut4_R0Cut20","RecoTrksR0ALL_EtaCut4_R0Cut20",400,0,40.)); 
 SecondaryRecoTrksR0_2Arms_EtaCut4_R0Cut20= std::auto_ptr<TH1D>(new TH1D("SecondaryRecoTrksR0_2Arms_EtaCut4_R0Cut20","SecondaryRecoTrksR0_2Arms_EtaCut4_R0Cut20",400,0,40.)); 
 PrimaryRecoTrksR0_2Arms_EtaCut4_R0Cut20= std::auto_ptr<TH1D>(new TH1D("PrimaryRecoTrksR0_2Arms_EtaCut4_R0Cut20","PrimaryRecoTrksR0_2Arms_EtaCut4_R0Cut20",400,0,40.));

 

 NumRecoSecondary_VsNumRecoPrimary=  std::auto_ptr<TH2D>(new TH2D("NumRecoSecondary_VsNumRecoPrimary","NumRecoSecondary_VsNumRecoPrimary",20, -0.5, 19.5,20, -0.5, 19.5)) ; 
 NumRecoSecondary_VsNumRecoPrimary->SetXTitle("Primary");
 NumRecoSecondary_VsNumRecoPrimary->SetYTitle("Secondary");

 NumRecoSecondary_VsNumRecoPrimary_EtaCut4_R0Cut20=  std::auto_ptr<TH2D>(new TH2D("NumRecoSecondary_VsNumRecoPrimary_EtaCut4_R0Cut20","NumRecoSecondary_VsNumRecoPrimary_EtaCut4_R0Cut20",20, -0.5, 19.5,20, -0.5, 19.5)) ;


 simhitUnitID= std::auto_ptr<TH1D>(new TH1D("simhitUnitID","simhitUnitID",2000, 100000,300000)) ; 
 rechitUnitID= std::auto_ptr<TH1D>(new TH1D("rechitUnitID","rechitUnitID",2000, 100000,300000)) ; 
 FractionOfPrimaryTrksInSecVtx = std::auto_ptr<TH1D>(new TH1D("FractionOfPrimaryTrksInSecVtx","FractionOfPrimaryTrksInSecVtx",20, 0, 1.1)) ; 
 
 FractionOfPrimaryTrksRecoAsSecondary = std::auto_ptr<TH1D>(new TH1D("FractionOfPrimaryTrksRecoAsSecondary","FractionOfPrimaryTrksRecoAsSecondary",20, 0, 1.1)) ; 


 FractionOfSecondaryTrksInPrimaryVtx = std::auto_ptr<TH1D>(new TH1D("FractionOfSecondaryTrksInPrimaryVtx","FractionOfSecTrksInPrimaryVtx",20, 0, 1.1)) ; 
 
 FractionOfSecondaryTrksRecoAsPrimary = std::auto_ptr<TH1D>(new TH1D("FractionOfSecondaryTrksRecoAsPrimary","FractionOfSecondaryTrksRecoAsPrimary",20, 0, 1.1)) ; 

 NumRecoPrimVtx= std::auto_ptr<TH1D>(new TH1D("NumRecoPrimVtx","NumRecoPrimVtx",10, -0.5, 9.5)) ; 

 Purity_andEfficiency= std::auto_ptr<TH1D>(new TH1D("Purity_andEfficiency","(-1):Geant Trk with <3 Hits. (0):Reco-Fail, (1):1 Reco-Trk 2:More Than One Reco-Trk",12, -1.5, 11.5)) ;
 
 Purity_andEfficiencyNoGhost= std::auto_ptr<TH1D>(new TH1D("Purity_andEfficiencyNoGhost","(-1):Geant Trk with <3 Hits. (0):Reco-Fail, (1):1 Reco-Trk 2:More Than One Reco-Trk",12, -1.5, 11.5)) ;

PrimaryGeanteta2= std::auto_ptr<TH1D> (new TH1D("PrimaryGeanteta2","Geant Primary Track #eta",150,-maxetafordndetahisto,maxetafordndetahisto));



 FakeSecondaryRecoDnDeta2FromVtx= std::auto_ptr<TProfile> (new TProfile("FakeSecondaryRecoDnDeta2FromVtx"," DN/D#eta of Reconstructed Secondary Trk Pointing to the Vtx vs #eta",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 FakeSecondaryRecoDnDeta2FromVtx->SetYTitle("dN/d#eta");
 FakeSecondaryRecoDnDeta2FromVtx->SetXTitle("#eta");

 PrimaryRecoDnDeta2NoCut= std::auto_ptr<TProfile> (new TProfile("PrimaryRecoDnDeta2NoCut"," DN/D#eta of Reconstructed Primary Trk (No cuts) vs #eta",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 PrimaryRecoDnDeta2NoCut->SetYTitle("dN/d#eta");
 PrimaryRecoDnDeta2NoCut->SetXTitle("#eta");

 
 PrimaryRecoDnDeta2FromVtx= std::auto_ptr<TProfile> (new TProfile("PrimaryRecoDnDeta2FromVtx"," DN/D#eta of Reconstructed Primary Trk Pointing to the Vtx vs #eta",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 PrimaryRecoDnDeta2FromVtx->SetYTitle("dN/d#eta");
 PrimaryRecoDnDeta2FromVtx->SetXTitle("#eta");

 RecoDnDeta2FromVtx= std::auto_ptr<TProfile> (new TProfile("RecoDnDeta2FromVtx","Reconstructed track DN/D#eta (Trk Pointing to the Vtx)",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 RecoDnDeta2FromVtx->SetYTitle("dN/d#eta");
 RecoDnDeta2FromVtx->SetXTitle("#eta");
  
RecoDnDeta2NoCut= std::auto_ptr<TProfile> (new TProfile("RecoDnDeta2NoCut","Reconstructed track DN/D#eta ",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 RecoDnDeta2NoCut->SetYTitle("dN/d#eta");
 RecoDnDeta2NoCut->SetXTitle("#eta");

 PrimaryGeantDNDeta2_IfOneRecoZImpInH0= std::auto_ptr<TProfile> (new TProfile("PrimaryGeantDNDeta2_IfOneRecoZImpInH0","Geant Primary Track DN/D#eta",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 PrimaryGeantDNDeta2_IfOneRecoZImpInH0->SetYTitle("dN/d#eta");
 PrimaryGeantDNDeta2_IfOneRecoZImpInH0->SetXTitle("#eta");


 sprintf(sZnamehist, "H0 Geant Primary Track DN/D#eta Ecut %f at least 1 in T2", PtCutinPrimaryEfficiency/*EnergyCutinPrimaryEfficiency*/); 
 PrimaryGeantDNDeta2_IfOneTrkInT2_H0= std::auto_ptr<TProfile> (new TProfile("PrimaryGeantDNDeta2_IfOneTrkInT2_H0",sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

 sprintf(sZnamehist, "H1 Geant Primary Track DN/D#eta Ecut %f at least 1 in T2", PtCutinPrimaryEfficiency/*EnergyCutinPrimaryEfficiency*/); 
 PrimaryGeantDNDeta2_IfOneTrkInT2_H1= std::auto_ptr<TProfile> (new TProfile("PrimaryGeantDNDeta2_IfOneTrkInT2_H1",sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

PrimaryGeantDNDeta2_IfOneTrkInT2_H1sometrkexcluding64= std::auto_ptr<TProfile> (new TProfile("PrimaryGeantDNDeta2_IfOneTrkInT2_H1sometrkexcluding64",sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

 sprintf(sZnamehist, "H0 Geant Primary Track DN/D#eta Ecut %f", PtCutinPrimaryEfficiency/*EnergyCutinPrimaryEfficiency*/); 
 PrimaryGeantDNDeta2= std::auto_ptr<TProfile> (new TProfile("PrimaryGeantDNDeta2",sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 PrimaryGeantDNDeta2->SetYTitle("dN/d#eta");
 PrimaryGeantDNDeta2->SetXTitle("#eta");


 PrimaryGeantDNDeta2H0_ZCut= std::auto_ptr<TProfile> (new TProfile("PrimaryGeantDNDeta2H0_ZCut","Geant Primary Track DN/D#eta Z cutted",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 PrimaryGeantDNDeta2H0_ZCut->SetYTitle("dN/d#eta");
 PrimaryGeantDNDeta2H0_ZCut->SetYTitle("#eta");

 
 sprintf(sZnamehist, "ChP  DN/D#eta Ecut %f at least 1 in T2", PtCutinPrimaryEfficiency/*EnergyCutinPrimaryEfficiency*/); 
 DNDetaMBALLMCGenerator_GeneratorTriggered= std::auto_ptr<TProfile> (new TProfile("DNDetaMBALLMCGenerator_GeneratorTriggered",sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));


 DNDetaMBALLMCGenerator_GeneratorTriggeredAtLeastOne5364= std::auto_ptr<TProfile> (new TProfile("DNDetaMBALLMCGenerator_GeneratorTriggeredAtLeastOne5364",sZnamehist,totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

 RecoDnDeta2UnfoldCut2= std::auto_ptr<TProfile> (new TProfile("RecoDnDeta2UnfoldCut2","Track DN/D#eta (UnfoldCut2)",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 RecoDnDeta2UnfoldCut2->SetYTitle("dN/d#eta");
 RecoDnDeta2UnfoldCut2->SetXTitle("#eta");


 RecoDnDeta2UnfoldCut2_Primary= std::auto_ptr<TProfile> (new TProfile("RecoDnDeta2UnfoldCut2_Primary","Track DN/D#eta (UnfoldCut2), primary",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 RecoDnDeta2UnfoldCut2_Primary->SetYTitle("dN/d#eta");
 RecoDnDeta2UnfoldCut2_Primary->SetXTitle("#eta");

  RecoDnDeta2UnfoldCut= std::auto_ptr<TProfile> (new TProfile("RecoDnDeta2UnfoldCut","Track DN/D#eta (UnfoldCut)",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 RecoDnDeta2UnfoldCut->SetYTitle("dN/d#eta");
 RecoDnDeta2UnfoldCut->SetXTitle("#eta");


 RecoDnDeta2UnfoldCut_Primary= std::auto_ptr<TProfile> (new TProfile("RecoDnDeta2UnfoldCut_Primary","Track DN/D#eta (UnfoldCut), primary",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
 RecoDnDeta2UnfoldCut_Primary->SetYTitle("dN/d#eta");
 RecoDnDeta2UnfoldCut_Primary->SetXTitle("#eta");


 TrkRecoEfficiencyProfilePlus= std::auto_ptr<TProfile> (new TProfile("TrkRecoEfficiencyProfilePlus","Trk Reco Efficiency as a function of # Geant Trk (Arm +)",25,-0.5,24.5));

 TrkRecoEfficiencyProfileMinus= std::auto_ptr<TProfile> (new TProfile("TrkRecoEfficiencyProfileMinus","Trk Reco Efficiency as a function of # Geant Trk (Arm -)",25,-0.5,24.5));

 GeantTrkMultiplicityfrequencyPlus= std::auto_ptr<TH1D>(new TH1D("GeantTrkMultiplicityfrequencyPlus","GeantTrkMultiplicityfrequencyPlus",30, -0.5, 29.5)) ; 


H0_TrkEtaEfficiencyCumulative= std::auto_ptr<TProfile> (new TProfile("H0_TrkEtaEfficiencyCumulative","H0_TrkEtaEfficiencyCumulative",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
H1_TrkEtaEfficiencyCumulative= std::auto_ptr<TProfile> (new TProfile("H1_TrkEtaEfficiencyCumulative","H1_TrkEtaEfficiencyCumulative",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
H2_TrkEtaEfficiencyCumulative= std::auto_ptr<TProfile> (new TProfile("H2_TrkEtaEfficiencyCumulative","H2_TrkEtaEfficiencyCumulative",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
H3_TrkEtaEfficiencyCumulative= std::auto_ptr<TProfile> (new TProfile("H3_TrkEtaEfficiencyCumulative","H3_TrkEtaEfficiencyCumulative",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

H0_TrkEtaEfficiencyCumulative_UnfoldCuts= std::auto_ptr<TProfile> (new TProfile("H0_TrkEtaEfficiencyCumulative_UnfoldCuts","H0_TrkEtaEfficiencyCumulative_UnfoldCuts",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
H1_TrkEtaEfficiencyCumulative_UnfoldCuts= std::auto_ptr<TProfile> (new TProfile("H1_TrkEtaEfficiencyCumulative_UnfoldCuts","H0_TrkEtaEfficiencyCumulative_UnfoldCuts",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
H2_TrkEtaEfficiencyCumulative_UnfoldCuts= std::auto_ptr<TProfile> (new TProfile("H2_TrkEtaEfficiencyCumulative_UnfoldCuts","H0_TrkEtaEfficiencyCumulative_UnfoldCuts",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
H3_TrkEtaEfficiencyCumulative_UnfoldCuts= std::auto_ptr<TProfile> (new TProfile("H3_TrkEtaEfficiencyCumulative_UnfoldCuts","H0_TrkEtaEfficiencyCumulative_UnfoldCuts",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));



//H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult




H0_TrkEtaEfficiencyCumulative_UnfoldCuts_2= std::auto_ptr<TProfile> (new TProfile("H0_TrkEtaEfficiencyCumulative_UnfoldCuts_2","H0_TrkEtaEfficiencyCumulative_UnfoldCuts_2",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
H1_TrkEtaEfficiencyCumulative_UnfoldCuts_2= std::auto_ptr<TProfile> (new TProfile("H1_TrkEtaEfficiencyCumulative_UnfoldCuts_2","H0_TrkEtaEfficiencyCumulative_UnfoldCuts_2",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
H2_TrkEtaEfficiencyCumulative_UnfoldCuts_2= std::auto_ptr<TProfile> (new TProfile("H2_TrkEtaEfficiencyCumulative_UnfoldCuts_2","H0_TrkEtaEfficiencyCumulative_UnfoldCuts_2",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));
H3_TrkEtaEfficiencyCumulative_UnfoldCuts_2= std::auto_ptr<TProfile> (new TProfile("H3_TrkEtaEfficiencyCumulative_UnfoldCuts_2","H0_TrkEtaEfficiencyCumulative_UnfoldCuts_2",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));


 ArmMinus_TrkEtaEfficiencyCumulative= std::auto_ptr<TProfile> (new TProfile("ArmMinus_TrkEtaEfficiencyCumulative","ArmMinus_TrkEtaEfficiencyCumulative",100,-0.5,99.5));
 ArmPlus_TrkEtaEfficiencyCumulative= std::auto_ptr<TProfile> (new TProfile("ArmPlus_TrkEtaEfficiencyCumulative","ArmPlus_TrkEtaEfficiencyCumulative",100,-0.5,99.5));
  
 ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts= std::auto_ptr<TProfile> (new TProfile("ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts","ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

 ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts= std::auto_ptr<TProfile> (new TProfile("ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts","ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

   
 ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts_2= std::auto_ptr<TProfile> (new TProfile("ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts_2","ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts_2",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

 ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts_2= std::auto_ptr<TProfile> (new TProfile("ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts_2","ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts_2",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));  



ArmPlus_TrkEtaEfficiencyCumulativeSecondary= std::auto_ptr<TProfile> (new TProfile("ArmPlus_TrkEtaEfficiencyCumulativeSecondary","ArmPlus_TrkEtaEfficiencyCumulativeSecondary",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

ArmMinus_TrkEtaEfficiencyCumulativeSecondary= std::auto_ptr<TProfile> (new TProfile("ArmMinus_TrkEtaEfficiencyCumulativeSecondary","ArmMinus_TrkEtaEfficiencyCumulativeSecondary",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCutsSecondary= std::auto_ptr<TProfile> (new TProfile("ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCutsSecondary","ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCutsSecondary",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCutsSecondary= std::auto_ptr<TProfile> (new TProfile("ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCutsSecondary","ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCutsSecondary",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts_2Secondary= std::auto_ptr<TProfile> (new TProfile("ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts_2Secondary","ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts_2Secondary",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));

ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts_2Secondary= std::auto_ptr<TProfile> (new TProfile("ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts_2Secondary","ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts_2Secondary",totnumberbinfordndeta,-maxetafordndetahisto,maxetafordndetahisto));


Geant_Prim_Or_Sec_Unfold= std::auto_ptr<TH1D> (new TH1D("Geant_Prim_Or_Sec_Unfold","Geant Trk Info 1:Primary 2:Secondary 3:Mix",4, -0.5, 3.5));

     ZImpOfAPrimaryGeantTrk_Cumulative= std::auto_ptr<TH1D> (new TH1D("ZImpOfAPrimaryGeantTrk_Cumulative","ZImpOfAPrimaryGeantTrk_Cumulative",3000, -15000, 15000));
     ZImpOfAPrimaryGeantTrk_CumulativeZMinCut= std::auto_ptr<TH1D> (new TH1D("ZImpOfAPrimaryGeantTrk_CumulativeZMinCut","ZImpOfAPrimaryGeantTrk_CumulativeZMinCut",3000, -15000, 15000));
     ZImpOfAPrimaryGeantTrk_CumulativeZMinCut_match= std::auto_ptr<TH1D> (new TH1D("ZImpOfAPrimaryGeantTrk_CumulativeZMinCut_match","ZImpOfAPrimaryGeantTrk_CumulativeZMinCut_match",3000, -15000, 15000));
     ZImpOfAPrimaryGeantTrk_Cumulative_match= std::auto_ptr<TH1D> (new TH1D("ZImpOfAPrimaryGeantTrk_Cumulative_match","ZImpOfAPrimaryGeantTrk_Cumulative_match",3000, -15000, 15000));

EtaResolGeant=std::auto_ptr<TProfile>(new TProfile("EtaResolGeant","EtaResolGeant",20,5.,7.,"")); 
EtaResolGenerator=std::auto_ptr<TProfile>(new TProfile("EtaResolGenerator","EtaResolGenerator",20,5.,7.,"")); 
EtaResolFreeGeant=std::auto_ptr<TProfile>(new TProfile("EtaResolFreeGeant","EtaResolFreeGeant",20,5.,7.,"")); 
EtaResolFreeGenerator=std::auto_ptr<TProfile>(new TProfile("EtaResolFreeGenerator","EtaResolFreeGenerator",20,5.,7.,"")); 





 H0TrkEffi_AxisAvgMult_ZImpactCut=  std::auto_ptr<TProfile>(new TProfile("H0TrkEffi_AxisAvgMult_ZImpactCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 
 H1TrkEffi_AxisAvgMult_ZImpactCut=  std::auto_ptr<TProfile>(new TProfile("H1TrkEffi_AxisAvgMult_ZImpactCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 
 H2TrkEffi_AxisAvgMult_ZImpactCut=  std::auto_ptr<TProfile>(new TProfile("H2TrkEffi_AxisAvgMult_ZImpactCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 
 H3TrkEffi_AxisAvgMult_ZImpactCut=  std::auto_ptr<TProfile>(new TProfile("H3TrkEffi_AxisAvgMult_ZImpactCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,"")); 

 H0TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut=  std::auto_ptr<TProfile>(new TProfile("H0TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,""));
 H1TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut=  std::auto_ptr<TProfile>(new TProfile("H1TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,""));
 H2TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut=  std::auto_ptr<TProfile>(new TProfile("H2TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,""));
 H3TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut=  std::auto_ptr<TProfile>(new TProfile("H3TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,""));


 H0TrkEffi_AxisAvgMult_NORecCut_GeantZCut=  std::auto_ptr<TProfile>(new TProfile("H0TrkEffi_AxisAvgMult_NORecCut_GeantZCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,""));
 H1TrkEffi_AxisAvgMult_NORecCut_GeantZCut=  std::auto_ptr<TProfile>(new TProfile("H1TrkEffi_AxisAvgMult_NORecCut_GeantZCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,""));
 H2TrkEffi_AxisAvgMult_NORecCut_GeantZCut=  std::auto_ptr<TProfile>(new TProfile("H2TrkEffi_AxisAvgMult_NORecCut_GeantZCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,""));
  H3TrkEffi_AxisAvgMult_NORecCut_GeantZCut=  std::auto_ptr<TProfile>(new TProfile("H3TrkEffi_AxisAvgMult_NORecCut_GeantZCut","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,""));
 
H0TrkEffi_AxisAvgMult_NORecCut_GeantZCut_MatchProj=  std::auto_ptr<TProfile>(new TProfile("H0TrkEffi_AxisAvgMult_NORecCut_GeantZCut_MatchProj","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,""));
H0TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut_MatchProj=  std::auto_ptr<TProfile>(new TProfile("H0TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut_MatchProj","<# matching Tracks HTest / # Tracks HRef_raw > vs # Tracks HRef_raw. ",20,-0.5,59.5,""));


 Cumulative_GeantTrk_ZImpactRefQ= std::auto_ptr<TH1D> (new TH1D("Cumulative_GeantTrk_ZImpactRefQ","Geant Trk ZImpact for reference quarter",1500, -15000, 15000));
 

 H0H1_PhiGeantTrkWhenEffiFail= std::auto_ptr<TH1D> (new TH1D("H0H1_PhiGeantTrkWhenEffiFail","--H0H1_PhiGeantTrkWhenEffiFailMultLess20",360, -0.5, 359.5));

 H0H1_EtaPhiGeantTrkWhenEffiFail= std::auto_ptr<TH2D> (new TH2D("H0H1_EtaPhiGeantTrkWhenEffiFail","H0H1_EtaPhiGeantTrkWhenEffiFail (H0)",90, -0.5, 359.5,21, 4.5, 8.5));


 H0_EtaPhiGeantTrkWhenEffiFailMultLess20= std::auto_ptr<TH2D> (new TH2D("H0_EtaPhiGeantTrkWhenEffiFailMultLess20","H0_EtaPhiGeantTrkWhenEffiFailMultLess20 (H0)",21, 4.5, 8.5,90, -0.5, 359.5));

 Numtrk_vs_NumtrkpasCutsWhenEffiFail= std::auto_ptr<TH2D> (new TH2D("Numtrk_vs_NumtrkpasCutsWhenEffiFail","Numtrk_vs_NumtrkpasCutsWhenEffiFail (H0?)",21, -0.5, 20.5,21, -0.5, 20.5));

 NumH0PlaneActiveWhenGeantEffiFail= std::auto_ptr<TH1D> (new TH1D("NumH0PlaneActiveWhenGeantEffiFail","NumH0PlaneActiveWhenGeantEffiFail",21, -0.5, 20.5));
 NumH0TrkWhenGeantEffiFail= std::auto_ptr<TH1D> (new TH1D("NumH0TrkWhenGeantEffiFail","NumH0TrkWhenGeantEffiFail",21, -0.5, 20.5));
 NumH0PlaneActive_vsMultipl_WhenGeantEffiFail= std::auto_ptr<TH2D> (new TH2D("NumH0PlaneActive_vsMultipl_WhenGeantEffiFail","NumH0PlaneActive_vsMultipl_WhenGeantEffiFail",21, -0.5, 20.5,56, -0.5, 55.5));
minDxDyWhenGeantEffiFailH1= std::auto_ptr<TH2D>(new TH2D("minDxDyWhenGeantEffiFailH1","minDxDyWhenGeantEffiFail (H1)",241, -120.5, 120.5,241, -120.5, 120.5)); 
 minDxDyWhenGeantEffiFail= std::auto_ptr<TH2D>(new TH2D("minDxDyWhenGeantEffiFail","minDxDyWhenGeantEffiFail (H0)",241, -120.5, 120.5,241, -120.5, 120.5)); 
 //minDphiWhenGeantEffiFail= std::auto_ptr<TH1D> (new TH1D("minDphiWhenGeantEffiFail","minDphiWhenGeantEffiFail",100, -0.5, 90.5));
 //minDrWhenGeantEffiFail= std::auto_ptr<TH1D> (new TH1D("minDrWhenGeantEffiFail","minDrWhenGeantEffiFail",21, -0.5, 20.5));
 ZImpOfARecTrkWhenEffiFail = std::auto_ptr<TH1D> (new TH1D("ZImpOfARecTrkWhenEffiFail","ZImpOfARecTrkWhenEffiFail",150, -15000, 15000));

 ZImpOfAGeantTrkWhenEffiFailH1= std::auto_ptr<TH1D> (new TH1D("ZImpOfAGeantTrkWhenEffiFailH1","ZImpOfAGeantTrkWhenEffiFailH1",150, -15000, 15000));
 ZImpOfAGeantTrkWhenEffiFailH0= std::auto_ptr<TH1D> (new TH1D("ZImpOfAGeantTrkWhenEffiFailH0","ZImpOfAGeantTrkWhenEffiFailH0",150, -15000, 15000));

 Cumulative_GeantTrk_ZImpactRefQ_PtAboveZ5= std::auto_ptr<TH1D> (new TH1D("Cumulative_GeantTrk_ZImpactRefQ_PtAboveZ5","Cumulative_GeantTrk_ZImpactRefQ_PtAboveZ5",1500, 0., 3000));
 Cumulative_GeantTrk_ZImpactRefQ_E= std::auto_ptr<TH1D> (new TH1D("Cumulative_GeantTrk_ZImpactRefQ_E","Cumulative_GeantTrk_ZImpactRefQ_E",1500, 0., 3000));
 Cumulative_GeantTrk_ZImpactRefQ_EAboveZ5= std::auto_ptr<TH1D> (new TH1D("Cumulative_GeantTrk_ZImpactRefQ_EAboveZ5","Cumulative_GeantTrk_ZImpactRefQ_EAboveZ5",1500, 0., 3000));

 numberOfPrimaryTracks = 0;
 numberOfSecondaryTracks = 0;



 PlaneClusterActiveCounter= std::auto_ptr<TH1D> (new TH1D("PlaneClusterActiveCounter","PlaneClusterActiveCounter",50, -0.5, 49.5));


 PrimTrkTestedQuarter= std::auto_ptr<TH1D> (new TH1D("PrimTrkTestedQuarter","PrimTrkTestedQuarter",4, -0.5, 3.5));
 PrimTrkTestedPlaneQuarter= std::auto_ptr<TH1D> (new TH1D("PrimTrkTestedPlaneQuarter","PrimTrkTestedPlaneQuarter",40, -0.5, 39.5));
 PrimTrkTestedQuarterUnamb= std::auto_ptr<TH1D> (new TH1D("PrimTrkTestedQuarterUnamb","PrimTrkTestedQuarterUnamb",4, -0.5, 3.5));

 

 eta2geant_forEffQ0= std::auto_ptr<TH1D> (new TH1D("eta2geant_forEffQ0","eta2geant_forEffQ0",400, -8.5, 8.5));
 eta2geant_forEffQ1= std::auto_ptr<TH1D> (new TH1D("eta2geant_forEffQ1","eta2geant_forEffQ1",400, -8.5, 8.5));
 eta2geant_forEffQ2= std::auto_ptr<TH1D> (new TH1D("eta2geant_forEffQ2","eta2geant_forEffQ2",400, -8.5, 8.5));
 eta2geant_forEffQ3= std::auto_ptr<TH1D> (new TH1D("eta2geant_forEffQ3","eta2geant_forEffQ3",400, -8.5, 8.5));

 eta2geant_forEffQ0eff= std::auto_ptr<TH1D> (new TH1D("eta2geant_forEffQ0eff","eta2geant_forEffQ0eff",400, -8.5, 8.5));
 eta2geant_forEffQ1eff= std::auto_ptr<TH1D> (new TH1D("eta2geant_forEffQ1eff","eta2geant_forEffQ1eff",400, -8.5, 8.5));
 eta2geant_forEffQ2eff= std::auto_ptr<TH1D> (new TH1D("eta2geant_forEffQ2eff","eta2geant_forEffQ2eff",400, -8.5, 8.5));
 eta2geant_forEffQ3eff= std::auto_ptr<TH1D> (new TH1D("eta2geant_forEffQ3eff","eta2geant_forEffQ3eff",400, -8.5, 8.5));

 H0_BinMultStat= std::auto_ptr<TH1D> (new TH1D("H0_BinMultStat","H0_BinMultStat",400, -0.5, 48.5));
 H1_BinMultStat= std::auto_ptr<TH1D> (new TH1D("H1_BinMultStat","H1_BinMultStat",400, -0.5, 48.5));
 H2_BinMultStat= std::auto_ptr<TH1D> (new TH1D("H2_BinMultStat","H2_BinMultStat",400, -0.5, 48.5));
 H3_BinMultStat= std::auto_ptr<TH1D> (new TH1D("H3_BinMultStat","H3_BinMultStat",400, -0.5, 48.5));
 H0_BinMultStatFail= std::auto_ptr<TH1D> (new TH1D("H0_BinMultStatFail","H0_BinMultStatFail",400, -0.5, 48.5));
 H1_BinMultStatFail= std::auto_ptr<TH1D> (new TH1D("H1_BinMultStatFail","H1_BinMultStatFail",400, -0.5, 48.5));
 H2_BinMultStatFail= std::auto_ptr<TH1D> (new TH1D("H2_BinMultStatFail","H2_BinMultStatFail",400, -0.5, 48.5));
 H3_BinMultStatFail= std::auto_ptr<TH1D> (new TH1D("H3_BinMultStatFail","H3_BinMultStatFail",400, -0.5, 48.5));

 NumHitInRecTrkH0= std::auto_ptr<TH1D>(new TH1D("NumHitInRecTrkH0","NumHitInRecTrkH0", 21,-0.5, 20.5));
 NumHitInRecTrkH1= std::auto_ptr<TH1D>(new TH1D("NumHitInRecTrkH1","NumHitInRecTrkH1", 21,-0.5, 20.5));
 NumHitInRecTrkH2= std::auto_ptr<TH1D>(new TH1D("NumHitInRecTrkH2","NumHitInRecTrkH2", 21,-0.5, 20.5));
 NumHitInRecTrkH3= std::auto_ptr<TH1D>(new TH1D("NumHitInRecTrkH3","NumHitInRecTrkH3", 21,-0.5, 20.5));

 AveragePadCLSDistrH1= std::auto_ptr<TH1D>(new TH1D("AveragePadCLSDistrH1","AveragePadCLSDistrH1", 400,-0.5, 399.5));
 AveragePadCLSDistrH0= std::auto_ptr<TH1D>(new TH1D("AveragePadCLSDistrH0","AveragePadCLSDistrH0", 400,-0.5, 399.5));
   
 corposprimaryXYQ0=std::auto_ptr<TH2D>(new TH2D("corposprimaryXYQ0","corposprimaryXYQ0",290, -145, 145,290, -145, 145)); 
 corposprimaryXYQ1=std::auto_ptr<TH2D>(new TH2D("corposprimaryXYQ1","corposprimaryXYQ1",290, -145, 145,290, -145, 145)); 

 Q0PadRCPrimary=std::auto_ptr<TH2D>(new TH2D("Q0PadRCPrimary","Q0PadRCPrimary",50, -0.5, 49.5,100, -0.5, 99.5));
 Q1PadRCPrimary=std::auto_ptr<TH2D>(new TH2D("Q1PadRCPrimary","Q1PadRCPrimary",50, -0.5, 49.5,100, -0.5, 99.5));
     
 PadRowGeantRecoH0= std::auto_ptr<TProfile>(new TProfile("PadRowGeantRecoH0","PadRowGeantRecoH0",20,-0.5,59.5,""));//std::auto_ptr<TH2D>(new TH2D("PadRowGeantRecoH0","PadRowGeantRecoH0",50, -0.5, 49.5,100, -0.5, 99.5)); 
 PadRowGeantRecoH0_bis= std::auto_ptr<TProfile>(new TProfile("PadRowGeantRecoH0_bis","PadRowGeantRecoH0_bis",20,-0.5,59.5,""));//std::auto_ptr<TH2D>(new TH2D("PadRowGeantRecoH0_bis","PadRowGeantRecoH0_bis",50, -0.5, 49.5,100, -0.5, 99.5)); 
 PadRowGeantRecoH1= std::auto_ptr<TProfile>(new TProfile("PadRowGeantRecoH1","PadRowGeantRecoH1",20,-0.5,59.5,""));//std::auto_ptr<TH2D>(new TH2D("PadRowGeantRecoH1","PadRowGeantRecoH1",50, -0.5, 49.5,100, -0.5, 99.5)); 
 PadRowGeantRecoH1_bis= std::auto_ptr<TProfile>(new TProfile("PadRowGeantRecoH1_bis","PadRowGeantRecoH1_bis ",20,-0.5,59.5,""));//std::auto_ptr<TH2D>(new TH2D("PadRowGeantRecoH1_bis","PadRowGeantRecoH1_bis",50, -0.5, 49.5,100, -0.5, 99.5));
 
 RPhiGeantHitinH0Col32=std::auto_ptr<TH2D>(new TH2D("RPhiGeantHitinH0Col32","RPhiGeantHitinH0Col32",50, 0., 150.5,361, -0.5, 360.5));
 XYGeantHitinH0Col32=std::auto_ptr<TH2D>(new TH2D("XYGeantHitinH0Col32","XYGeantHitinH0Col32",100, -150., 150.,100, -150, 150));

  H1PrimaryTrkPdgID=std::auto_ptr<TH1D>(new TH1D("H1PrimaryTrkPdgID","H1PrimaryTrkPdgID",2000,-3500,3500)) ;
  H0PrimaryTrkPdgID=std::auto_ptr<TH1D>(new TH1D("H0PrimaryTrkPdgID","H0PrimaryTrkPdgID",2000,-3500,3500)) ;
  H0PrimaryTrkE=std::auto_ptr<TH1D>(new TH1D("H0PrimaryTrkE","H0PrimaryTrkE",300,-0.5,299.5)) ;
  H1PrimaryTrkE=std::auto_ptr<TH1D>(new TH1D("H1PrimaryTrkE","H1PrimaryTrkE",300,-0.5,299.5)) ;

  for(unsigned int m=0;m<=22; m++){
    //    (int)(fabs(the_eta2)-5.35)/0.05;//21 Bin
    


    sprintf(sZname2, "Associated_PrimaryTrkZImpact_H0 %d", m);
    sprintf(sZnamehist, "Associated_PrimaryTrkZImpact_H0-%d", m); 
    
   Associated_PrimaryTrkZImpact_H0[m]= std::auto_ptr<TH1D>(new TH1D(sZname2,sZnamehist,300,-15000., 15000.));


   sprintf(sZname2, "Associated_PrimaryTrkZImpact_H1 %d", m);
   sprintf(sZnamehist, "Associated_PrimaryTrkZImpact_H1-%d", m); 
 
   Associated_PrimaryTrkZImpact_H1[m]= std::auto_ptr<TH1D>(new TH1D(sZname2,sZnamehist,300,-15000., 15000));

   sprintf(sZname2, "Associated_PrimaryTrkZImpact_H2 %d", m);
   sprintf(sZnamehist, "Associated_PrimaryTrkZImpact_H2-%d", m); 
  
   Associated_PrimaryTrkZImpact_H2[m]= std::auto_ptr<TH1D>(new TH1D(sZname2,sZnamehist,300,-15000., 15000));

   sprintf(sZname2, "Associated_PrimaryTrkZImpact_H3 %d", m);
   sprintf(sZnamehist, "Associated_PrimaryTrkZImpact_H3-%d", m); 
 
   Associated_PrimaryTrkZImpact_H3[m]= std::auto_ptr<TH1D>(new TH1D(sZname2,sZnamehist,300,-15000., 15000));
}


 for(unsigned int m=0;m<40; m++){

   sprintf(sZname2, "PlaneGeantPadColXY %d", m);
   sprintf(sZnamehist, "PlaneGeantPadColXY-%d", m); 
   PlaneGeantPadColXY[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,70,-0.5,69.5,70,-0.5,69.5));

   sprintf(sZname2, "PlaneGeantToRecXY %d", m);
   sprintf(sZnamehist, "PlaneGeantToRecXY-%d", m); 
   PlaneGeantToRecXY[m]= std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,300,-150.,150.,300,-140.,140.));

 }
 
 CutZimpLeftH0= std::auto_ptr<TProfile>(new TProfile("CutZimpLeftH0","CutZimpLeftH0",100,5.,7.));
 CutZimpRightH0= std::auto_ptr<TProfile>(new TProfile("CutZimpRightH0","CutZimpRightH0",100,5.,7.));

 CutZimpLeftH1= std::auto_ptr<TProfile>(new TProfile("CutZimpLeftH1","CutZimpLeftH1",100,5.,7.));
 CutZimpRightH1= std::auto_ptr<TProfile>(new TProfile("CutZimpRightH1","CutZimpRightH1",100,5.,7.));

 CutZimpLeftH2= std::auto_ptr<TProfile>(new TProfile("CutZimpLeftH2","CutZimpLeftH2",100,5.,7.));
 CutZimpRightH2= std::auto_ptr<TProfile>(new TProfile("CutZimpRightH2","CutZimpRightH2",100,5.,7.));

 CutZimpLeftH3= std::auto_ptr<TProfile>(new TProfile("CutZimpLeftH3","CutZimpLeftH3",100,5.,7.));
 CutZimpRightH3= std::auto_ptr<TProfile>(new TProfile("CutZimpRightH3","CutZimpRightH3",100,5.,7.));





}

// ------------ method called once each job just after ending the event loop  ------------
void T2BackgroundAn::endJob()
{
  if(verbosity>0)
    std::cout<<"Begin endjob"<<std::endl;
  T2Geometry t2geompl;
  T2GeometryUtil convpl;  
  /*
	double X=PrimaryGeantHitPos.at(i).x();
	double Y=PrimaryGeantHitPos.at(i).y();
	double Z=theZ.at(i);
	double EX=2.;
	double EY=2.;
	double EZ=1.;
	uint32_t mycmsswid=thecmsswId.at(i);
	T2Hit arechit(X,  Y,  Z,  EX,  EY, EZ, mycmsswid);*/
  /*
  for (unsigned int symb=0; symb<40; symb++){
    T2GeometryUtil::T2DetInfo detInfo = convpl.GetT2Info(symb);
    t2geompl.setPlane(detInfo.cmsswid);
    for(int xx=-150;xx<150; xx++)
      for(int yy=-150;yy<150; yy++){
	//Hit creation
	double zz=detInfo.Zdet;	
	double EX=2.;
	double EY=2.;
	double EZ=1.;
	T2Hit arechit(xx,  yy,  zz,  EX,  EY, EZ, detInfo.cmsswid);
	PlaneGeantToRecXY[symb]->Fill(arechit.GetHitX(),arechit.GetHitY());
	Local3DPoint corPos(xx,yy,zz);
	int padRow = t2geompl.getNearestPadRow2_(&corPos,detInfo.cmsswid);
	int padCol = t2geompl.getNearestPadCol_(arechit.GetHitR(),arechit.GetHitPhi(),detInfo.cmsswid);//getNearestPadCol
	PlaneGeantPadColXY[symb]->Fill(padRow,padCol);
      }
  } 
*/
 
  /*
  //Loop 0-40
  //Convert to cmsswid
  //loop -150+150mm XY. Create Hit(x,y,z)
  T2Hit hit = recTrack->GetHitT2(h);
  int detUnitId = hit.GetHitDetRawId();

  rechitUnitID->Fill(hit.GetHitDetRawId()/10000);
  T2GeometryUtil::T2DetInfo detInfo2 = conv.GetT2Info(detUnitId);
      
  t2geom.setPlane(detUnitId);
  Local3DPoint corPos(hit.GetHitX(),hit.GetHitY(),hit.GetHitZ());
  int padRow = t2geom.getNearestPadRow(&corPos);
  int padCol = t2geom.getNearestPadCol(&corPos);
  */



  
  if (verbosity > 1){
    std::cout << "Reconstructed primary tracks: " << numberOfPrimaryTracks << std::endl
	      << "Reconstructed secondary tracks: " << numberOfSecondaryTracks << std::endl
	      << std::endl;
  }

  
  TFile *f = TFile::Open(outputFileName.c_str(), "recreate");
 if( !f || !f->IsWritable() ){
   std::cout << "Output file not opened correctly !!" << std::endl;
 }


  Histonumtoeventgen->Write();
  Histonumshowevent->Write();
  Histonumcleanevent->Write();
   
 DNDetaT2Clean->Write();
 DNDetaT2Show->Write();
 DNDetaT2ALL->Write();
 
 CutZimpLeftH0->Write();
 CutZimpRightH0->Write();

 CutZimpLeftH1->Write();
 CutZimpRightH1->Write();

 CutZimpLeftH2->Write();
 CutZimpRightH2->Write();

 CutZimpLeftH3->Write();
 CutZimpRightH3->Write();



 ArmPlus_TrkEtaEfficiencyCumulative->Write(); 
 ArmMinus_TrkEtaEfficiencyCumulative->Write(); 

 ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts->Write(); 
 ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts->Write(); 

 EtaResolFreeGenerator->Write(); 
 EtaResolFreeGeant->Write(); 
 EtaResolGenerator->Write(); 
 EtaResolGeant->Write();

 /*
 Pythia_UnstablePDG->Write(); Pythia_PDG_HavingStableDesc->Write();
 PID_ofSecRecoTrack->Write(); Pythia_PDG_Status2_Theta003->Write(); 
 Pythia_PDG_HavingStableDesc_Theta003->Write();
 PID_GeantTrk_atleast4HitInT2->Write();
 PID_GeantTrk_atleast4HitInT2_StableAtIP->Write();
 PID_GeantTrk_atleast4HitInT2_ConsideredAsPrimary->Write();

 PID_DircetMotherForAFakeRecoPrimaryTrk->Write();
 PID_OldestMotherForAFakeRecoPrimaryTrk->Write();
 PID_ofSecRecoTrackAsaFakePrimary->Write();
 
 
 Pythia_PDG_HavingStableDesc_Theta001->Write();
 SecondaryTrksZvsEta->Write(); 
 SecondaryTrksZvsPId->Write(); 
 SecondaryTrksEtaFromIonPump->Write();
 SimuVtxPositionForAFakeRecoPrimaryTrk->Write();
 
 PdgSecondaryInT2EtaAcceptance->Write();
 PrimaryTrksEta->Write(); AllTrksEta->Write();SecondaryTrksEta->Write();
 AllRecoVtxYZ->Write();
 AllRecoVtxZ->Write();

 FractionOfPrimaryTrksRecoAsSecondary->Write(); 
 FractionOfPrimaryTrksInSecVtx->Write();
 FractionOfSecondaryTrksRecoAsPrimary->Write();
 FractionOfSecondaryTrksInPrimaryVtx->Write();
 NumRecoPrimVtx->Write();

 RecoVtx_ZvsR->Write();
 SecondaryVtx_ZvsR_ContributingInT2->Write();

 PrimaryVtxPositionContributingInT2->Write();
 SecondaryVtxPositionContributingInT2->Write();
 SecondaryVtxPositionContributingInT2MultiCount->Write();
 ALLVertexPositionYZ->Write(); 
 SecondaryVertexPositionYZ->Write(); 
 PrimaryVertexPositionYZ->Write(); 
 
 SecondaryVtx_R_ContributingInT2->Write();
 SecondaryVtx_Z_ContributingInT2->Write();

 PrimaryTrksEtavsZatRmin->Write(); 
 ParticlePdgvsEinT2->Write(); 
 ParticlePdgvsENotinT2->Write();

 PrimaryRecoTrksEta->Write();

 RecoTrksEtaALL->Write();RecoTrksZ0ALL->Write();
 SecondaryRecoTrksEta->Write();
 SecondaryRecoTrksEtaPointingToVtx ->Write();
 SecondaryRecoTrksR0PointingToVtx ->Write();
 SecondaryRecoTrksZ0PointingToVtx->Write();

 SecondaryRecoTrksR0_2Arms->Write();
 PrimaryRecoTrksR0_2Arms->Write();

 SecondaryRecoTrksR02PointingToVtx->Write();
 PrimaryRecoTrksZ0PointingToVtx->Write();
 PrimaryRecoTrksrR02PointingToVtx->Write();
 PrimaryRecoTrksrR0PointingToVtx->Write();
 SecondaryRecoTrksR02_2Arms->Write();
 SecondaryRecoTrksZ0_2Arms->Write();
 PrimaryRecoTrksR02_2Arms->Write();
 PrimaryRecoTrksZ0_2Arms ->Write(); 
 PrimaryRecoTrksZ0_2ArmsEsmall5->Write(); 
 PrimaryRecoTrksZ0_2ArmsEbig5->Write(); 


 PrimaryRecoTrksZ0_2ArmsMixedEnergy->Write(); 
 H0_PrimaryRecoTrksZ0_2Arms->Write(); 
 H0_PrimaryRecoTrksZ0_2ArmsEbig5->Write();

 H0_PrimaryRecoTrksBX_2Arms->Write();
 H0_PrimaryRecoTrksBY_2Arms->Write();
 H0_SecondaryRecoTrksBX_2Arms->Write();
 H0_SecondaryRecoTrksBY_2Arms->Write();

 H0_PrimaryRecoTrksR0_2ArmsEbig5->Write();  
 H0_PrimaryRecoTrksEta_2ArmsEbig5R0Less20->Write();  
 H0_PrimaryRecoTrksEta_2ArmsEbig5->Write();  
 H0_PrimaryRecoTrksZ0_2ArmsEbig5R0Less20->Write(); 

 H0_RecoTrksEta->Write();
 H0_SecondaryRecoTrksEta->Write();
 H0_RecoTrksEtaR0Small20->Write();
 H0_SecondaryRecoTrksEtaR0Small20->Write();

 RecoTrkPrimaryMostProbEnergy->Write();
 RecoTrkPrimaryNumHitEnergy->Write();

 // H0_PrimaryZ0_BisMatch->Write();
 H0_PrimaryRecoTrksZ0_2Arms_BisMatch->Write();
 H0_PrimaryZ0_OrtogImpact_BisMatch->Write();  
 H0_PrimaryXY_OrtogImpact_BisMatch->Write();
 
 

 H0_PrimaryZ0_OrtogImpact->Write();
 H0_Z0_OrtogImpact->Write();  
 H0_RecoTrksZ0->Write();  
 H0_PrimaryZ0_OrtogImpactEbig5->Write();  

 H0_XY_OrtogImpact->Write();  
 H0_PrimaryXY_OrtogImpact->Write();  
 H0_PrimaryXY_OrtogImpactEbig5->Write();  
   
 PrimaryRecoTrksZ0_2Arms_EtaCut4_R0Cut20->Write(); 

 NumRecoSecondary_VsNumRecoPrimary->Write(); 
 NumRecoSecondary_VsNumRecoPrimary_EtaCut4_R0Cut20->Write(); 
 
 RecoTrksR0ALL->Write();  

 RecoTrksEtaALL_EtaCut4_R0Cut20->Write();  
 SecondaryRecoTrksEta_2Arms_EtaCut4_R0Cut20->Write();  
 PrimaryRecoTrksEta_2Arms_EtaCut4_R0Cut20->Write(); 

 Reco_Trk_DeltaThetaX_vs_DeltaThetaY_FromSecondary->Write();   
 Reco_Trk_DeltaThetaX_vs_DeltaThetaY_FromPrimary->Write();   

 Reco_TrkZ0_vs_R0_AllH0_FromPrimary->Write();   
 Reco_TrkZ0_vs_R0_AllH0_FromSecondary->Write();   
 Reco_TrkZ0_vs_Eta_AllH0_FromPrimary->Write();   
 Reco_TrkZ0_vs_Eta_AllH0_FromSecondary->Write();  
 Reco_TrkZ0_vs_Rmin_AllH0_FromSecondary->Write(); 
 Reco_TrkZ0_vs_Rmin_AllH0_FromPrimary->Write(); 

 RecoTrksR0ALL_EtaCut4_R0Cut20->Write();  
 SecondaryRecoTrksR0_2Arms_EtaCut4_R0Cut20->Write();  
 PrimaryRecoTrksR0_2Arms_EtaCut4_R0Cut20->Write(); 


 RecoTrksZ0ALL_R0Cut10->Write();
 RecoTrksZ0ALL_R0Cut20->Write();
 SecondaryRecoTrksZ0_2Arms_R0Cut20->Write();
 SecondarySmallZ0_EtavsR0->Write(); 
 RecoTrksZ0ALL_EtaCut4->Write(); 

 SecondaryRecoTrksZ0_2Arms_EtaCut4_R0Cut20->Write(); 
 RecoTrksZ0ALL_EtaCut4_R0Cut20->Write(); 

 SecondaryBigZ0_EtavsR0->Write(); 

 //ProbabilityToRecoFakePrimary_vs_Eta->Write(); 
 //ProbabilityToRecoRealPrimary_vs_Eta->Write(); 
 
 RecoDnDeta2FromVtx->Write();  RecoDnDeta2NoCut->Write(); 
 PrimaryRecoDnDeta2FromVtx->Write();  PrimaryRecoDnDeta2NoCut->Write(); 
 FakeSecondaryRecoDnDeta2FromVtx->Write();
 RecoDnDeta2UnfoldCut2->Write();
 RecoDnDeta2UnfoldCut2_Primary->Write();
 RecoDnDeta2UnfoldCut->Write();
 RecoDnDeta2UnfoldCut_Primary->Write();
 

  
  H0_sameTrk_PadCol->Write(); H0_sameTrk_PadRow->Write(); H0_sameTrk_StripNumb->Write(); H0_sameTrk_StripColMaxDist->Write(); 
  H0_sameTrk_PadRowMaxDist->Write();  H0_sameTrk_PadColMaxDist->Write();
  H0_sameTrk_StripRowMaxDist->Write();H0_sameTrk_PadColMaxDist->Write(); 
  
  CumulativeGeantTrkHadronR->Write();
  CumulativeGeantTrkElectronR->Write(); 
  CumulativeGeantTrkElectronXY->Write(); 
  DigiPadOccupancy->Write(); 

  GeantThetaXofPrimary->Write(); 
  GeantThetaYofPrimary->Write(); 

  TrkRecoEfficiencyProfilePlus->Write(); 
  TrkRecoEfficiencyProfileMinus->Write();
  GeantTrkMultiplicityfrequencyPlus->Write();
  SelectedEventGeantPad->Write(); SelectedEventGeantExitPad->Write();  SelectedEventGeantEntryPad->Write();
  SelectedEventGeantPadTrkWithMoreThan3Hit->Write(); SelectedEventGeantPadTrkWithMoreThan3HitEntry_Weight->Write();
  SelectedEventGeantPadTrkWithMoreThan3HitExit_Weight->Write();SelectedEventGeantPadTrkWithMoreThan3Tracing->Write();
  SelectedEventGeantPadTrkWithMoreThan3TracingMultiCount->Write();
  SelectedEventGeantPadTrkWithMoreThan3TracingGEO->Write();SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY->Write();
  SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY_Noelectron->Write();SelectedEventGeantPadTrkWithMoreThan3_PythiaChinIP->Write();
  SelectedEventGeantPadTrkWithMoreThan3_PythiaChinIPGEOXY->Write(); 
 
  Purity_andEfficiency->Write(); 
  Purity_andEfficiencyNoGhost->Write(); 
  RecoTrkEta_MatchingPrimary->Write(); 

 
  
  ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts_2->Write(); 
  ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts_2->Write();   
  
  ArmPlus_TrkEtaEfficiencyCumulativeSecondary->Write();   
  ArmMinus_TrkEtaEfficiencyCumulativeSecondary->Write();   
  ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCutsSecondary->Write();   
  ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCutsSecondary->Write();   
  ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts_2Secondary->Write();   
  ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts_2Secondary->Write();  


  
  H0_TrkEtaEfficiencyCumulative->Write();  
  H1_TrkEtaEfficiencyCumulative->Write();  
  H2_TrkEtaEfficiencyCumulative->Write();  
  H3_TrkEtaEfficiencyCumulative->Write();  

  H0_TrkEtaEfficiencyCumulative_UnfoldCuts->Write();  
  H1_TrkEtaEfficiencyCumulative_UnfoldCuts->Write();  
  H2_TrkEtaEfficiencyCumulative_UnfoldCuts->Write();  
  H3_TrkEtaEfficiencyCumulative_UnfoldCuts->Write();  

  H0_TrkEtaEfficiencyCumulative_UnfoldCuts_2->Write();  
  H1_TrkEtaEfficiencyCumulative_UnfoldCuts_2->Write();  
  H2_TrkEtaEfficiencyCumulative_UnfoldCuts_2->Write();  
  H3_TrkEtaEfficiencyCumulative_UnfoldCuts_2->Write();  
  
  
  H0TrkEffi_AxisAvgMult_ZImpactCut->Write();  
  H1TrkEffi_AxisAvgMult_ZImpactCut->Write();  
  H2TrkEffi_AxisAvgMult_ZImpactCut->Write();  
  H3TrkEffi_AxisAvgMult_ZImpactCut->Write();  
  H0TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut->Write(); 
  H1TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut->Write(); 
  H2TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut->Write(); 
  H3TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut->Write();

  H0TrkEffi_AxisAvgMult_NORecCut_GeantZCut->Write();
  H1TrkEffi_AxisAvgMult_NORecCut_GeantZCut->Write();
  H2TrkEffi_AxisAvgMult_NORecCut_GeantZCut->Write();
  H3TrkEffi_AxisAvgMult_NORecCut_GeantZCut->Write();
 
  H0TrkEffi_AxisAvgMult_NORecCut_GeantZCut_MatchProj->Write();
  H0TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut_MatchProj->Write();
  
  
  NumH0PlaneActiveWhenGeantEffiFail->Write();
  minDxDyWhenGeantEffiFail->Write();
  Numtrk_vs_NumtrkpasCutsWhenEffiFail->Write();
  ZImpOfARecTrkWhenEffiFail->Write();
  //minDphiWhenGeantEffiFail->Write();
  //minDrWhenGeantEffiFail->Write();
  NumH0PlaneActive_vsMultipl_WhenGeantEffiFail->Write();
  NumH0TrkWhenGeantEffiFail->Write();
  NumHitInRecTrk->Write();
  Cumulative_GeantTrk_ZImpactRefQ->Write(); 
  Geant_Prim_Or_Sec_Unfold->Write();   
  Cumulative_GeantTrk_ZImpactRefQ_PtAboveZ5->Write(); 
  Cumulative_GeantTrk_ZImpactRefQ_E->Write(); 
  Cumulative_GeantTrk_ZImpactRefQ_EAboveZ5->Write(); 
  PSimHitT2_EnergyLoss->Write();

 */


   
 Histonumtoeventgen->Write();
  Histonumshowevent->Write();  




  RPhiGeantHitinH0Col32->Write();
  XYGeantHitinH0Col32->Write();

 PrimaryGeantDNDeta2_IfOneTrkInT2_H0->Write(); 
 PrimaryGeantDNDeta2_IfOneTrkInT2_H1->Write(); corposprimaryXYQ0->Write();corposprimaryXYQ1->Write();
PrimaryGeantDNDeta2_IfOneTrkInT2_H1sometrkexcluding64->Write();
 ZImpOfARecTrkWhenEffiFail->Write();
 Numtrk_vs_NumtrkpasCutsWhenEffiFail->Write(); minDxDyWhenGeantEffiFail->Write(); minDxDyWhenGeantEffiFailH1->Write();
 NumH0TrkWhenGeantEffiFail->Write();
   NumH0PlaneActiveWhenGeantEffiFail->Write(); NumH0PlaneActive_vsMultipl_WhenGeantEffiFail->Write();

 H0H1_PhiGeantTrkWhenEffiFail->Write();
  H0H1_EtaPhiGeantTrkWhenEffiFail->Write();

ZImpOfAGeantTrkWhenEffiFailH1->Write();
ZImpOfAGeantTrkWhenEffiFailH0->Write();

  //H0_EtaPhiGeantTrkWhenEffiFail->Write();
 H1PrimaryTrkPdgID->Write(); H0PrimaryTrkPdgID->Write(); H0PrimaryTrkE->Write(); H1PrimaryTrkE->Write();
  
 /*
 gDirectory = f->mkdir("RotturaDiPalle");
 for(unsigned int m=0;m<40; m++){
   PlaneGeantPadColXY[m]->Write();
   PlaneGeantToRecXY[m]->Write();
 }
 */

 PrimaryGeantDNDeta2->Write();  PrimaryGeanteta2->Write(); 
 PrimaryGeantDNDeta2_IfOneRecoZImpInH0->Write(); 
 PrimaryGeantDNDeta2H0_ZCut->Write();
 DNDetaMBALLMCGenerator_GeneratorTriggered->Write(); 
 DNDetaMBALLMCGenerator_GeneratorTriggeredAtLeastOne5364->Write(); 
  gDirectory = f->mkdir("BinnedPrimaryEfficiency");

  H0_TrkPhiEfficiencyCumulative_UnfoldCutsCumul->Write();
  H1_TrkPhiEfficiencyCumulative_UnfoldCutsCumul->Write();
 

  PlaneClusterActiveCounter->Write();PrimTrkTestedQuarter->Write();
  PrimTrkTestedQuarterUnamb->Write();PrimTrkTestedPlaneQuarter->Write();

 
  for(unsigned int m=0;m<10; m++)
    {
      H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]->Write();
      //H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_56Division[m]->Write();
      //    H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_ZcutOnGeant[m]->Write();
      H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]->Write();
      H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]->Write();H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_PhiCut[m]->Write();
      H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[m]->Write();

      H0_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult[m]->Write();
      H1_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult[m]->Write();
    }
  
  for(unsigned int m=0;m<=22; m++){
    Associated_PrimaryTrkZImpact_H0[m]->Write();
    Associated_PrimaryTrkZImpact_H1[m]->Write();
    Associated_PrimaryTrkZImpact_H2[m]->Write();
    Associated_PrimaryTrkZImpact_H3[m]->Write();
  }

  /*
  gDirectory = f->mkdir("Unfolding");

  
  for(unsigned int m=0;m<50;m++)
    {ArmPlus_TrkEtaEfficiency_FixedGeantMult[m]->Write();}

  for(unsigned int m=0;m<50;m++)
    {
      ArmPlus_TrkRecoCutMult_FixedGeantMult[m]->Write();
      ArmPlus_TrkRecoCutMult_FixedGeantMult2[m]->Write();
    }
  
  */

  /*
  gDirectory = f->mkdir("K0");

  for(unsigned int m=0;m<50;m++){SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadColvsPadRow[m]->Write();}
  for(unsigned int m=0;m<50;m++){SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOXY[m]->Write();}     
  for(unsigned int m=0;m<50;m++){SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEORPhi[m]->Write();}     
  for(unsigned int m=0;m<50;m++){SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOYZ[m]->Write();}
  for(unsigned int m=0;m<50;m++){SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOXZ[m]->Write();}  
  for(unsigned int m=0;m<50;m++){SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadColvsZ[m]->Write();}  
  for(unsigned int m=0;m<50;m++){SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadRowvsZ[m]->Write();}



  gDirectory = f->mkdir("ChParticleIP");

  for(unsigned int m=0;m<50;m++){SelectedEventGeantPadTrkWithMoreThan3TracingGEOXZ_NoelectronPrimary[m]->Write();}
  for(unsigned int m=0;m<50;m++){SelectedEventGeantPadTrkWithMoreThan3TracingGEOYZ_NoelectronPrimary[m]->Write();}
  for(unsigned int m=0;m<50;m++){SelectedEventGeantPadTrkWithMoreThan3TracingPadColvsZ_NoelectronPrimary[m]->Write();}
  for(unsigned int m=0;m<50;m++){SelectedEventGeantPadTrkWithMoreThan3TracingPadRowvsZ_NoelectronPrimary[m]->Write();} 
  for(unsigned int m=0;m<50;m++){SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY_NoelectronPrimary[m]->Write();}



  gDirectory = f->mkdir("DIGI");
  
  for(unsigned int m=0;m<50;m++){SelectedEventDIGIPad_GTrkWithMoreThan3_GEOYZ[m]->Write(); }
  for(unsigned int m=0;m<50;m++){SelectedEventDIGIPad_GTrkWithMoreThan3_GEOXZ[m]->Write(); }
  for(unsigned int m=0;m<50;m++){SelectedEventDIGIPad_GTrkWithMoreThan3_GEOXY[m]->Write(); }
  for(unsigned int m=0;m<50;m++){SelectedEventDIGIPad_GTrkWithMoreThan3_PadColvsPadRow[m]->Write();}
  for(unsigned int m=0;m<50;m++){SelectedEventDIGIPad_GTrkWithMoreThan3_PadColvsZ[m]->Write(); }
  for(unsigned int m=0;m<50;m++){SelectedEventDIGIPad_GTrkWithMoreThan3_PadRowvsZ[m]->Write(); }


  gDirectory = f->mkdir("NEW_ROADFINDER");
  for(unsigned int m=0;m<50;m++){SelectedEventRoadPadFinderClu_PadRowvsZ[m]->Write(); }
  for(unsigned int m=0;m<50;m++){SelectedEventRoadPadFinderClu_PadColvsZ[m]->Write(); }
  for(unsigned int m=0;m<50;m++){SelectedEventRoadPadFinderClu_GEOXZ[m]->Write(); }
  for(unsigned int m=0;m<50;m++){SelectedEventRoadPadFinderClu_GEOYZ[m]->Write(); }

   
  gDirectory = f->mkdir("NEW_ROADFINDER_FullRoad");
  
  for(unsigned int m=0;m<50;m++){SelectedEventRoadPadFinder_FullR_Clu_PadRowvsZ[m]->Write(); }
  for(unsigned int m=0;m<50;m++){SelectedEventRoadPadFinder_FullR_Clu_PadColvsZ[m]->Write(); }
  for(unsigned int m=0;m<50;m++){SelectedEventRoadPadFinder_FullR_Clu_GEOXZ[m]->Write(); }
  for(unsigned int m=0;m<50;m++){SelectedEventRoadPadFinder_FullR_Clu_GEOYZ[m]->Write(); }
  
  */


  NumEventTrackTrigger->Write();
  NumEventGeneratorTrigger->Write();
    eta2geant_forEffHisto3->Write(); eta2geant_forEffHisto1->Write(); eta2geant_forEffHisto2->Write();

   H0_PrimaryRecoTrksZ0_Vs_Z0OrtogImpact_BisMatch->Write();
   H0_SecondaryRecoTrksZ0_Vs_Z0OrtogImpact_BisMatch->Write();
   H0_PrimaryRecoTrksZ0Ortog_Vs_R0Ortog_BisMatch->Write();
   H0_SecondaryRecoTrksZ0Ortog_Vs_R0Ortog_BisMatch->Write();
   H0_PrimaryRecoTrksX0Ortog_Vs_Y0Ortog_BisMatch->Write();
   H0_SecondaryRecoTrksX0Ortog_Vs_Y0Ortog_BisMatch->Write(); 
   H0_PrimaryRecoTrksZ0Ortog_Vs_X0Ortog_BisMatch->Write();
   H0_SecondaryRecoTrksZ0Ortog_Vs_X0Ortog_BisMatch->Write(); 
   H0_PrimaryRecoTrksZ0Ortog_Vs_Y0Ortog_BisMatch->Write();
   H0_SecondaryRecoTrksZ0Ortog_Vs_Y0Ortog_BisMatch->Write();

    CumulativeZimpH1->Write(); CumulativeZimpH0->Write();    CumulativeZimpH2->Write(); CumulativeZimpH3->Write();
   eta2geant_forEffQ0->Write();
   eta2geant_forEffQ1->Write();
   eta2geant_forEffQ2->Write();
   eta2geant_forEffQ3->Write();
   eta2geant_forEffQ0eff->Write();
   eta2geant_forEffQ1eff->Write();
   eta2geant_forEffQ2eff->Write();
   eta2geant_forEffQ3eff->Write();
   H0_BinMultStatFail->Write();H1_BinMultStatFail->Write();H2_BinMultStatFail->Write();H3_BinMultStatFail->Write();
   H0_BinMultStat->Write();H1_BinMultStat->Write();H2_BinMultStat->Write();H3_BinMultStat->Write();
    NumHitInRecTrkH0->Write();   NumHitInRecTrkH1->Write();   NumHitInRecTrkH2->Write();   NumHitInRecTrkH3->Write();

    AveragePadCLSDistrH1->Write();    AveragePadCLSDistrH0->Write();   
    rechitUnitID->Write();
    simhitUnitID->Write();
    Q1PadRCPrimary->Write();
    
    Q0PadRCPrimary->Write();
    

    H0_TrkPhiEfficiencyCumulative_UnfoldCutsCumul->Write();
    H1_TrkPhiEfficiencyCumulative_UnfoldCutsCumul->Write();

     PadRowGeantRecoH0->Write();
     PadRowGeantRecoH0_bis->Write();
     PadRowGeantRecoH1->Write();
     PadRowGeantRecoH1_bis->Write();
 
     ZImpOfAPrimaryGeantTrk_Cumulative->Write();
     ZImpOfAPrimaryGeantTrk_Cumulative_match->Write();
     ZImpOfAPrimaryGeantTrk_CumulativeZMinCut->Write();
     ZImpOfAPrimaryGeantTrk_CumulativeZMinCut_match->Write();

     NumChPartInT2PerEvt->Write(); 
     NumTrkInT2PerEvt7mcut->Write();
     NumTrkInT2PerEvt->Write();

     f->SaveSelf(kTRUE);
     
     f->Close();

   std::cout<<"Job end. Last call"<<std::endl;

}







DEFINE_FWK_MODULE(T2BackgroundAn);
