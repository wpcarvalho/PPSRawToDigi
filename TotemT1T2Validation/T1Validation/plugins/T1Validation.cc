#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TotemT1T2Validation/T1Validation/interface/T1Validation.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "HepMC/GenEvent.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "TFile.h"
#include "TH1.h"

//#define _PRINT_
#define _ALSO_DIGIS_
#define _ALSO_RECO_
#define _PRINT_HITS_

T1Validation::T1Validation(const edm::ParameterSet& iConfig):_SIM_(1),_DIGI_(1),_RECO_(1)
{
  t1ClusterCollectionLabel = iConfig.getParameter<edm::InputTag>("T1ClusterCollectionLabel");
  t1RoadCollectionLabel = iConfig.getParameter<edm::InputTag>("T1RoadCollectionLabel");
  t1DigiWireCollectionLabel = iConfig.getParameter<edm::InputTag>("T1DigiWireCollectionLabel");
  t1DigiVfatCollectionLabel = iConfig.getParameter<edm::InputTag>("T1DigiVfatCollectionLabel");
//now do what ever initialization is needed
  std::string thegeometryfile = "Geometry/TotemGeometry/data/T1_data_geometry.dat";
  layer = std::auto_ptr<T1Geometry>(new T1Geometry(thegeometryfile));
  _SIM_ = iConfig.getParameter<double>("SIM");
  _DIGI_ = iConfig.getParameter<double>("DIGI");
  _RECO_ = iConfig.getParameter<double>("RECO");
  trackLabel= iConfig.getParameter<std::string>("TrackLabel");
  outputFileName = iConfig.getParameter<std::string>("OutputFile");

  MaxTracks = 10000;
  if(iConfig.exists("MaxTracks")){
    MaxTracks = iConfig.getParameter<int> ("MaxTracks");
  }
  MinTracks = 0;
  if(iConfig.exists("MinTracks")){
    MinTracks = iConfig.getParameter<int> ("MinTracks");
  }



  copy_of_hits = new std::vector<PSimHit>();

}


T1Validation::~T1Validation()
{
// do anything here that needs to be done at destruction time
// (e.g. close files, deallocate resources etc.)
  delete copy_of_hits;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
T1Validation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;


  int _NWireCh0Arm0Pl0 = 0;
  int _NWireCh0Arm0Pl1 = 0;
  int _NWireCh0Arm0Pl2 = 0;
  int _NWireCh0Arm0Pl3 = 0;
  int _NWireCh0Arm0Pl4 = 0;

  int _NWireCh0Arm1Pl0 = 0;
  int _NWireCh0Arm1Pl1 = 0;
  int _NWireCh0Arm1Pl2 = 0;
  int _NWireCh0Arm1Pl3 = 0;
  int _NWireCh0Arm1Pl4 = 0;

  int _NStripACh0Arm0Pl0 = 0;
  int _NStripACh0Arm0Pl1 = 0;
  int _NStripACh0Arm0Pl2 = 0;
  int _NStripACh0Arm0Pl3 = 0;
  int _NStripACh0Arm0Pl4 = 0;

  int _NStripACh0Arm1Pl0 = 0;
  int _NStripACh0Arm1Pl1 = 0;
  int _NStripACh0Arm1Pl2 = 0;
  int _NStripACh0Arm1Pl3 = 0;
  int _NStripACh0Arm1Pl4 = 0;

  int _NStripBCh0Arm0Pl0 = 0;
  int _NStripBCh0Arm0Pl1 = 0;
  int _NStripBCh0Arm0Pl2 = 0;
  int _NStripBCh0Arm0Pl3 = 0;
  int _NStripBCh0Arm0Pl4 = 0;

  int _NStripBCh0Arm1Pl0 = 0;
  int _NStripBCh0Arm1Pl1 = 0;
  int _NStripBCh0Arm1Pl2 = 0;
  int _NStripBCh0Arm1Pl3 = 0;
  int _NStripBCh0Arm1Pl4 = 0;


  int _NWireCh1Arm0Pl0 = 0;
  int _NWireCh1Arm0Pl1 = 0;
  int _NWireCh1Arm0Pl2 = 0;
  int _NWireCh1Arm0Pl3 = 0;
  int _NWireCh1Arm0Pl4 = 0;

  int _NWireCh1Arm1Pl0 = 0;
  int _NWireCh1Arm1Pl1 = 0;
  int _NWireCh1Arm1Pl2 = 0;
  int _NWireCh1Arm1Pl3 = 0;
  int _NWireCh1Arm1Pl4 = 0;

  int _NStripACh1Arm0Pl0 = 0;
  int _NStripACh1Arm0Pl1 = 0;
  int _NStripACh1Arm0Pl2 = 0;
  int _NStripACh1Arm0Pl3 = 0;
  int _NStripACh1Arm0Pl4 = 0;

  int _NStripACh1Arm1Pl0 = 0;
  int _NStripACh1Arm1Pl1 = 0;
  int _NStripACh1Arm1Pl2 = 0;
  int _NStripACh1Arm1Pl3 = 0;
  int _NStripACh1Arm1Pl4 = 0;

  int _NStripBCh1Arm0Pl0 = 0;
  int _NStripBCh1Arm0Pl1 = 0;
  int _NStripBCh1Arm0Pl2 = 0;
  int _NStripBCh1Arm0Pl3 = 0;
  int _NStripBCh1Arm0Pl4 = 0;

  int _NStripBCh1Arm1Pl0 = 0;
  int _NStripBCh1Arm1Pl1 = 0;
  int _NStripBCh1Arm1Pl2 = 0;
  int _NStripBCh1Arm1Pl3 = 0;
  int _NStripBCh1Arm1Pl4 = 0;


  int _NWireCh2Arm0Pl0 = 0;
  int _NWireCh2Arm0Pl1 = 0;
  int _NWireCh2Arm0Pl2 = 0;
  int _NWireCh2Arm0Pl3 = 0;
  int _NWireCh2Arm0Pl4 = 0;

  int _NWireCh2Arm1Pl0 = 0;
  int _NWireCh2Arm1Pl1 = 0;
  int _NWireCh2Arm1Pl2 = 0;
  int _NWireCh2Arm1Pl3 = 0;
  int _NWireCh2Arm1Pl4 = 0;

  int _NStripACh2Arm0Pl0 = 0;
  int _NStripACh2Arm0Pl1 = 0;
  int _NStripACh2Arm0Pl2 = 0;
  int _NStripACh2Arm0Pl3 = 0;
  int _NStripACh2Arm0Pl4 = 0;

  int _NStripACh2Arm1Pl0 = 0;
  int _NStripACh2Arm1Pl1 = 0;
  int _NStripACh2Arm1Pl2 = 0;
  int _NStripACh2Arm1Pl3 = 0;
  int _NStripACh2Arm1Pl4 = 0;

  int _NStripBCh2Arm0Pl0 = 0;
  int _NStripBCh2Arm0Pl1 = 0;
  int _NStripBCh2Arm0Pl2 = 0;
  int _NStripBCh2Arm0Pl3 = 0;
  int _NStripBCh2Arm0Pl4 = 0;

  int _NStripBCh2Arm1Pl0 = 0;
  int _NStripBCh2Arm1Pl1 = 0;
  int _NStripBCh2Arm1Pl2 = 0;
  int _NStripBCh2Arm1Pl3 = 0;
  int _NStripBCh2Arm1Pl4 = 0;


  int _NWireCh3Arm0Pl0 = 0;
  int _NWireCh3Arm0Pl1 = 0;
  int _NWireCh3Arm0Pl2 = 0;
  int _NWireCh3Arm0Pl3 = 0;
  int _NWireCh3Arm0Pl4 = 0;

  int _NWireCh3Arm1Pl0 = 0;
  int _NWireCh3Arm1Pl1 = 0;
  int _NWireCh3Arm1Pl2 = 0;
  int _NWireCh3Arm1Pl3 = 0;
  int _NWireCh3Arm1Pl4 = 0;

  int _NStripACh3Arm0Pl0 = 0;
  int _NStripACh3Arm0Pl1 = 0;
  int _NStripACh3Arm0Pl2 = 0;
  int _NStripACh3Arm0Pl3 = 0;
  int _NStripACh3Arm0Pl4 = 0;

  int _NStripACh3Arm1Pl0 = 0;
  int _NStripACh3Arm1Pl1 = 0;
  int _NStripACh3Arm1Pl2 = 0;
  int _NStripACh3Arm1Pl3 = 0;
  int _NStripACh3Arm1Pl4 = 0;

  int _NStripBCh3Arm0Pl0 = 0;
  int _NStripBCh3Arm0Pl1 = 0;
  int _NStripBCh3Arm0Pl2 = 0;
  int _NStripBCh3Arm0Pl3 = 0;
  int _NStripBCh3Arm0Pl4 = 0;

  int _NStripBCh3Arm1Pl0 = 0;
  int _NStripBCh3Arm1Pl1 = 0;
  int _NStripBCh3Arm1Pl2 = 0;
  int _NStripBCh3Arm1Pl3 = 0;
  int _NStripBCh3Arm1Pl4 = 0;


  int _NWireCh4Arm0Pl0 = 0;
  int _NWireCh4Arm0Pl1 = 0;
  int _NWireCh4Arm0Pl2 = 0;
  int _NWireCh4Arm0Pl3 = 0;
  int _NWireCh4Arm0Pl4 = 0;

  int _NWireCh4Arm1Pl0 = 0;
  int _NWireCh4Arm1Pl1 = 0;
  int _NWireCh4Arm1Pl2 = 0;
  int _NWireCh4Arm1Pl3 = 0;
  int _NWireCh4Arm1Pl4 = 0;

  int _NStripACh4Arm0Pl0 = 0;
  int _NStripACh4Arm0Pl1 = 0;
  int _NStripACh4Arm0Pl2 = 0;
  int _NStripACh4Arm0Pl3 = 0;
  int _NStripACh4Arm0Pl4 = 0;

  int _NStripACh4Arm1Pl0 = 0;
  int _NStripACh4Arm1Pl1 = 0;
  int _NStripACh4Arm1Pl2 = 0;
  int _NStripACh4Arm1Pl3 = 0;
  int _NStripACh4Arm1Pl4 = 0;

  int _NStripBCh4Arm0Pl0 = 0;
  int _NStripBCh4Arm0Pl1 = 0;
  int _NStripBCh4Arm0Pl2 = 0;
  int _NStripBCh4Arm0Pl3 = 0;
  int _NStripBCh4Arm0Pl4 = 0;

  int _NStripBCh4Arm1Pl0 = 0;
  int _NStripBCh4Arm1Pl1 = 0;
  int _NStripBCh4Arm1Pl2 = 0;
  int _NStripBCh4Arm1Pl3 = 0;
  int _NStripBCh4Arm1Pl4 = 0;


  int _NWireCh5Arm0Pl0 = 0;
  int _NWireCh5Arm0Pl1 = 0;
  int _NWireCh5Arm0Pl2 = 0;
  int _NWireCh5Arm0Pl3 = 0;
  int _NWireCh5Arm0Pl4 = 0;

  int _NWireCh5Arm1Pl0 = 0;
  int _NWireCh5Arm1Pl1 = 0;
  int _NWireCh5Arm1Pl2 = 0;
  int _NWireCh5Arm1Pl3 = 0;
  int _NWireCh5Arm1Pl4 = 0;

  int _NStripACh5Arm0Pl0 = 0;
  int _NStripACh5Arm0Pl1 = 0;
  int _NStripACh5Arm0Pl2 = 0;
  int _NStripACh5Arm0Pl3 = 0;
  int _NStripACh5Arm0Pl4 = 0;

  int _NStripACh5Arm1Pl0 = 0;
  int _NStripACh5Arm1Pl1 = 0;
  int _NStripACh5Arm1Pl2 = 0;
  int _NStripACh5Arm1Pl3 = 0;
  int _NStripACh5Arm1Pl4 = 0;

  int _NStripBCh5Arm0Pl0 = 0;
  int _NStripBCh5Arm0Pl1 = 0;
  int _NStripBCh5Arm0Pl2 = 0;
  int _NStripBCh5Arm0Pl3 = 0;
  int _NStripBCh5Arm0Pl4 = 0;

  int _NStripBCh5Arm1Pl0 = 0;
  int _NStripBCh5Arm1Pl1 = 0;
  int _NStripBCh5Arm1Pl2 = 0;
  int _NStripBCh5Arm1Pl3 = 0;
  int _NStripBCh5Arm1Pl4 = 0;


  bool bool_MinMaxTracks=true;


  if(_RECO_>0){

#ifdef _ALSO_RECO_

    edm::Handle<T1T2TrackCollection> trackCollection;
    iEvent.getByLabel(trackLabel,"T1TrackColl",trackCollection);

    T1T2TrackCollection::const_iterator TC_it;


    if(trackCollection->size() > (unsigned int)MaxTracks || trackCollection->size() < (unsigned int)MinTracks) bool_MinMaxTracks=false;

    if(bool_MinMaxTracks){

      for(TC_it=trackCollection->begin(); TC_it!=trackCollection->end(); TC_it++){
	hAllEtaRec -> Fill((*TC_it).Eta() );
	hAllPhiRec -> Fill((*TC_it).Phi() );
	hChiSquaredOverN->Fill((*TC_it).ChiSquaredOverN() );
	hChiSquaredXOverN->Fill((*TC_it).ChiSquaredXOverN() );

	hChiSquaredYOverN->Fill((*TC_it).ChiSquaredYOverN() );
	hChiSquaredProb->Fill( TMath::Prob( (*TC_it).ChiSquared(), 2*( (*TC_it).GetHitEntries()-2 ) ));

	hChiSquaredXProb->Fill( TMath::Prob( (*TC_it).ChiSquaredX(), ( (*TC_it).GetHitEntries()-2 ) ));

	hChiSquaredYProb->Fill( TMath::Prob( (*TC_it).ChiSquaredY(), ( (*TC_it).GetHitEntries()-2 ) ));

      }




      edm::Handle<T1RecHit2DCollection>myRecoColl;
      iEvent.getByLabel("t1rechit",myRecoColl);

 

      T1RecHit2DCollection::const_iterator T1Reco_it;

      unsigned int NHITS[2][5][6] = { { { 0 } } };
  
      for(T1Reco_it=myRecoColl->begin();T1Reco_it!=myRecoColl->end();T1Reco_it++){

	T1DetId t1Id((*T1Reco_it).t1DetId());

	int ARM = t1Id.Arm();
	int PLANE = t1Id.Plane();
	int CSC = t1Id.CSC();
	float LX = (*T1Reco_it).localPosition().x();
	float LY = (*T1Reco_it).localPosition().y();


	switch(ARM){
	case 0:
	  switch(PLANE){
	  case 0:
	    switch(CSC){
	    case 0:
	      RecoLX_Arm0Plane0CSC0->Fill(LX);
	      RecoLY_Arm0Plane0CSC0->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 1:
	      RecoLX_Arm0Plane0CSC1->Fill(LX);
	      RecoLY_Arm0Plane0CSC1->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 2:
	      RecoLX_Arm0Plane0CSC2->Fill(LX);
	      RecoLY_Arm0Plane0CSC2->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 3:
	      RecoLX_Arm0Plane0CSC3->Fill(LX);
	      RecoLY_Arm0Plane0CSC3->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 4:
	      RecoLX_Arm0Plane0CSC4->Fill(LX);
	      RecoLY_Arm0Plane0CSC4->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 5:
	      RecoLX_Arm0Plane0CSC5->Fill(LX);
	      RecoLY_Arm0Plane0CSC5->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    }
	    break;
	  case 1:
	    switch(CSC){
	    case 0:
	      RecoLX_Arm0Plane1CSC0->Fill(LX);
	      RecoLY_Arm0Plane1CSC0->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 1:
	      RecoLX_Arm0Plane1CSC1->Fill(LX);
	      RecoLY_Arm0Plane1CSC1->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 2:
	      RecoLX_Arm0Plane1CSC2->Fill(LX);
	      RecoLY_Arm0Plane1CSC2->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 3:
	      RecoLX_Arm0Plane1CSC3->Fill(LX);
	      RecoLY_Arm0Plane1CSC3->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 4:
	      RecoLX_Arm0Plane1CSC4->Fill(LX);
	      RecoLY_Arm0Plane1CSC4->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 5:
	      RecoLX_Arm0Plane1CSC5->Fill(LX);
	      RecoLY_Arm0Plane1CSC5->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    }
	    break;
	  case 2:
	    switch(CSC){
	    case 0:
	      RecoLX_Arm0Plane2CSC0->Fill(LX);
	      RecoLY_Arm0Plane2CSC0->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 1:
	      RecoLX_Arm0Plane2CSC1->Fill(LX);
	      RecoLY_Arm0Plane2CSC1->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 2:
	      RecoLX_Arm0Plane2CSC2->Fill(LX);
	      RecoLY_Arm0Plane2CSC2->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 3:
	      RecoLX_Arm0Plane2CSC3->Fill(LX);
	      RecoLY_Arm0Plane2CSC3->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 4:
	      RecoLX_Arm0Plane2CSC4->Fill(LX);
	      RecoLY_Arm0Plane2CSC4->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 5:
	      RecoLX_Arm0Plane2CSC5->Fill(LX);
	      RecoLY_Arm0Plane2CSC5->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    }
	    break;
	  case 3:
	    switch(CSC){
	    case 0:
	      RecoLX_Arm0Plane3CSC0->Fill(LX);
	      RecoLY_Arm0Plane3CSC0->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 1:
	      RecoLX_Arm0Plane3CSC1->Fill(LX);
	      RecoLY_Arm0Plane3CSC1->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 2:
	      RecoLX_Arm0Plane3CSC2->Fill(LX);
	      RecoLY_Arm0Plane3CSC2->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 3:
	      RecoLX_Arm0Plane3CSC3->Fill(LX);
	      RecoLY_Arm0Plane3CSC3->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 4:
	      RecoLX_Arm0Plane3CSC4->Fill(LX);
	      RecoLY_Arm0Plane3CSC4->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 5:
	      RecoLX_Arm0Plane3CSC5->Fill(LX);
	      RecoLY_Arm0Plane3CSC5->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    }
	    break;
	  case 4:
	    switch(CSC){
	    case 0:
	      RecoLX_Arm0Plane4CSC0->Fill(LX);
	      RecoLY_Arm0Plane4CSC0->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 1:
	      RecoLX_Arm0Plane4CSC1->Fill(LX);
	      RecoLY_Arm0Plane4CSC1->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 2:
	      RecoLX_Arm0Plane4CSC2->Fill(LX);
	      RecoLY_Arm0Plane4CSC2->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 3:
	      RecoLX_Arm0Plane4CSC3->Fill(LX);
	      RecoLY_Arm0Plane4CSC3->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 4:
	      RecoLX_Arm0Plane4CSC4->Fill(LX);
	      RecoLY_Arm0Plane4CSC4->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 5:
	      RecoLX_Arm0Plane4CSC5->Fill(LX);
	      RecoLY_Arm0Plane4CSC5->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    }
	    break;
	  }
	  break;
	case 1:
	  switch(PLANE){
	  case 0:
	    switch(CSC){
	    case 0:
	      RecoLX_Arm1Plane0CSC0->Fill(LX);
	      RecoLY_Arm1Plane0CSC0->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 1:
	      RecoLX_Arm1Plane0CSC1->Fill(LX);
	      RecoLY_Arm1Plane0CSC1->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 2:
	      RecoLX_Arm1Plane0CSC2->Fill(LX);
	      RecoLY_Arm1Plane0CSC2->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 3:
	      RecoLX_Arm1Plane0CSC3->Fill(LX);
	      RecoLY_Arm1Plane0CSC3->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 4:
	      RecoLX_Arm1Plane0CSC4->Fill(LX);
	      RecoLY_Arm1Plane0CSC4->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 5:
	      RecoLX_Arm1Plane0CSC5->Fill(LX);
	      RecoLY_Arm1Plane0CSC5->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    }
	    break;
	  case 1:
	    switch(CSC){
	    case 0:
	      RecoLX_Arm1Plane1CSC0->Fill(LX);
	      RecoLY_Arm1Plane1CSC0->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 1:
	      RecoLX_Arm1Plane1CSC1->Fill(LX);
	      RecoLY_Arm1Plane1CSC1->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 2:
	      RecoLX_Arm1Plane1CSC2->Fill(LX);
	      RecoLY_Arm1Plane1CSC2->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 3:
	      RecoLX_Arm1Plane1CSC3->Fill(LX);
	      RecoLY_Arm1Plane1CSC3->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 4:
	      RecoLX_Arm1Plane1CSC4->Fill(LX);
	      RecoLY_Arm1Plane1CSC4->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 5:
	      RecoLX_Arm1Plane1CSC5->Fill(LX);
	      RecoLY_Arm1Plane1CSC5->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    }
	    break;
	  case 2:
	    switch(CSC){
	    case 0:
	      RecoLX_Arm1Plane2CSC0->Fill(LX);
	      RecoLY_Arm1Plane2CSC0->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 1:
	      RecoLX_Arm1Plane2CSC1->Fill(LX);
	      RecoLY_Arm1Plane2CSC1->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 2:
	      RecoLX_Arm1Plane2CSC2->Fill(LX);
	      RecoLY_Arm1Plane2CSC2->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 3:
	      RecoLX_Arm1Plane2CSC3->Fill(LX);
	      RecoLY_Arm1Plane2CSC3->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 4:
	      RecoLX_Arm1Plane2CSC4->Fill(LX);
	      RecoLY_Arm1Plane2CSC4->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 5:
	      RecoLX_Arm1Plane2CSC5->Fill(LX);
	      RecoLY_Arm1Plane2CSC5->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    }
	    break;
	  case 3:
	    switch(CSC){
	    case 0:
	      RecoLX_Arm1Plane3CSC0->Fill(LX);
	      RecoLY_Arm1Plane3CSC0->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 1:
	      RecoLX_Arm1Plane3CSC1->Fill(LX);
	      RecoLY_Arm1Plane3CSC1->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 2:
	      RecoLX_Arm1Plane3CSC2->Fill(LX);
	      RecoLY_Arm1Plane3CSC2->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 3:
	      RecoLX_Arm1Plane3CSC3->Fill(LX);
	      RecoLY_Arm1Plane3CSC3->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 4:
	      RecoLX_Arm1Plane3CSC4->Fill(LX);
	      RecoLY_Arm1Plane3CSC4->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 5:
	      RecoLX_Arm1Plane3CSC5->Fill(LX);
	      RecoLY_Arm1Plane3CSC5->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    }
	    break;
	  case 4:
	    switch(CSC){
	    case 0:
	      RecoLX_Arm1Plane4CSC0->Fill(LX);
	      RecoLY_Arm1Plane4CSC0->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 1:
	      RecoLX_Arm1Plane4CSC1->Fill(LX);
	      RecoLY_Arm1Plane4CSC1->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 2:
	      RecoLX_Arm1Plane4CSC2->Fill(LX);
	      RecoLY_Arm1Plane4CSC2->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 3:
	      RecoLX_Arm1Plane4CSC3->Fill(LX);
	      RecoLY_Arm1Plane4CSC3->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 4:
	      RecoLX_Arm1Plane4CSC4->Fill(LX);
	      RecoLY_Arm1Plane4CSC4->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    case 5:
	      RecoLX_Arm1Plane4CSC5->Fill(LX);
	      RecoLY_Arm1Plane4CSC5->Fill(LY);
	      NHITS[ARM][PLANE][CSC]++;    break;
	    }
	    break;
	  }
	  break;
	}






	float xxxxR=0;
	float yyyyR=0;
	float zzzzR = 0;

	layer->Local2Global((*T1Reco_it).t1DetId().rawId(),(*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y(),xxxxR,yyyyR,zzzzR);

	float R = 10000;
	float Dx = 10000;
	float Dy = 10000;
	float Dz = 10000;


	if(_SIM_>0 && bool_MinMaxTracks){
	  for(std::vector<PSimHit>::iterator hitItr = copy_of_hits->begin();
	      hitItr != copy_of_hits->end(); ++hitItr)
	    {


	      float xxxS = layer->xFromLocal2BeamSystem((*hitItr).detUnitId(),  (*hitItr).localPosition().x() );
	      float yyyS = layer->yFromLocal2BeamSystem( (*hitItr).detUnitId(),(*hitItr).localPosition().y());

	      float xxxxS=0;
	      float yyyyS=0;
	      layer->RotationLocal2Global((*hitItr).detUnitId(),xxxS,yyyS,xxxxS,yyyyS);
	      float zzzzS = layer->Zeta((*hitItr).detUnitId());
/*
  float Dx_temp = xxxxS - xxxxR;
  float Dy_temp = yyyyS - yyyyR;
  float Dz_temp = zzzzS - zzzzR;
*/
	      float Dx_temp = (*hitItr).localPosition().x() - (*T1Reco_it).localPosition().x();
	      float Dy_temp = (*hitItr).localPosition().y() - (*T1Reco_it).localPosition().y();
	      float Dz_temp = zzzzS - zzzzR;

	      float R_temp = sqrt(Dx_temp*Dx_temp + Dy_temp*Dy_temp + Dz_temp*Dz_temp);
	      if(R_temp < R && Dz_temp == 0){
		R =R_temp; 
		Dx = Dx_temp;
		Dy = Dy_temp;
		Dz = Dz_temp;
	      }
      
	    }
	}
#ifdef _PRINT_HITS_
//  std::cout << " REC HIT: " << xxxxR << " " << yyyyR << std::endl;
#endif

	if(t1Id.Arm()==0 && t1Id.Plane()==0){
	  RecoX_Arm0Pl0->Fill(xxxxR);
	  RecoY_Arm0Pl0->Fill(yyyyR);
	  RecoXY_00->Fill(xxxxR,yyyyR);
	}
	if(t1Id.Arm()==0 && t1Id.Plane()==1){
	  RecoX_Arm0Pl1->Fill(xxxxR);
	  RecoY_Arm0Pl1->Fill(yyyyR);
	  RecoXY_01->Fill(xxxxR,yyyyR);
	}
	if(t1Id.Arm()==0 && t1Id.Plane()==2){
	  RecoX_Arm0Pl2->Fill(xxxxR);
	  RecoY_Arm0Pl2->Fill(yyyyR);
	  RecoXY_02->Fill(xxxxR,yyyyR);
	}
	if(t1Id.Arm()==0 && t1Id.Plane()==3){
	  RecoX_Arm0Pl3->Fill(xxxxR);
	  RecoY_Arm0Pl3->Fill(yyyyR);
	  RecoXY_03->Fill(xxxxR,yyyyR);
	}
	if(t1Id.Arm()==0 && t1Id.Plane()==4){
	  RecoX_Arm0Pl4->Fill(xxxxR);
	  RecoY_Arm0Pl4->Fill(yyyyR);
	  RecoXY_04->Fill(xxxxR,yyyyR);
	}
	if(t1Id.Arm()==1 && t1Id.Plane()==0){
	  RecoX_Arm1Pl0->Fill(xxxxR);
	  RecoY_Arm1Pl0->Fill(yyyyR);
	  RecoXY_10->Fill(xxxxR,yyyyR);
	}
	if(t1Id.Arm()==1 && t1Id.Plane()==1){
	  RecoX_Arm1Pl1->Fill(xxxxR);
	  RecoY_Arm1Pl1->Fill(yyyyR);
	  RecoXY_11->Fill(xxxxR,yyyyR);
	}
	if(t1Id.Arm()==1 && t1Id.Plane()==2){
	  RecoX_Arm1Pl2->Fill(xxxxR);
	  RecoY_Arm1Pl2->Fill(yyyyR);
	  RecoXY_12->Fill(xxxxR,yyyyR);
	}
	if(t1Id.Arm()==1 && t1Id.Plane()==3){
	  RecoX_Arm1Pl3->Fill(xxxxR);
	  RecoY_Arm1Pl3->Fill(yyyyR);
	  RecoXY_13->Fill(xxxxR,yyyyR);
	}
	if(t1Id.Arm()==1 && t1Id.Plane()==4){
	  RecoX_Arm1Pl4->Fill(xxxxR);
	  RecoY_Arm1Pl4->Fill(yyyyR);
	  RecoXY_14->Fill(xxxxR,yyyyR);
	}

// ------------------------------

	if(Dz<10000){
	  if(t1Id.Arm()==0 && t1Id.Plane()==0 && t1Id.CSC() == 0){
	    DeltaX_Arm0Pl0->Fill(Dx);
	    DeltaY_Arm0Pl0->Fill(Dy);
	  }
	  if(t1Id.Arm()==0 && t1Id.Plane()==1 && t1Id.CSC() == 0){
	    DeltaX_Arm0Pl1->Fill(Dx);
	    DeltaY_Arm0Pl1->Fill(Dy);
	  }
	  if(t1Id.Arm()==0 && t1Id.Plane()==2 && t1Id.CSC() == 0){
	    DeltaX_Arm0Pl2->Fill(Dx);
	    DeltaY_Arm0Pl2->Fill(Dy);
	  }
	  if(t1Id.Arm()==0 && t1Id.Plane()==3 && t1Id.CSC() == 0){
	    DeltaX_Arm0Pl3->Fill(Dx);
	    DeltaY_Arm0Pl3->Fill(Dy);
	  }
	  if(t1Id.Arm()==0 && t1Id.Plane()==4 && t1Id.CSC() == 0){
	    DeltaX_Arm0Pl4->Fill(Dx);
	    DeltaY_Arm0Pl4->Fill(Dy);
	  }

	}








      }



	NReco_Arm0Plane0CSC0->Fill(NHITS[0][0][0]);
	NReco_Arm0Plane0CSC1->Fill(NHITS[0][0][1]);
	NReco_Arm0Plane0CSC2->Fill(NHITS[0][0][2]);
	NReco_Arm0Plane0CSC3->Fill(NHITS[0][0][3]);
	NReco_Arm0Plane0CSC4->Fill(NHITS[0][0][4]);
	NReco_Arm0Plane0CSC5->Fill(NHITS[0][0][5]);
	NReco_Arm0Plane1CSC0->Fill(NHITS[0][1][0]);
	NReco_Arm0Plane1CSC1->Fill(NHITS[0][1][1]);
	NReco_Arm0Plane1CSC2->Fill(NHITS[0][1][2]);
	NReco_Arm0Plane1CSC3->Fill(NHITS[0][1][3]);
	NReco_Arm0Plane1CSC4->Fill(NHITS[0][1][4]);
	NReco_Arm0Plane1CSC5->Fill(NHITS[0][1][5]);
	NReco_Arm0Plane2CSC0->Fill(NHITS[0][2][0]);
	NReco_Arm0Plane2CSC1->Fill(NHITS[0][2][1]);
	NReco_Arm0Plane2CSC2->Fill(NHITS[0][2][2]);
	NReco_Arm0Plane2CSC3->Fill(NHITS[0][2][3]);
	NReco_Arm0Plane2CSC4->Fill(NHITS[0][2][4]);
	NReco_Arm0Plane2CSC5->Fill(NHITS[0][2][5]);
	NReco_Arm0Plane3CSC0->Fill(NHITS[0][3][0]);
	NReco_Arm0Plane3CSC1->Fill(NHITS[0][3][1]);
	NReco_Arm0Plane3CSC2->Fill(NHITS[0][3][2]);
	NReco_Arm0Plane3CSC3->Fill(NHITS[0][3][3]);
	NReco_Arm0Plane3CSC4->Fill(NHITS[0][3][4]);
	NReco_Arm0Plane3CSC5->Fill(NHITS[0][3][5]);
	NReco_Arm0Plane4CSC0->Fill(NHITS[0][4][0]);
	NReco_Arm0Plane4CSC1->Fill(NHITS[0][4][1]);
	NReco_Arm0Plane4CSC2->Fill(NHITS[0][4][2]);
	NReco_Arm0Plane4CSC3->Fill(NHITS[0][4][3]);
	NReco_Arm0Plane4CSC4->Fill(NHITS[0][4][4]);
	NReco_Arm0Plane4CSC5->Fill(NHITS[0][4][5]);
	NReco_Arm1Plane0CSC0->Fill(NHITS[1][0][0]);
	NReco_Arm1Plane0CSC1->Fill(NHITS[1][0][1]);
	NReco_Arm1Plane0CSC2->Fill(NHITS[1][0][2]);
	NReco_Arm1Plane0CSC3->Fill(NHITS[1][0][3]);
	NReco_Arm1Plane0CSC4->Fill(NHITS[1][0][4]);
	NReco_Arm1Plane0CSC5->Fill(NHITS[1][0][5]);
	NReco_Arm1Plane1CSC0->Fill(NHITS[1][1][0]);
	NReco_Arm1Plane1CSC1->Fill(NHITS[1][1][1]);
	NReco_Arm1Plane1CSC2->Fill(NHITS[1][1][2]);
	NReco_Arm1Plane1CSC3->Fill(NHITS[1][1][3]);
	NReco_Arm1Plane1CSC4->Fill(NHITS[1][1][4]);
	NReco_Arm1Plane1CSC5->Fill(NHITS[1][1][5]);
	NReco_Arm1Plane2CSC0->Fill(NHITS[1][2][0]);
	NReco_Arm1Plane2CSC1->Fill(NHITS[1][2][1]);
	NReco_Arm1Plane2CSC2->Fill(NHITS[1][2][2]);
	NReco_Arm1Plane2CSC3->Fill(NHITS[1][2][3]);
	NReco_Arm1Plane2CSC4->Fill(NHITS[1][2][4]);
	NReco_Arm1Plane2CSC5->Fill(NHITS[1][2][5]);
	NReco_Arm1Plane3CSC0->Fill(NHITS[1][3][0]);
	NReco_Arm1Plane3CSC1->Fill(NHITS[1][3][1]);
	NReco_Arm1Plane3CSC2->Fill(NHITS[1][3][2]);
	NReco_Arm1Plane3CSC3->Fill(NHITS[1][3][3]);
	NReco_Arm1Plane3CSC4->Fill(NHITS[1][3][4]);
	NReco_Arm1Plane3CSC5->Fill(NHITS[1][3][5]);
	NReco_Arm1Plane4CSC0->Fill(NHITS[1][4][0]);
	NReco_Arm1Plane4CSC1->Fill(NHITS[1][4][1]);
	NReco_Arm1Plane4CSC2->Fill(NHITS[1][4][2]);
	NReco_Arm1Plane4CSC3->Fill(NHITS[1][4][3]);
	NReco_Arm1Plane4CSC4->Fill(NHITS[1][4][4]);
	NReco_Arm1Plane4CSC5->Fill(NHITS[1][4][5]);



      unsigned int NClusterCSC[2][5][6];
      bool Layer1[2][5][6];
      bool Layer2[2][5][6];
      int Width1[2][5][6];
      int Width2[2][5][6];
      bool BothLayers[2][5][6];
      unsigned int NClusterSex[2][6];
      float ClRMS[2][6];
      float ClAvg[2][6];

      for(int ai=0; ai<2; ai++)
	for(int aj=0; aj <5; aj++)
	  for(int ak=0; ak<6; ak++){
	    NClusterCSC[ai][aj][ak]=0;
	    NClusterSex[ai][ak]=0;
	    Layer1[ai][aj][ak]=false;
	    Layer2[ai][aj][ak]=false;
	    BothLayers[ai][aj][ak]=false;
	    Width1[ai][aj][ak]=0;
	    Width2[ai][aj][ak]=0;
	    ClRMS[ai][ak]=0;
	    ClAvg[ai][ak]=0;
	  }

      Handle<T1ClusterCollection> Sdigis;
      iEvent.getByLabel(t1ClusterCollectionLabel, Sdigis);

      T1ClusterCollection::DigiRangeIterator t1sdgIt;

      for (t1sdgIt = Sdigis->begin(); t1sdgIt != Sdigis->end();
	   ++t1sdgIt){
      //         std::cout << t1Id << " ------ " << (*t1sdgIt).first << std::endl;
      //	std::cout << t1Id.rawId()<< std::endl;

	int braccio = (*t1sdgIt).first.Arm();
	int piano = (*t1sdgIt).first.Plane();
	int camera = (*t1sdgIt).first.CSC(); 
	int strato = (*t1sdgIt).first.Layer(); 
	const T1ClusterCollection::Range& srange = (*t1sdgIt).second;
	  
	for (T1ClusterCollection::const_iterator sdigiIt = srange.first;
	     sdigiIt!=srange.second;++sdigiIt){
	  if(braccio == 0 && piano == 0 && camera == 0){   ilClWidthCh0Arm0Pl0->Fill( (*sdigiIt).Width() ); NClusterCSC[0][0][0]++; NClusterSex[0][0]++; if(strato==1)Layer1[0][0][0] =true; if(strato==2)Layer2[0][0][0] =true; BothLayers[0][0][0] = Layer1[0][0][0] && Layer2[0][0][0];  if(strato == 1)Width1[0][0][0] =(*sdigiIt).Width();  if(strato==2)Width2[0][0][0] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 0 && camera == 1){   ilClWidthCh1Arm0Pl0->Fill( (*sdigiIt).Width() ); NClusterCSC[0][0][1]++; NClusterSex[0][1]++; if(strato==1)Layer1[0][0][1] =true; if(strato==2)Layer2[0][0][1] =true; BothLayers[0][0][1] = Layer1[0][0][1] && Layer2[0][0][1];  if(strato == 1)Width1[0][0][1] =(*sdigiIt).Width();  if(strato==2)Width2[0][0][1] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 0 && camera == 2){   ilClWidthCh2Arm0Pl0->Fill( (*sdigiIt).Width() ); NClusterCSC[0][0][2]++; NClusterSex[0][2]++; if(strato==1)Layer1[0][0][2] =true; if(strato==2)Layer2[0][0][2] =true; BothLayers[0][0][2] = Layer1[0][0][2] && Layer2[0][0][2];  if(strato == 1)Width1[0][0][2] =(*sdigiIt).Width();  if(strato==2)Width2[0][0][2] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 0 && camera == 3){   ilClWidthCh3Arm0Pl0->Fill( (*sdigiIt).Width() ); NClusterCSC[0][0][3]++; NClusterSex[0][3]++; if(strato==1)Layer1[0][0][3] =true; if(strato==2)Layer2[0][0][3] =true; BothLayers[0][0][3] = Layer1[0][0][3] && Layer2[0][0][3];  if(strato == 1)Width1[0][0][3] =(*sdigiIt).Width();  if(strato==2)Width2[0][0][3] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 0 && camera == 4){   ilClWidthCh4Arm0Pl0->Fill( (*sdigiIt).Width() ); NClusterCSC[0][0][4]++; NClusterSex[0][4]++; if(strato==1)Layer1[0][0][4] =true; if(strato==2)Layer2[0][0][4] =true; BothLayers[0][0][4] = Layer1[0][0][4] && Layer2[0][0][4];  if(strato == 1)Width1[0][0][4] =(*sdigiIt).Width();  if(strato==2)Width2[0][0][4] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 0 && camera == 5){   ilClWidthCh5Arm0Pl0->Fill( (*sdigiIt).Width() ); NClusterCSC[0][0][5]++; NClusterSex[0][5]++; if(strato==1)Layer1[0][0][5] =true; if(strato==2)Layer2[0][0][5] =true; BothLayers[0][0][5] = Layer1[0][0][5] && Layer2[0][0][5];  if(strato == 1)Width1[0][0][5] =(*sdigiIt).Width();  if(strato==2)Width2[0][0][5] =(*sdigiIt).Width()    ;}
														       							    								   		     							 						    		
														       							    								   		     							 						    		
	  if(braccio == 0 && piano == 1 && camera == 0){   ilClWidthCh0Arm0Pl1->Fill( (*sdigiIt).Width() ); NClusterCSC[0][1][0]++; NClusterSex[0][0]++; if(strato==1)Layer1[0][1][0] =true; if(strato==2)Layer2[0][1][0] =true; BothLayers[0][1][0] = Layer1[0][1][0] && Layer2[0][1][0];  if(strato == 1)Width1[0][1][0] =(*sdigiIt).Width();  if(strato==2)Width2[0][1][0] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 1 && camera == 1){   ilClWidthCh1Arm0Pl1->Fill( (*sdigiIt).Width() ); NClusterCSC[0][1][1]++; NClusterSex[0][1]++; if(strato==1)Layer1[0][1][1] =true; if(strato==2)Layer2[0][1][1] =true; BothLayers[0][1][1] = Layer1[0][1][1] && Layer2[0][1][1];  if(strato == 1)Width1[0][1][1] =(*sdigiIt).Width();  if(strato==2)Width2[0][1][1] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 1 && camera == 2){   ilClWidthCh2Arm0Pl1->Fill( (*sdigiIt).Width() ); NClusterCSC[0][1][2]++; NClusterSex[0][2]++; if(strato==1)Layer1[0][1][2] =true; if(strato==2)Layer2[0][1][2] =true; BothLayers[0][1][2] = Layer1[0][1][2] && Layer2[0][1][2];  if(strato == 1)Width1[0][1][2] =(*sdigiIt).Width();  if(strato==2)Width2[0][1][2] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 1 && camera == 3){   ilClWidthCh3Arm0Pl1->Fill( (*sdigiIt).Width() ); NClusterCSC[0][1][3]++; NClusterSex[0][3]++; if(strato==1)Layer1[0][1][3] =true; if(strato==2)Layer2[0][1][3] =true; BothLayers[0][1][3] = Layer1[0][1][3] && Layer2[0][1][3];  if(strato == 1)Width1[0][1][3] =(*sdigiIt).Width();  if(strato==2)Width2[0][1][3] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 1 && camera == 4){   ilClWidthCh4Arm0Pl1->Fill( (*sdigiIt).Width() ); NClusterCSC[0][1][4]++; NClusterSex[0][4]++; if(strato==1)Layer1[0][1][4] =true; if(strato==2)Layer2[0][1][4] =true; BothLayers[0][1][4] = Layer1[0][1][4] && Layer2[0][1][4];  if(strato == 1)Width1[0][1][4] =(*sdigiIt).Width();  if(strato==2)Width2[0][1][4] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 1 && camera == 5){   ilClWidthCh5Arm0Pl1->Fill( (*sdigiIt).Width() ); NClusterCSC[0][1][5]++; NClusterSex[0][5]++; if(strato==1)Layer1[0][1][5] =true; if(strato==2)Layer2[0][1][5] =true; BothLayers[0][1][5] = Layer1[0][1][5] && Layer2[0][1][5];  if(strato == 1)Width1[0][1][5] =(*sdigiIt).Width();  if(strato==2)Width2[0][1][5] =(*sdigiIt).Width()    ;}
														       							    								   		     							 						    		
														       							    								   		     							 						    		
	  if(braccio == 0 && piano == 2 && camera == 0){   ilClWidthCh0Arm0Pl2->Fill( (*sdigiIt).Width() ); NClusterCSC[0][2][0]++; NClusterSex[0][0]++; if(strato==1)Layer1[0][2][0] =true; if(strato==2)Layer2[0][2][0] =true; BothLayers[0][2][0] = Layer1[0][2][0] && Layer2[0][2][0];  if(strato == 1)Width1[0][2][0] =(*sdigiIt).Width();  if(strato==2)Width2[0][2][0] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 2 && camera == 1){   ilClWidthCh1Arm0Pl2->Fill( (*sdigiIt).Width() ); NClusterCSC[0][2][1]++; NClusterSex[0][1]++; if(strato==1)Layer1[0][2][1] =true; if(strato==2)Layer2[0][2][1] =true; BothLayers[0][2][1] = Layer1[0][2][1] && Layer2[0][2][1];  if(strato == 1)Width1[0][2][1] =(*sdigiIt).Width();  if(strato==2)Width2[0][2][1] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 2 && camera == 2){   ilClWidthCh2Arm0Pl2->Fill( (*sdigiIt).Width() ); NClusterCSC[0][2][2]++; NClusterSex[0][2]++; if(strato==1)Layer1[0][2][2] =true; if(strato==2)Layer2[0][2][2] =true; BothLayers[0][2][2] = Layer1[0][2][2] && Layer2[0][2][2];  if(strato == 1)Width1[0][2][2] =(*sdigiIt).Width();  if(strato==2)Width2[0][2][2] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 2 && camera == 3){   ilClWidthCh3Arm0Pl2->Fill( (*sdigiIt).Width() ); NClusterCSC[0][2][3]++; NClusterSex[0][3]++; if(strato==1)Layer1[0][2][3] =true; if(strato==2)Layer2[0][2][3] =true; BothLayers[0][2][3] = Layer1[0][2][3] && Layer2[0][2][3];  if(strato == 1)Width1[0][2][3] =(*sdigiIt).Width();  if(strato==2)Width2[0][2][3] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 2 && camera == 4){   ilClWidthCh4Arm0Pl2->Fill( (*sdigiIt).Width() ); NClusterCSC[0][2][4]++; NClusterSex[0][4]++; if(strato==1)Layer1[0][2][4] =true; if(strato==2)Layer2[0][2][4] =true; BothLayers[0][2][4] = Layer1[0][2][4] && Layer2[0][2][4];  if(strato == 1)Width1[0][2][4] =(*sdigiIt).Width();  if(strato==2)Width2[0][2][4] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 2 && camera == 5){   ilClWidthCh5Arm0Pl2->Fill( (*sdigiIt).Width() ); NClusterCSC[0][2][5]++; NClusterSex[0][5]++; if(strato==1)Layer1[0][2][5] =true; if(strato==2)Layer2[0][2][5] =true; BothLayers[0][2][5] = Layer1[0][2][5] && Layer2[0][2][5];  if(strato == 1)Width1[0][2][5] =(*sdigiIt).Width();  if(strato==2)Width2[0][2][5] =(*sdigiIt).Width()    ;}
														       							    								   		     							 						    		
	  if(braccio == 0 && piano == 3 && camera == 0){   ilClWidthCh0Arm0Pl3->Fill( (*sdigiIt).Width() ); NClusterCSC[0][3][0]++; NClusterSex[0][0]++; if(strato==1)Layer1[0][3][0] =true; if(strato==2)Layer2[0][3][0] =true; BothLayers[0][3][0] = Layer1[0][3][0] && Layer2[0][3][0];  if(strato == 1)Width1[0][3][0] =(*sdigiIt).Width();  if(strato==2)Width2[0][3][0] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 3 && camera == 1){   ilClWidthCh1Arm0Pl3->Fill( (*sdigiIt).Width() ); NClusterCSC[0][3][1]++; NClusterSex[0][1]++; if(strato==1)Layer1[0][3][1] =true; if(strato==2)Layer2[0][3][1] =true; BothLayers[0][3][1] = Layer1[0][3][1] && Layer2[0][3][1];  if(strato == 1)Width1[0][3][1] =(*sdigiIt).Width();  if(strato==2)Width2[0][3][1] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 3 && camera == 2){   ilClWidthCh2Arm0Pl3->Fill( (*sdigiIt).Width() ); NClusterCSC[0][3][2]++; NClusterSex[0][2]++; if(strato==1)Layer1[0][3][2] =true; if(strato==2)Layer2[0][3][2] =true; BothLayers[0][3][2] = Layer1[0][3][2] && Layer2[0][3][2];  if(strato == 1)Width1[0][3][2] =(*sdigiIt).Width();  if(strato==2)Width2[0][3][2] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 3 && camera == 3){   ilClWidthCh3Arm0Pl3->Fill( (*sdigiIt).Width() ); NClusterCSC[0][3][3]++; NClusterSex[0][3]++; if(strato==1)Layer1[0][3][3] =true; if(strato==2)Layer2[0][3][3] =true; BothLayers[0][3][3] = Layer1[0][3][3] && Layer2[0][3][3];  if(strato == 1)Width1[0][3][3] =(*sdigiIt).Width();  if(strato==2)Width2[0][3][3] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 3 && camera == 4){   ilClWidthCh4Arm0Pl3->Fill( (*sdigiIt).Width() ); NClusterCSC[0][3][4]++; NClusterSex[0][4]++; if(strato==1)Layer1[0][3][4] =true; if(strato==2)Layer2[0][3][4] =true; BothLayers[0][3][4] = Layer1[0][3][4] && Layer2[0][3][4];  if(strato == 1)Width1[0][3][4] =(*sdigiIt).Width();  if(strato==2)Width2[0][3][4] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 3 && camera == 5){   ilClWidthCh5Arm0Pl3->Fill( (*sdigiIt).Width() ); NClusterCSC[0][3][5]++; NClusterSex[0][5]++; if(strato==1)Layer1[0][3][5] =true; if(strato==2)Layer2[0][3][5] =true; BothLayers[0][3][5] = Layer1[0][3][5] && Layer2[0][3][5];  if(strato == 1)Width1[0][3][5] =(*sdigiIt).Width();  if(strato==2)Width2[0][3][5] =(*sdigiIt).Width()    ;}
														       							    								   		     							 						    		
	  if(braccio == 0 && piano == 4 && camera == 0){   ilClWidthCh0Arm0Pl4->Fill( (*sdigiIt).Width() ); NClusterCSC[0][4][0]++; NClusterSex[0][0]++; if(strato==1)Layer1[0][4][0] =true; if(strato==2)Layer2[0][4][0] =true; BothLayers[0][4][0] = Layer1[0][4][0] && Layer2[0][4][0];  if(strato == 1)Width1[0][4][0] =(*sdigiIt).Width();  if(strato==2)Width2[0][4][0] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 4 && camera == 1){   ilClWidthCh1Arm0Pl4->Fill( (*sdigiIt).Width() ); NClusterCSC[0][4][1]++; NClusterSex[0][1]++; if(strato==1)Layer1[0][4][1] =true; if(strato==2)Layer2[0][4][1] =true; BothLayers[0][4][1] = Layer1[0][4][1] && Layer2[0][4][1];  if(strato == 1)Width1[0][4][1] =(*sdigiIt).Width();  if(strato==2)Width2[0][4][1] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 4 && camera == 2){   ilClWidthCh2Arm0Pl4->Fill( (*sdigiIt).Width() ); NClusterCSC[0][4][2]++; NClusterSex[0][2]++; if(strato==1)Layer1[0][4][2] =true; if(strato==2)Layer2[0][4][2] =true; BothLayers[0][4][2] = Layer1[0][4][2] && Layer2[0][4][2];  if(strato == 1)Width1[0][4][2] =(*sdigiIt).Width();  if(strato==2)Width2[0][4][2] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 4 && camera == 3){   ilClWidthCh3Arm0Pl4->Fill( (*sdigiIt).Width() ); NClusterCSC[0][4][3]++; NClusterSex[0][3]++; if(strato==1)Layer1[0][4][3] =true; if(strato==2)Layer2[0][4][3] =true; BothLayers[0][4][3] = Layer1[0][4][3] && Layer2[0][4][3];  if(strato == 1)Width1[0][4][3] =(*sdigiIt).Width();  if(strato==2)Width2[0][4][3] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 4 && camera == 4){   ilClWidthCh4Arm0Pl4->Fill( (*sdigiIt).Width() ); NClusterCSC[0][4][4]++; NClusterSex[0][4]++; if(strato==1)Layer1[0][4][4] =true; if(strato==2)Layer2[0][4][4] =true; BothLayers[0][4][4] = Layer1[0][4][4] && Layer2[0][4][4];  if(strato == 1)Width1[0][4][4] =(*sdigiIt).Width();  if(strato==2)Width2[0][4][4] =(*sdigiIt).Width()    ;}
	  if(braccio == 0 && piano == 4 && camera == 5){   ilClWidthCh5Arm0Pl4->Fill( (*sdigiIt).Width() ); NClusterCSC[0][4][5]++; NClusterSex[0][5]++; if(strato==1)Layer1[0][4][5] =true; if(strato==2)Layer2[0][4][5] =true; BothLayers[0][4][5] = Layer1[0][4][5] && Layer2[0][4][5];  if(strato == 1)Width1[0][4][5] =(*sdigiIt).Width();  if(strato==2)Width2[0][4][5] =(*sdigiIt).Width()    ;}
														       							    								   		     							 						    																																																																						
														       							    								   		     							 						    																																																																						
														       							    								   		     							 						    																																																																						
														       							    								   		     							 						    																																																																						
	  if(braccio == 1 && piano == 0 && camera == 0){   ilClWidthCh0Arm1Pl0->Fill( (*sdigiIt).Width() ); NClusterCSC[1][0][0]++; NClusterSex[1][0]++; if(strato==1)Layer1[1][0][0] =true; if(strato==2)Layer2[1][0][0] =true; BothLayers[1][0][0] = Layer1[1][0][0] && Layer2[1][0][0];  if(strato == 1)Width1[1][0][0] =(*sdigiIt).Width();  if(strato==2)Width2[1][0][0] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 0 && camera == 1){   ilClWidthCh1Arm1Pl0->Fill( (*sdigiIt).Width() ); NClusterCSC[1][0][1]++; NClusterSex[1][1]++; if(strato==1)Layer1[1][0][1] =true; if(strato==2)Layer2[1][0][1] =true; BothLayers[1][0][1] = Layer1[1][0][1] && Layer2[1][0][1];  if(strato == 1)Width1[1][0][1] =(*sdigiIt).Width();  if(strato==2)Width2[1][0][1] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 0 && camera == 2){   ilClWidthCh2Arm1Pl0->Fill( (*sdigiIt).Width() ); NClusterCSC[1][0][2]++; NClusterSex[1][2]++; if(strato==1)Layer1[1][0][2] =true; if(strato==2)Layer2[1][0][2] =true; BothLayers[1][0][2] = Layer1[1][0][2] && Layer2[1][0][2];  if(strato == 1)Width1[1][0][2] =(*sdigiIt).Width();  if(strato==2)Width2[1][0][2] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 0 && camera == 3){   ilClWidthCh3Arm1Pl0->Fill( (*sdigiIt).Width() ); NClusterCSC[1][0][3]++; NClusterSex[1][3]++; if(strato==1)Layer1[1][0][3] =true; if(strato==2)Layer2[1][0][3] =true; BothLayers[1][0][3] = Layer1[1][0][3] && Layer2[1][0][3];  if(strato == 1)Width1[1][0][3] =(*sdigiIt).Width();  if(strato==2)Width2[1][0][3] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 0 && camera == 4){   ilClWidthCh4Arm1Pl0->Fill( (*sdigiIt).Width() ); NClusterCSC[1][0][4]++; NClusterSex[1][4]++; if(strato==1)Layer1[1][0][4] =true; if(strato==2)Layer2[1][0][4] =true; BothLayers[1][0][4] = Layer1[1][0][4] && Layer2[1][0][4];  if(strato == 1)Width1[1][0][4] =(*sdigiIt).Width();  if(strato==2)Width2[1][0][4] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 0 && camera == 5){   ilClWidthCh5Arm1Pl0->Fill( (*sdigiIt).Width() ); NClusterCSC[1][0][5]++; NClusterSex[1][5]++; if(strato==1)Layer1[1][0][5] =true; if(strato==2)Layer2[1][0][5] =true; BothLayers[1][0][5] = Layer1[1][0][5] && Layer2[1][0][5];  if(strato == 1)Width1[1][0][5] =(*sdigiIt).Width();  if(strato==2)Width2[1][0][5] =(*sdigiIt).Width()    ;}
														       							    								   		     							 						    																																																																						
														       							    								   		     							 						    																																																																						
	  if(braccio == 1 && piano == 1 && camera == 0){   ilClWidthCh0Arm1Pl1->Fill( (*sdigiIt).Width() ); NClusterCSC[1][1][0]++; NClusterSex[1][0]++; if(strato==1)Layer1[1][1][0] =true; if(strato==2)Layer2[1][1][0] =true; BothLayers[1][1][0] = Layer1[1][1][0] && Layer2[1][1][0];  if(strato == 1)Width1[1][1][0] =(*sdigiIt).Width();  if(strato==2)Width2[1][1][0] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 1 && camera == 1){   ilClWidthCh1Arm1Pl1->Fill( (*sdigiIt).Width() ); NClusterCSC[1][1][1]++; NClusterSex[1][1]++; if(strato==1)Layer1[1][1][1] =true; if(strato==2)Layer2[1][1][1] =true; BothLayers[1][1][1] = Layer1[1][1][1] && Layer2[1][1][1];  if(strato == 1)Width1[1][1][1] =(*sdigiIt).Width();  if(strato==2)Width2[1][1][1] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 1 && camera == 2){   ilClWidthCh2Arm1Pl1->Fill( (*sdigiIt).Width() ); NClusterCSC[1][1][2]++; NClusterSex[1][2]++; if(strato==1)Layer1[1][1][2] =true; if(strato==2)Layer2[1][1][2] =true; BothLayers[1][1][2] = Layer1[1][1][2] && Layer2[1][1][2];  if(strato == 1)Width1[1][1][2] =(*sdigiIt).Width();  if(strato==2)Width2[1][1][2] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 1 && camera == 3){   ilClWidthCh3Arm1Pl1->Fill( (*sdigiIt).Width() ); NClusterCSC[1][1][3]++; NClusterSex[1][3]++; if(strato==1)Layer1[1][1][3] =true; if(strato==2)Layer2[1][1][3] =true; BothLayers[1][1][3] = Layer1[1][1][3] && Layer2[1][1][3];  if(strato == 1)Width1[1][1][3] =(*sdigiIt).Width();  if(strato==2)Width2[1][1][3] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 1 && camera == 4){   ilClWidthCh4Arm1Pl1->Fill( (*sdigiIt).Width() ); NClusterCSC[1][1][4]++; NClusterSex[1][4]++; if(strato==1)Layer1[1][1][4] =true; if(strato==2)Layer2[1][1][4] =true; BothLayers[1][1][4] = Layer1[1][1][4] && Layer2[1][1][4];  if(strato == 1)Width1[1][1][4] =(*sdigiIt).Width();  if(strato==2)Width2[1][1][4] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 1 && camera == 5){   ilClWidthCh5Arm1Pl1->Fill( (*sdigiIt).Width() ); NClusterCSC[1][1][5]++; NClusterSex[1][5]++; if(strato==1)Layer1[1][1][5] =true; if(strato==2)Layer2[1][1][5] =true; BothLayers[1][1][5] = Layer1[1][1][5] && Layer2[1][1][5];  if(strato == 1)Width1[1][1][5] =(*sdigiIt).Width();  if(strato==2)Width2[1][1][5] =(*sdigiIt).Width()    ;}
														       							    								   		     							 						    																																																																						
														       							    								   		     							 						    																																																																						
	  if(braccio == 1 && piano == 2 && camera == 0){   ilClWidthCh0Arm1Pl2->Fill( (*sdigiIt).Width() ); NClusterCSC[1][2][0]++; NClusterSex[1][0]++; if(strato==1)Layer1[1][2][0] =true; if(strato==2)Layer2[1][2][0] =true; BothLayers[1][2][0] = Layer1[1][2][0] && Layer2[1][2][0];  if(strato == 1)Width1[1][2][0] =(*sdigiIt).Width();  if(strato==2)Width2[1][2][0] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 2 && camera == 1){   ilClWidthCh1Arm1Pl2->Fill( (*sdigiIt).Width() ); NClusterCSC[1][2][1]++; NClusterSex[1][1]++; if(strato==1)Layer1[1][2][1] =true; if(strato==2)Layer2[1][2][1] =true; BothLayers[1][2][1] = Layer1[1][2][1] && Layer2[1][2][1];  if(strato == 1)Width1[1][2][1] =(*sdigiIt).Width();  if(strato==2)Width2[1][2][1] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 2 && camera == 2){   ilClWidthCh2Arm1Pl2->Fill( (*sdigiIt).Width() ); NClusterCSC[1][2][2]++; NClusterSex[1][2]++; if(strato==1)Layer1[1][2][2] =true; if(strato==2)Layer2[1][2][2] =true; BothLayers[1][2][2] = Layer1[1][2][2] && Layer2[1][2][2];  if(strato == 1)Width1[1][2][2] =(*sdigiIt).Width();  if(strato==2)Width2[1][2][2] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 2 && camera == 3){   ilClWidthCh3Arm1Pl2->Fill( (*sdigiIt).Width() ); NClusterCSC[1][2][3]++; NClusterSex[1][3]++; if(strato==1)Layer1[1][2][3] =true; if(strato==2)Layer2[1][2][3] =true; BothLayers[1][2][3] = Layer1[1][2][3] && Layer2[1][2][3];  if(strato == 1)Width1[1][2][3] =(*sdigiIt).Width();  if(strato==2)Width2[1][2][3] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 2 && camera == 4){   ilClWidthCh4Arm1Pl2->Fill( (*sdigiIt).Width() ); NClusterCSC[1][2][4]++; NClusterSex[1][4]++; if(strato==1)Layer1[1][2][4] =true; if(strato==2)Layer2[1][2][4] =true; BothLayers[1][2][4] = Layer1[1][2][4] && Layer2[1][2][4];  if(strato == 1)Width1[1][2][4] =(*sdigiIt).Width();  if(strato==2)Width2[1][2][4] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 2 && camera == 5){   ilClWidthCh5Arm1Pl2->Fill( (*sdigiIt).Width() ); NClusterCSC[1][2][5]++; NClusterSex[1][5]++; if(strato==1)Layer1[1][2][5] =true; if(strato==2)Layer2[1][2][5] =true; BothLayers[1][2][5] = Layer1[1][2][5] && Layer2[1][2][5];  if(strato == 1)Width1[1][2][5] =(*sdigiIt).Width();  if(strato==2)Width2[1][2][5] =(*sdigiIt).Width()    ;}
														       							    								   		     							 						    																																																																						
	  if(braccio == 1 && piano == 3 && camera == 0){   ilClWidthCh0Arm1Pl3->Fill( (*sdigiIt).Width() ); NClusterCSC[1][3][0]++; NClusterSex[1][0]++; if(strato==1)Layer1[1][3][0] =true; if(strato==2)Layer2[1][3][0] =true; BothLayers[1][3][0] = Layer1[1][3][0] && Layer2[1][3][0];  if(strato == 1)Width1[1][3][0] =(*sdigiIt).Width();  if(strato==2)Width2[1][3][0] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 3 && camera == 1){   ilClWidthCh1Arm1Pl3->Fill( (*sdigiIt).Width() ); NClusterCSC[1][3][1]++; NClusterSex[1][1]++; if(strato==1)Layer1[1][3][1] =true; if(strato==2)Layer2[1][3][1] =true; BothLayers[1][3][1] = Layer1[1][3][1] && Layer2[1][3][1];  if(strato == 1)Width1[1][3][1] =(*sdigiIt).Width();  if(strato==2)Width2[1][3][1] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 3 && camera == 2){   ilClWidthCh2Arm1Pl3->Fill( (*sdigiIt).Width() ); NClusterCSC[1][3][2]++; NClusterSex[1][2]++; if(strato==1)Layer1[1][3][2] =true; if(strato==2)Layer2[1][3][2] =true; BothLayers[1][3][2] = Layer1[1][3][2] && Layer2[1][3][2];  if(strato == 1)Width1[1][3][2] =(*sdigiIt).Width();  if(strato==2)Width2[1][3][2] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 3 && camera == 3){   ilClWidthCh3Arm1Pl3->Fill( (*sdigiIt).Width() ); NClusterCSC[1][3][3]++; NClusterSex[1][3]++; if(strato==1)Layer1[1][3][3] =true; if(strato==2)Layer2[1][3][3] =true; BothLayers[1][3][3] = Layer1[1][3][3] && Layer2[1][3][3];  if(strato == 1)Width1[1][3][3] =(*sdigiIt).Width();  if(strato==2)Width2[1][3][3] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 3 && camera == 4){   ilClWidthCh4Arm1Pl3->Fill( (*sdigiIt).Width() ); NClusterCSC[1][3][4]++; NClusterSex[1][4]++; if(strato==1)Layer1[1][3][4] =true; if(strato==2)Layer2[1][3][4] =true; BothLayers[1][3][4] = Layer1[1][3][4] && Layer2[1][3][4];  if(strato == 1)Width1[1][3][4] =(*sdigiIt).Width();  if(strato==2)Width2[1][3][4] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 3 && camera == 5){   ilClWidthCh5Arm1Pl3->Fill( (*sdigiIt).Width() ); NClusterCSC[1][3][5]++; NClusterSex[1][5]++; if(strato==1)Layer1[1][3][5] =true; if(strato==2)Layer2[1][3][5] =true; BothLayers[1][3][5] = Layer1[1][3][5] && Layer2[1][3][5];  if(strato == 1)Width1[1][3][5] =(*sdigiIt).Width();  if(strato==2)Width2[1][3][5] =(*sdigiIt).Width()    ;}
														       							    								   		     							 						    																																																																						
	  if(braccio == 1 && piano == 4 && camera == 0){   ilClWidthCh0Arm1Pl4->Fill( (*sdigiIt).Width() ); NClusterCSC[1][4][0]++; NClusterSex[1][0]++; if(strato==1)Layer1[1][4][0] =true; if(strato==2)Layer2[1][4][0] =true; BothLayers[1][4][0] = Layer1[1][4][0] && Layer2[1][4][0];  if(strato == 1)Width1[1][4][0] =(*sdigiIt).Width();  if(strato==2)Width2[1][4][0] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 4 && camera == 1){   ilClWidthCh1Arm1Pl4->Fill( (*sdigiIt).Width() ); NClusterCSC[1][4][1]++; NClusterSex[1][1]++; if(strato==1)Layer1[1][4][1] =true; if(strato==2)Layer2[1][4][1] =true; BothLayers[1][4][1] = Layer1[1][4][1] && Layer2[1][4][1];  if(strato == 1)Width1[1][4][1] =(*sdigiIt).Width();  if(strato==2)Width2[1][4][1] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 4 && camera == 2){   ilClWidthCh2Arm1Pl4->Fill( (*sdigiIt).Width() ); NClusterCSC[1][4][2]++; NClusterSex[1][2]++; if(strato==1)Layer1[1][4][2] =true; if(strato==2)Layer2[1][4][2] =true; BothLayers[1][4][2] = Layer1[1][4][2] && Layer2[1][4][2];  if(strato == 1)Width1[1][4][2] =(*sdigiIt).Width();  if(strato==2)Width2[1][4][2] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 4 && camera == 3){   ilClWidthCh3Arm1Pl4->Fill( (*sdigiIt).Width() ); NClusterCSC[1][4][3]++; NClusterSex[1][3]++; if(strato==1)Layer1[1][4][3] =true; if(strato==2)Layer2[1][4][3] =true; BothLayers[1][4][3] = Layer1[1][4][3] && Layer2[1][4][3];  if(strato == 1)Width1[1][4][3] =(*sdigiIt).Width();  if(strato==2)Width2[1][4][3] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 4 && camera == 4){   ilClWidthCh4Arm1Pl4->Fill( (*sdigiIt).Width() ); NClusterCSC[1][4][4]++; NClusterSex[1][4]++; if(strato==1)Layer1[1][4][4] =true; if(strato==2)Layer2[1][4][4] =true; BothLayers[1][4][4] = Layer1[1][4][4] && Layer2[1][4][4];  if(strato == 1)Width1[1][4][4] =(*sdigiIt).Width();  if(strato==2)Width2[1][4][4] =(*sdigiIt).Width()    ;}
	  if(braccio == 1 && piano == 4 && camera == 5){   ilClWidthCh5Arm1Pl4->Fill( (*sdigiIt).Width() ); NClusterCSC[1][4][5]++; NClusterSex[1][5]++; if(strato==1)Layer1[1][4][5] =true; if(strato==2)Layer2[1][4][5] =true; BothLayers[1][4][5] = Layer1[1][4][5] && Layer2[1][4][5];  if(strato == 1)Width1[1][4][5] =(*sdigiIt).Width();  if(strato==2)Width2[1][4][5] =(*sdigiIt).Width()    ;}







	}
      }

// compare cluster width of different layers when only 1 cluster per layer
      for(int ai=0; ai<2; ai++)
	for(int aj=0; aj <5; aj++)
	  for(int ak=0; ak<6; ak++){
	    if( NClusterCSC[ai][aj][ak]==2 && BothLayers[ai][aj][ak] && Width1[ai][aj][ak]>0 && Width2[ai][aj][ak]>0 )hClDiffBetweenLayers->Fill(Width1[ai][aj][ak]-Width2[ai][aj][ak]);

	  }
      for(int ai=0; ai<2; ai++)
	for(int ak=0; ak<6; ak++){
	  if(NClusterSex[ai][ak]>=6 && NClusterCSC[ai][0][ak]<=2 && NClusterCSC[ai][1][ak]<=2 && NClusterCSC[ai][2][ak]<=2 && NClusterCSC[ai][3][ak]<=2 && NClusterCSC[ai][4][ak]<=2 ){

	    ClAvg[ai][ak] = (float)(Width1[ai][0][ak] +Width1[ai][1][ak] +Width1[ai][2][ak] +Width1[ai][3][ak] +Width1[ai][4][ak] + Width2[ai][0][ak] +Width2[ai][1][ak] +Width2[ai][2][ak] +Width2[ai][3][ak] +Width2[ai][4][ak])/(float)(NClusterCSC[ai][0][ak] +NClusterCSC[ai][1][ak] +NClusterCSC[ai][2][ak] +NClusterCSC[ai][3][ak] +NClusterCSC[ai][4][ak]);

	    if(Layer1[ai][0][ak])ClRMS[ai][ak] += (float)(Width1[ai][0][ak] - ClAvg[ai][ak])*(Width1[ai][0][ak] - ClAvg[ai][ak])/(float)(NClusterCSC[ai][0][ak] +NClusterCSC[ai][1][ak] +NClusterCSC[ai][2][ak] +NClusterCSC[ai][3][ak] +NClusterCSC[ai][4][ak]);
	    if(Layer1[ai][1][ak])ClRMS[ai][ak] += (float)(Width1[ai][1][ak] - ClAvg[ai][ak])*(Width1[ai][1][ak] - ClAvg[ai][ak])/(float)(NClusterCSC[ai][0][ak] +NClusterCSC[ai][1][ak] +NClusterCSC[ai][2][ak] +NClusterCSC[ai][3][ak] +NClusterCSC[ai][4][ak]);
	    if(Layer1[ai][2][ak])ClRMS[ai][ak] += (float)(Width1[ai][2][ak] - ClAvg[ai][ak])*(Width1[ai][2][ak] - ClAvg[ai][ak])/(float)(NClusterCSC[ai][0][ak] +NClusterCSC[ai][1][ak] +NClusterCSC[ai][2][ak] +NClusterCSC[ai][3][ak] +NClusterCSC[ai][4][ak]);
	    if(Layer1[ai][3][ak])ClRMS[ai][ak] += (float)(Width1[ai][3][ak] - ClAvg[ai][ak])*(Width1[ai][3][ak] - ClAvg[ai][ak])/(float)(NClusterCSC[ai][0][ak] +NClusterCSC[ai][1][ak] +NClusterCSC[ai][2][ak] +NClusterCSC[ai][3][ak] +NClusterCSC[ai][4][ak]);
	    if(Layer1[ai][4][ak])ClRMS[ai][ak] += (float)(Width1[ai][4][ak] - ClAvg[ai][ak])*(Width1[ai][4][ak] - ClAvg[ai][ak])/(float)(NClusterCSC[ai][0][ak] +NClusterCSC[ai][1][ak] +NClusterCSC[ai][2][ak] +NClusterCSC[ai][3][ak] +NClusterCSC[ai][4][ak]);

	    if(Layer2[ai][0][ak])ClRMS[ai][ak] += (float)(Width2[ai][0][ak] - ClAvg[ai][ak])*(Width2[ai][0][ak] - ClAvg[ai][ak])/(float)(NClusterCSC[ai][0][ak] +NClusterCSC[ai][1][ak] +NClusterCSC[ai][2][ak] +NClusterCSC[ai][3][ak] +NClusterCSC[ai][4][ak]);
	    if(Layer2[ai][1][ak])ClRMS[ai][ak] += (float)(Width2[ai][1][ak] - ClAvg[ai][ak])*(Width2[ai][1][ak] - ClAvg[ai][ak])/(float)(NClusterCSC[ai][0][ak] +NClusterCSC[ai][1][ak] +NClusterCSC[ai][2][ak] +NClusterCSC[ai][3][ak] +NClusterCSC[ai][4][ak]);
	    if(Layer2[ai][2][ak])ClRMS[ai][ak] += (float)(Width2[ai][2][ak] - ClAvg[ai][ak])*(Width2[ai][2][ak] - ClAvg[ai][ak])/(float)(NClusterCSC[ai][0][ak] +NClusterCSC[ai][1][ak] +NClusterCSC[ai][2][ak] +NClusterCSC[ai][3][ak] +NClusterCSC[ai][4][ak]);
	    if(Layer2[ai][3][ak])ClRMS[ai][ak] += (float)(Width2[ai][3][ak] - ClAvg[ai][ak])*(Width2[ai][3][ak] - ClAvg[ai][ak])/(float)(NClusterCSC[ai][0][ak] +NClusterCSC[ai][1][ak] +NClusterCSC[ai][2][ak] +NClusterCSC[ai][3][ak] +NClusterCSC[ai][4][ak]);
	    if(Layer2[ai][4][ak])ClRMS[ai][ak] += (float)(Width2[ai][4][ak] - ClAvg[ai][ak])*(Width2[ai][4][ak] - ClAvg[ai][ak])/(float)(NClusterCSC[ai][0][ak] +NClusterCSC[ai][1][ak] +NClusterCSC[ai][2][ak] +NClusterCSC[ai][3][ak] +NClusterCSC[ai][4][ak]);




	    ClRMS[ai][ak] = sqrt(ClRMS[ai][ak]);


	    hClRMSSex -> Fill(	ClRMS[ai][ak] );

	  }
	}











      edm::Handle<T1RoadCollection> roadCollection;
      iEvent.getByLabel(t1RoadCollectionLabel, roadCollection);
      hRoadNumber->Fill(roadCollection->size());
      T1RoadCollection::const_iterator RC_it;

      for(RC_it=roadCollection->begin(); RC_it!=roadCollection->end(); RC_it++){
	hRoadSize->Fill((*RC_it).size() );
	for(unsigned int i=0; i<(*RC_it).size(); i++){

	  hRoadHitSigmaX->Fill(sqrt((*RC_it)[i].GlobalPositionError().cxx())); 
	  hRoadHitSigmaY->Fill(sqrt((*RC_it)[i].GlobalPositionError().cyy())); 
	}
      }


    }

#endif
  }











  if(_SIM_>0 && bool_MinMaxTracks){

    int i_myColl=0;

    std::vector<edm::Handle<std::vector<PSimHit> > > resultsim;
    iEvent.getManyByType(resultsim);

    int ss=resultsim.size();
    for (int ii=0;ii<ss;ii++) {
      edm::BranchDescription desc = resultsim[ii].provenance()->product();
      edm::LogInfo("T1DP") <<"For "<<desc.productInstanceName()<<" " <<resultsim[ii].product()->size()<<" Simhits added";
      if(desc.productInstanceName()=="TotemHitsT1")
        i_myColl=ii;
    }

    const std::vector<PSimHit> * simHits=resultsim[i_myColl].product();

    for(std::vector<PSimHit>::const_iterator hitItr = simHits->begin();
        hitItr != simHits->end(); ++hitItr)
      {

	copy_of_hits->push_back(  (*hitItr)  );

        ilMomento->Fill( (*hitItr).pabs() );
        ilMomentoLog->Fill( log10( (*hitItr).pabs() ) );
        lEnergia->Fill( (*hitItr).energyLoss() );
        ilTipo->Fill( (*hitItr).particleType() );
        ilDetector->Fill( (*hitItr).detUnitId() );

        float px = (*hitItr).localDirection().x() * (*hitItr).pabs();
        float py = (*hitItr).localDirection().y() * (*hitItr).pabs();
        float pt = sqrt(px*px+py*py);

        if((*hitItr).detUnitId()!=0){
          T1DetId oggetto((uint32_t)(*hitItr).detUnitId());
        }else{
          std::cout << " WRONG DETID " << std::endl;
        }

      // calcolo l'Eta di tutte le particelle
        lEta->Fill(layer->eta( (*hitItr).detUnitId(),  (*hitItr).localPosition().x(),(*hitItr).localPosition().y() ));

        ilPT -> Fill(pt);


#ifdef _PRINT_HITS_
	float xxxR = layer->xFromLocal2BeamSystem((*hitItr).detUnitId(),  (*hitItr).localPosition().x() );
	float yyyR = layer->yFromLocal2BeamSystem( (*hitItr).detUnitId(),(*hitItr).localPosition().y());

	float xxxxR=0;
	float yyyyR=0;
	layer->RotationLocal2Global((*hitItr).detUnitId(),xxxR,yyyR,xxxxR,yyyyR);

      // std::cout << " SIM HIT: " << xxxxR << " " << yyyyR << std::endl;
#endif



      }




    lEta->GetXaxis()->SetTitle("#eta");
    ilPT->GetXaxis()->SetTitle("p_T");
  }

#ifdef _ALSO_DIGIS_
  if(_DIGI_>0 && bool_MinMaxTracks){


    edm::Handle<T1DigiWireCollection> myDigiColl;
    iEvent.getByLabel(t1DigiWireCollectionLabel, myDigiColl);

    edm::Handle<T1DigiVfatCollection> myDigiCollStrip;
    iEvent.getByLabel(t1DigiVfatCollectionLabel, myDigiCollStrip);

    T1DigiWireCollection::DigiRangeIterator T1Digi_it;
    for(T1Digi_it=myDigiColl->begin();T1Digi_it!=myDigiColl->end();T1Digi_it++){

      int braccio = (*T1Digi_it).first.Arm();
      int piano = (*T1Digi_it).first.Plane();
      int camera = (*T1Digi_it).first.CSC();

      const T1DigiWireCollection::Range& range = (*T1Digi_it).second;
      for (T1DigiWireCollection::const_iterator digiIt = range.first;
           digiIt!=range.second;++digiIt){
#ifdef _PRINT_
        (*digiIt).print();
#endif

        if(piano == 0 && braccio == 0){
	  switch (camera){
	  case 0:
	    ilWireCh0Arm0Pl0->Fill(digiIt->wire());
	    ilWireCh0Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    ilWireCh1Arm0Pl0->Fill(digiIt->wire());
	    ilWireCh1Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    ilWireCh2Arm0Pl0->Fill(digiIt->wire());
	    ilWireCh2Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    ilWireCh3Arm0Pl0->Fill(digiIt->wire());
	    ilWireCh3Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    ilWireCh4Arm0Pl0->Fill(digiIt->wire());
	    ilWireCh4Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    ilWireCh5Arm0Pl0->Fill(digiIt->wire());
	    ilWireCh5Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 1 && braccio == 0){
	  switch (camera){
	  case 0:
	    ilWireCh0Arm0Pl1->Fill(digiIt->wire());
	    ilWireCh0Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    ilWireCh1Arm0Pl1->Fill(digiIt->wire());
	    ilWireCh1Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    ilWireCh2Arm0Pl1->Fill(digiIt->wire());
	    ilWireCh2Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    ilWireCh3Arm0Pl1->Fill(digiIt->wire());
	    ilWireCh3Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    ilWireCh4Arm0Pl1->Fill(digiIt->wire());
	    ilWireCh4Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    ilWireCh5Arm0Pl1->Fill(digiIt->wire());
	    ilWireCh5Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 2 && braccio == 0){
	  switch (camera){
	  case 0:
	    ilWireCh0Arm0Pl2->Fill(digiIt->wire());
	    ilWireCh0Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    ilWireCh1Arm0Pl2->Fill(digiIt->wire());
	    ilWireCh1Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    ilWireCh2Arm0Pl2->Fill(digiIt->wire());
	    ilWireCh2Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    ilWireCh3Arm0Pl2->Fill(digiIt->wire());
	    ilWireCh3Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    ilWireCh4Arm0Pl2->Fill(digiIt->wire());
	    ilWireCh4Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    ilWireCh5Arm0Pl2->Fill(digiIt->wire());
	    ilWireCh5Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 3 && braccio == 0){
	  switch (camera){
	  case 0:
	    ilWireCh0Arm0Pl3->Fill(digiIt->wire());
	    ilWireCh0Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    ilWireCh1Arm0Pl3->Fill(digiIt->wire());
	    ilWireCh1Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    ilWireCh2Arm0Pl3->Fill(digiIt->wire());
	    ilWireCh2Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    ilWireCh3Arm0Pl3->Fill(digiIt->wire());
	    ilWireCh3Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    ilWireCh4Arm0Pl3->Fill(digiIt->wire());
	    ilWireCh4Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    ilWireCh5Arm0Pl3->Fill(digiIt->wire());
	    ilWireCh5Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 4 && braccio == 0){
	  switch (camera){
	  case 0:
	    ilWireCh0Arm0Pl4->Fill(digiIt->wire());
	    ilWireCh0Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    ilWireCh1Arm0Pl4->Fill(digiIt->wire());
	    ilWireCh1Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    ilWireCh2Arm0Pl4->Fill(digiIt->wire());
	    ilWireCh2Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    ilWireCh3Arm0Pl4->Fill(digiIt->wire());
	    ilWireCh3Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    ilWireCh4Arm0Pl4->Fill(digiIt->wire());
	    ilWireCh4Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    ilWireCh5Arm0Pl4->Fill(digiIt->wire());
	    ilWireCh5Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

      ///////

        if(piano == 0 && braccio == 1){
	  switch (camera){
	  case 0:
	    ilWireCh0Arm1Pl0->Fill(digiIt->wire());
	    ilWireCh0Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    ilWireCh1Arm1Pl0->Fill(digiIt->wire());
	    ilWireCh1Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    ilWireCh2Arm1Pl0->Fill(digiIt->wire());
	    ilWireCh2Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    ilWireCh3Arm1Pl0->Fill(digiIt->wire());
	    ilWireCh3Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    ilWireCh4Arm1Pl0->Fill(digiIt->wire());
	    ilWireCh4Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    ilWireCh5Arm1Pl0->Fill(digiIt->wire());
	    ilWireCh5Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 1 && braccio == 1){
	  switch (camera){
	  case 0:
	    ilWireCh0Arm1Pl1->Fill(digiIt->wire());
	    ilWireCh0Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    ilWireCh1Arm1Pl1->Fill(digiIt->wire());
	    ilWireCh1Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    ilWireCh2Arm1Pl1->Fill(digiIt->wire());
	    ilWireCh2Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    ilWireCh3Arm1Pl1->Fill(digiIt->wire());
	    ilWireCh3Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    ilWireCh4Arm1Pl1->Fill(digiIt->wire());
	    ilWireCh4Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    ilWireCh5Arm1Pl1->Fill(digiIt->wire());
	    ilWireCh5Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 2 && braccio == 1){
	  switch (camera){
	  case 0:
	    ilWireCh0Arm1Pl2->Fill(digiIt->wire());
	    ilWireCh0Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    ilWireCh1Arm1Pl2->Fill(digiIt->wire());
	    ilWireCh1Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    ilWireCh2Arm1Pl2->Fill(digiIt->wire());
	    ilWireCh2Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    ilWireCh3Arm1Pl2->Fill(digiIt->wire());
	    ilWireCh3Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    ilWireCh4Arm1Pl2->Fill(digiIt->wire());
	    ilWireCh4Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    ilWireCh5Arm1Pl2->Fill(digiIt->wire());
	    ilWireCh5Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 3 && braccio == 1){
	  switch (camera){
	  case 0:
	    ilWireCh0Arm1Pl3->Fill(digiIt->wire());
	    ilWireCh0Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    ilWireCh1Arm1Pl3->Fill(digiIt->wire());
	    ilWireCh1Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    ilWireCh2Arm1Pl3->Fill(digiIt->wire());
	    ilWireCh2Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    ilWireCh3Arm1Pl3->Fill(digiIt->wire());
	    ilWireCh3Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    ilWireCh4Arm1Pl3->Fill(digiIt->wire());
	    ilWireCh4Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    ilWireCh5Arm1Pl3->Fill(digiIt->wire());
	    ilWireCh5Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 4 && braccio == 1){
	  switch (camera){
	  case 0:
	    ilWireCh0Arm1Pl4->Fill(digiIt->wire());
	    ilWireCh0Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    ilWireCh1Arm1Pl4->Fill(digiIt->wire());
	    ilWireCh1Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    ilWireCh2Arm1Pl4->Fill(digiIt->wire());
	    ilWireCh2Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    ilWireCh3Arm1Pl4->Fill(digiIt->wire());
	    ilWireCh3Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    ilWireCh4Arm1Pl4->Fill(digiIt->wire());
	    ilWireCh4Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    ilWireCh5Arm1Pl4->Fill(digiIt->wire());
	    ilWireCh5Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;

	  }
        }
/////////////////////////////////////////////////

        if(piano == 0 && braccio == 0){
	  switch (camera){
	  case 0:
          // _NWireCh0Arm0Pl0++ ; //(digiIt++ ; //wire());
	    _NWireCh0Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NWireCh1Arm0Pl0++ ; //(digiIt++ ; //wire());
	    _NWireCh1Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NWireCh2Arm0Pl0++ ; //(digiIt++ ; //wire());
	    _NWireCh2Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NWireCh3Arm0Pl0++ ; //(digiIt++ ; //wire());
	    _NWireCh3Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NWireCh4Arm0Pl0++ ; //(digiIt++ ; //wire());
	    _NWireCh4Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NWireCh5Arm0Pl0++ ; //(digiIt++ ; //wire());
	    _NWireCh5Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 1 && braccio == 0){
	  switch (camera){
	  case 0:
          // _NWireCh0Arm0Pl1++ ; //(digiIt++ ; //wire());
	    _NWireCh0Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NWireCh1Arm0Pl1++ ; //(digiIt++ ; //wire());
	    _NWireCh1Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NWireCh2Arm0Pl1++ ; //(digiIt++ ; //wire());
	    _NWireCh2Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NWireCh3Arm0Pl1++ ; //(digiIt++ ; //wire());
	    _NWireCh3Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NWireCh4Arm0Pl1++ ; //(digiIt++ ; //wire());
	    _NWireCh4Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NWireCh5Arm0Pl1++ ; //(digiIt++ ; //wire());
	    _NWireCh5Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 2 && braccio == 0){
	  switch (camera){
	  case 0:
          // _NWireCh0Arm0Pl2++ ; //(digiIt++ ; //wire());
	    _NWireCh0Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NWireCh1Arm0Pl2++ ; //(digiIt++ ; //wire());
	    _NWireCh1Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NWireCh2Arm0Pl2++ ; //(digiIt++ ; //wire());
	    _NWireCh2Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NWireCh3Arm0Pl2++ ; //(digiIt++ ; //wire());
	    _NWireCh3Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NWireCh4Arm0Pl2++ ; //(digiIt++ ; //wire());
	    _NWireCh4Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NWireCh5Arm0Pl2++ ; //(digiIt++ ; //wire());
	    _NWireCh5Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 3 && braccio == 0){
	  switch (camera){
	  case 0:
          // _NWireCh0Arm0Pl3++ ; //(digiIt++ ; //wire());
	    _NWireCh0Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NWireCh1Arm0Pl3++ ; //(digiIt++ ; //wire());
	    _NWireCh1Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NWireCh2Arm0Pl3++ ; //(digiIt++ ; //wire());
	    _NWireCh2Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NWireCh3Arm0Pl3++ ; //(digiIt++ ; //wire());
	    _NWireCh3Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NWireCh4Arm0Pl3++ ; //(digiIt++ ; //wire());
	    _NWireCh4Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NWireCh5Arm0Pl3++ ; //(digiIt++ ; //wire());
	    _NWireCh5Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 4 && braccio == 0){
	  switch (camera){
	  case 0:
          // _NWireCh0Arm0Pl4++ ; //(digiIt++ ; //wire());
	    _NWireCh0Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NWireCh1Arm0Pl4++ ; //(digiIt++ ; //wire());
	    _NWireCh1Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NWireCh2Arm0Pl4++ ; //(digiIt++ ; //wire());
	    _NWireCh2Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NWireCh3Arm0Pl4++ ; //(digiIt++ ; //wire());
	    _NWireCh3Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NWireCh4Arm0Pl4++ ; //(digiIt++ ; //wire());
	    _NWireCh4Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NWireCh5Arm0Pl4++ ; //(digiIt++ ; //wire());
	    _NWireCh5Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

      ///////

        if(piano == 0 && braccio == 1){
	  switch (camera){
	  case 0:
          // _NWireCh0Arm1Pl0++ ; //(digiIt++ ; //wire());
	    _NWireCh0Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NWireCh1Arm1Pl0++ ; //(digiIt++ ; //wire());
	    _NWireCh1Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NWireCh2Arm1Pl0++ ; //(digiIt++ ; //wire());
	    _NWireCh2Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NWireCh3Arm1Pl0++ ; //(digiIt++ ; //wire());
	    _NWireCh3Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NWireCh4Arm1Pl0++ ; //(digiIt++ ; //wire());
	    _NWireCh4Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NWireCh5Arm1Pl0++ ; //(digiIt++ ; //wire());
	    _NWireCh5Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 1 && braccio == 1){
	  switch (camera){
	  case 0:
          // _NWireCh0Arm1Pl1++ ; //(digiIt++ ; //wire());
	    _NWireCh0Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NWireCh1Arm1Pl1++ ; //(digiIt++ ; //wire());
	    _NWireCh1Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NWireCh2Arm1Pl1++ ; //(digiIt++ ; //wire());
	    _NWireCh2Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NWireCh3Arm1Pl1++ ; //(digiIt++ ; //wire());
	    _NWireCh3Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NWireCh4Arm1Pl1++ ; //(digiIt++ ; //wire());
	    _NWireCh4Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NWireCh5Arm1Pl1++ ; //(digiIt++ ; //wire());
	    _NWireCh5Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 2 && braccio == 1){
	  switch (camera){
	  case 0:
          // _NWireCh0Arm1Pl2++ ; //(digiIt++ ; //wire());
	    _NWireCh0Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NWireCh1Arm1Pl2++ ; //(digiIt++ ; //wire());
	    _NWireCh1Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NWireCh2Arm1Pl2++ ; //(digiIt++ ; //wire());
	    _NWireCh2Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NWireCh3Arm1Pl2++ ; //(digiIt++ ; //wire());
	    _NWireCh3Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NWireCh4Arm1Pl2++ ; //(digiIt++ ; //wire());
	    _NWireCh4Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NWireCh5Arm1Pl2++ ; //(digiIt++ ; //wire());
	    _NWireCh5Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 3 && braccio == 1){
	  switch (camera){
	  case 0:
          // _NWireCh0Arm1Pl3++ ; //(digiIt++ ; //wire());
	    _NWireCh0Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NWireCh1Arm1Pl3++ ; //(digiIt++ ; //wire());
	    _NWireCh1Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NWireCh2Arm1Pl3++ ; //(digiIt++ ; //wire());
	    _NWireCh2Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NWireCh3Arm1Pl3++ ; //(digiIt++ ; //wire());
	    _NWireCh3Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NWireCh4Arm1Pl3++ ; //(digiIt++ ; //wire());
	    _NWireCh4Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NWireCh5Arm1Pl3++ ; //(digiIt++ ; //wire());
	    _NWireCh5Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 4 && braccio == 1){
	  switch (camera){
	  case 0:
          // _NWireCh0Arm1Pl4++ ; //(digiIt++ ; //wire());
	    _NWireCh0Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NWireCh1Arm1Pl4++ ; //(digiIt++ ; //wire());
	    _NWireCh1Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NWireCh2Arm1Pl4++ ; //(digiIt++ ; //wire());
	    _NWireCh2Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NWireCh3Arm1Pl4++ ; //(digiIt++ ; //wire());
	    _NWireCh3Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NWireCh4Arm1Pl4++ ; //(digiIt++ ; //wire());
	    _NWireCh4Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NWireCh5Arm1Pl4++ ; //(digiIt++ ; //wire());
	    _NWireCh5Arm1Pl4++ ; //GetXaxis()->SetTitle("N");
	    break;

	  }
        }

      
      }
    }
  


    T1DigiVfatCollection::DigiRangeIterator T1DigiS_it;
    for(T1DigiS_it=myDigiCollStrip->begin();T1DigiS_it!=myDigiCollStrip->end();T1DigiS_it++){
      int piano = (*T1DigiS_it).first.Plane();
      int braccio = (*T1DigiS_it).first.Arm();
      int camera = (*T1DigiS_it).first.CSC();
      int strato = (*T1DigiS_it).first.Layer();

      const T1DigiVfatCollection::Range& range = (*T1DigiS_it).second;
      for (T1DigiVfatCollection::const_iterator digiIt = range.first;
           digiIt!=range.second;++digiIt){
#ifdef _PRINT_
        (*digiIt).print();
#endif

       
      //       std::cout <<braccio<<piano<<camera<<strato << std::endl;
        if(piano == 0 && braccio == 0 && strato == 1){
	  switch (camera){

	  case 0:
	    laStripACh0Arm0Pl0->Fill(digiIt->strip() );
	    laStripACh0Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripACh1Arm0Pl0->Fill(digiIt->strip() );
	    laStripACh1Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripACh2Arm0Pl0->Fill(digiIt->strip() );
	    laStripACh2Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripACh3Arm0Pl0->Fill(digiIt->strip() );
	    laStripACh3Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripACh4Arm0Pl0->Fill(digiIt->strip() );
	    laStripACh4Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripACh5Arm0Pl0->Fill(digiIt->strip() );
	    laStripACh5Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;

	  }
	}

        if(piano == 1 && braccio == 0 && strato == 1){
	  switch (camera){
	  case 0:
	    laStripACh0Arm0Pl1->Fill(digiIt->strip() );
	    laStripACh0Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripACh1Arm0Pl1->Fill(digiIt->strip() );
	    laStripACh1Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripACh2Arm0Pl1->Fill(digiIt->strip() );
	    laStripACh2Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripACh3Arm0Pl1->Fill(digiIt->strip() );
	    laStripACh3Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripACh4Arm0Pl1->Fill(digiIt->strip() );
	    laStripACh4Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripACh5Arm0Pl1->Fill(digiIt->strip() );
	    laStripACh5Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 2 && braccio == 0 && strato == 1){
	  switch (camera){
	  case 0:
	    laStripACh0Arm0Pl2->Fill(digiIt->strip() );
	    laStripACh0Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripACh1Arm0Pl2->Fill(digiIt->strip() );
	    laStripACh1Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripACh2Arm0Pl2->Fill(digiIt->strip() );
	    laStripACh2Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripACh3Arm0Pl2->Fill(digiIt->strip() );
	    laStripACh3Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripACh4Arm0Pl2->Fill(digiIt->strip() );
	    laStripACh4Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripACh5Arm0Pl2->Fill(digiIt->strip() );
	    laStripACh5Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 3 && braccio == 0 && strato == 1){
	  switch (camera){
	  case 0:
	    laStripACh0Arm0Pl3->Fill(digiIt->strip() );
	    laStripACh0Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripACh1Arm0Pl3->Fill(digiIt->strip() );
	    laStripACh1Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripACh2Arm0Pl3->Fill(digiIt->strip() );
	    laStripACh2Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripACh3Arm0Pl3->Fill(digiIt->strip() );
	    laStripACh3Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripACh4Arm0Pl3->Fill(digiIt->strip() );
	    laStripACh4Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripACh5Arm0Pl3->Fill(digiIt->strip() );
	    laStripACh5Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 4 && braccio == 0 && strato == 1){
	  switch (camera){
	  case 0:
	    laStripACh0Arm0Pl4->Fill(digiIt->strip() );
	    laStripACh0Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripACh1Arm0Pl4->Fill(digiIt->strip() );
	    laStripACh1Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripACh2Arm0Pl4->Fill(digiIt->strip() );
	    laStripACh2Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripACh3Arm0Pl4->Fill(digiIt->strip() );
	    laStripACh3Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripACh4Arm0Pl4->Fill(digiIt->strip() );
	    laStripACh4Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripACh5Arm0Pl4->Fill(digiIt->strip() );
	    laStripACh5Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;

	  }
        }



      ///////



        if(piano == 0 && braccio == 1 && strato == 1){
	  switch (camera){
	  case 0:
	    laStripACh0Arm1Pl0->Fill(digiIt->strip() );
	    laStripACh0Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripACh1Arm1Pl0->Fill(digiIt->strip() );
	    laStripACh1Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripACh2Arm1Pl0->Fill(digiIt->strip() );
	    laStripACh2Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripACh3Arm1Pl0->Fill(digiIt->strip() );
	    laStripACh3Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripACh4Arm1Pl0->Fill(digiIt->strip() );
	    laStripACh4Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripACh5Arm1Pl0->Fill(digiIt->strip() );
	    laStripACh5Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 1 && braccio == 1 && strato == 1){
	  switch (camera){
	  case 0:
	    laStripACh0Arm1Pl1->Fill(digiIt->strip() );
	    laStripACh0Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripACh1Arm1Pl1->Fill(digiIt->strip() );
	    laStripACh1Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripACh2Arm1Pl1->Fill(digiIt->strip() );
	    laStripACh2Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripACh3Arm1Pl1->Fill(digiIt->strip() );
	    laStripACh3Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripACh4Arm1Pl1->Fill(digiIt->strip() );
	    laStripACh4Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripACh5Arm1Pl1->Fill(digiIt->strip() );
	    laStripACh5Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 2 && braccio == 1 && strato == 1){
	  switch (camera){
	  case 0:
	    laStripACh0Arm1Pl2->Fill(digiIt->strip() );
	    laStripACh0Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripACh1Arm1Pl2->Fill(digiIt->strip() );
	    laStripACh1Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripACh2Arm1Pl2->Fill(digiIt->strip() );
	    laStripACh2Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripACh3Arm1Pl2->Fill(digiIt->strip() );
	    laStripACh3Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripACh4Arm1Pl2->Fill(digiIt->strip() );
	    laStripACh4Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripACh5Arm1Pl2->Fill(digiIt->strip() );
	    laStripACh5Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;

	  }
	}


        if(piano == 3 && braccio == 1 && strato == 1){
	  switch (camera){
	  case 0:
	    laStripACh0Arm1Pl3->Fill(digiIt->strip() );
	    laStripACh0Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripACh1Arm1Pl3->Fill(digiIt->strip() );
	    laStripACh1Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripACh2Arm1Pl3->Fill(digiIt->strip() );
	    laStripACh2Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripACh3Arm1Pl3->Fill(digiIt->strip() );
	    laStripACh3Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripACh4Arm1Pl3->Fill(digiIt->strip() );
	    laStripACh4Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripACh5Arm1Pl3->Fill(digiIt->strip() );
	    laStripACh5Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 4 && braccio == 1 && strato == 1){
	  switch (camera){
	  case 0:
	    laStripACh0Arm1Pl4->Fill(digiIt->strip() );
	    laStripACh0Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripACh1Arm1Pl4->Fill(digiIt->strip() );
	    laStripACh1Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripACh2Arm1Pl4->Fill(digiIt->strip() );
	    laStripACh2Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripACh3Arm1Pl4->Fill(digiIt->strip() );
	    laStripACh3Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripACh4Arm1Pl4->Fill(digiIt->strip() );
	    laStripACh4Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripACh5Arm1Pl4->Fill(digiIt->strip() );
	    laStripACh5Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

//////


        if(piano == 0 && braccio == 0 && strato ==2){
	  switch (camera){
	  case 0:
	    laStripBCh0Arm0Pl0->Fill(digiIt->strip() );
	    laStripBCh0Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripBCh1Arm0Pl0->Fill(digiIt->strip() );
	    laStripBCh1Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripBCh2Arm0Pl0->Fill(digiIt->strip() );
	    laStripBCh2Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripBCh3Arm0Pl0->Fill(digiIt->strip() );
	    laStripBCh3Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripBCh4Arm0Pl0->Fill(digiIt->strip() );
	    laStripBCh4Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripBCh5Arm0Pl0->Fill(digiIt->strip() );
	    laStripBCh5Arm0Pl0->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 1 && braccio == 0 && strato ==2){
	  switch (camera){
	  case 0:
	    laStripBCh0Arm0Pl1->Fill(digiIt->strip() );
	    laStripBCh0Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripBCh1Arm0Pl1->Fill(digiIt->strip() );
	    laStripBCh1Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripBCh2Arm0Pl1->Fill(digiIt->strip() );
	    laStripBCh2Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripBCh3Arm0Pl1->Fill(digiIt->strip() );
	    laStripBCh3Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripBCh4Arm0Pl1->Fill(digiIt->strip() );
	    laStripBCh4Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripBCh5Arm0Pl1->Fill(digiIt->strip() );
	    laStripBCh5Arm0Pl1->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 2 && braccio == 0 && strato ==2){
	  switch (camera){
	  case 0:
	    laStripBCh0Arm0Pl2->Fill(digiIt->strip() );
	    laStripBCh0Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripBCh1Arm0Pl2->Fill(digiIt->strip() );
	    laStripBCh1Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripBCh2Arm0Pl2->Fill(digiIt->strip() );
	    laStripBCh2Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripBCh3Arm0Pl2->Fill(digiIt->strip() );
	    laStripBCh3Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripBCh4Arm0Pl2->Fill(digiIt->strip() );
	    laStripBCh4Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripBCh5Arm0Pl2->Fill(digiIt->strip() );
	    laStripBCh5Arm0Pl2->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 3 && braccio == 0 && strato ==2){
	  switch (camera){
	  case 0:
	    laStripBCh0Arm0Pl3->Fill(digiIt->strip() );
	    laStripBCh0Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripBCh1Arm0Pl3->Fill(digiIt->strip() );
	    laStripBCh1Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripBCh2Arm0Pl3->Fill(digiIt->strip() );
	    laStripBCh2Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripBCh3Arm0Pl3->Fill(digiIt->strip() );
	    laStripBCh3Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripBCh4Arm0Pl3->Fill(digiIt->strip() );
	    laStripBCh4Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripBCh5Arm0Pl3->Fill(digiIt->strip() );
	    laStripBCh5Arm0Pl3->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 4 && braccio == 0 && strato ==2){
	  switch (camera){
	  case 0:
	    laStripBCh0Arm0Pl4->Fill(digiIt->strip() );
	    laStripBCh0Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripBCh1Arm0Pl4->Fill(digiIt->strip() );
	    laStripBCh1Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripBCh2Arm0Pl4->Fill(digiIt->strip() );
	    laStripBCh2Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripBCh3Arm0Pl4->Fill(digiIt->strip() );
	    laStripBCh3Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripBCh4Arm0Pl4->Fill(digiIt->strip() );
	    laStripBCh4Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripBCh5Arm0Pl4->Fill(digiIt->strip() );
	    laStripBCh5Arm0Pl4->GetXaxis()->SetTitle("N");
	    break;

	  }
	}



      ///////



        if(piano == 0 && braccio == 1 && strato ==2){
	  switch (camera){
	  case 0:
	    laStripBCh0Arm1Pl0->Fill(digiIt->strip() );
	    laStripBCh0Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripBCh1Arm1Pl0->Fill(digiIt->strip() );
	    laStripBCh1Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripBCh2Arm1Pl0->Fill(digiIt->strip() );
	    laStripBCh2Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripBCh3Arm1Pl0->Fill(digiIt->strip() );
	    laStripBCh3Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripBCh4Arm1Pl0->Fill(digiIt->strip() );
	    laStripBCh4Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripBCh5Arm1Pl0->Fill(digiIt->strip() );
	    laStripBCh5Arm1Pl0->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 1 && braccio == 1 && strato ==2){
	  switch (camera){
	  case 0:
	    laStripBCh0Arm1Pl1->Fill(digiIt->strip() );
	    laStripBCh0Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripBCh1Arm1Pl1->Fill(digiIt->strip() );
	    laStripBCh1Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripBCh2Arm1Pl1->Fill(digiIt->strip() );
	    laStripBCh2Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripBCh3Arm1Pl1->Fill(digiIt->strip() );
	    laStripBCh3Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripBCh4Arm1Pl1->Fill(digiIt->strip() );
	    laStripBCh4Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripBCh5Arm1Pl1->Fill(digiIt->strip() );
	    laStripBCh5Arm1Pl1->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 2 && braccio == 1 && strato ==2){
	  switch (camera){
	  case 0:
	    laStripBCh0Arm1Pl2->Fill(digiIt->strip() );
	    laStripBCh0Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripBCh1Arm1Pl2->Fill(digiIt->strip() );
	    laStripBCh1Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripBCh2Arm1Pl2->Fill(digiIt->strip() );
	    laStripBCh2Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripBCh3Arm1Pl2->Fill(digiIt->strip() );
	    laStripBCh3Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripBCh4Arm1Pl2->Fill(digiIt->strip() );
	    laStripBCh4Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripBCh5Arm1Pl2->Fill(digiIt->strip() );
	    laStripBCh5Arm1Pl2->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 3 && braccio == 1 && strato ==2){
	  switch (camera){
	  case 0:
	    laStripBCh0Arm1Pl3->Fill(digiIt->strip() );
	    laStripBCh0Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripBCh1Arm1Pl3->Fill(digiIt->strip() );
	    laStripBCh1Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripBCh2Arm1Pl3->Fill(digiIt->strip() );
	    laStripBCh2Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripBCh3Arm1Pl3->Fill(digiIt->strip() );
	    laStripBCh3Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripBCh4Arm1Pl3->Fill(digiIt->strip() );
	    laStripBCh4Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripBCh5Arm1Pl3->Fill(digiIt->strip() );
	    laStripBCh5Arm1Pl3->GetXaxis()->SetTitle("N");
	    break;

	  }
        }

        if(piano == 4 && braccio == 1 && strato ==2){
	  switch (camera){
	  case 0:
	    laStripBCh0Arm1Pl4->Fill(digiIt->strip() );
	    laStripBCh0Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 1:
	    laStripBCh1Arm1Pl4->Fill(digiIt->strip() );
	    laStripBCh1Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 2:
	    laStripBCh2Arm1Pl4->Fill(digiIt->strip() );
	    laStripBCh2Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 3:
	    laStripBCh3Arm1Pl4->Fill(digiIt->strip() );
	    laStripBCh3Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;
	  case 4:
	    laStripBCh4Arm1Pl4->Fill(digiIt->strip() );
	    laStripBCh4Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;

	  case 5:
	    laStripBCh5Arm1Pl4->Fill(digiIt->strip() );
	    laStripBCh5Arm1Pl4->GetXaxis()->SetTitle("N");
	    break;

	  }
	}

//////////////////////////////////////////////////////////

       
      //       std::cout <<braccio<<piano<<camera<<strato << std::endl;
        if(piano == 0 && braccio == 0 && strato == 1){
	  switch (camera){

	  case 0:
          // _NStripACh0Arm0Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh0Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripACh1Arm0Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh1Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripACh2Arm0Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh2Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripACh3Arm0Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh3Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripACh4Arm0Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh4Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripACh5Arm0Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh5Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
	}

        if(piano == 1 && braccio == 0 && strato == 1){
	  switch (camera){
	  case 0:
          // _NStripACh0Arm0Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh0Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripACh1Arm0Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh1Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripACh2Arm0Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh2Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripACh3Arm0Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh3Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripACh4Arm0Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh4Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripACh5Arm0Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh5Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 2 && braccio == 0 && strato == 1){
	  switch (camera){
	  case 0:
          // _NStripACh0Arm0Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh0Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripACh1Arm0Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh1Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripACh2Arm0Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh2Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripACh3Arm0Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh3Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripACh4Arm0Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh4Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripACh5Arm0Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh5Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 3 && braccio == 0 && strato == 1){
	  switch (camera){
	  case 0:
          // _NStripACh0Arm0Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh0Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripACh1Arm0Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh1Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripACh2Arm0Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh2Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripACh3Arm0Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh3Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripACh4Arm0Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh4Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripACh5Arm0Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh5Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 4 && braccio == 0 && strato == 1){
	  switch (camera){
	  case 0:
          // _NStripACh0Arm0Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh0Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripACh1Arm0Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh1Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripACh2Arm0Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh2Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripACh3Arm0Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh3Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripACh4Arm0Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh4Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripACh5Arm0Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh5Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }



      ///////



        if(piano == 0 && braccio == 1 && strato == 1){
	  switch (camera){
	  case 0:
          // _NStripACh0Arm1Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh0Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripACh1Arm1Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh1Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripACh2Arm1Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh2Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripACh3Arm1Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh3Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripACh4Arm1Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh4Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripACh5Arm1Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh5Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 1 && braccio == 1 && strato == 1){
	  switch (camera){
	  case 0:
          // _NStripACh0Arm1Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh0Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripACh1Arm1Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh1Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripACh2Arm1Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh2Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripACh3Arm1Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh3Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripACh4Arm1Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh4Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripACh5Arm1Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh5Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 2 && braccio == 1 && strato == 1){
	  switch (camera){
	  case 0:
          // _NStripACh0Arm1Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh0Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripACh1Arm1Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh1Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripACh2Arm1Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh2Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripACh3Arm1Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh3Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripACh4Arm1Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh4Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripACh5Arm1Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh5Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
	}


        if(piano == 3 && braccio == 1 && strato == 1){
	  switch (camera){
	  case 0:
          // _NStripACh0Arm1Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh0Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripACh1Arm1Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh1Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripACh2Arm1Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh2Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripACh3Arm1Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh3Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripACh4Arm1Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh4Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripACh5Arm1Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh5Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 4 && braccio == 1 && strato == 1){
	  switch (camera){
	  case 0:
          // _NStripACh0Arm1Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh0Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripACh1Arm1Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh1Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripACh2Arm1Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh2Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripACh3Arm1Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh3Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripACh4Arm1Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh4Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripACh5Arm1Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripACh5Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

//////


        if(piano == 0 && braccio == 0 && strato ==2){
	  switch (camera){
	  case 0:
          // _NStripBCh0Arm0Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh0Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripBCh1Arm0Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh1Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripBCh2Arm0Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh2Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripBCh3Arm0Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh3Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripBCh4Arm0Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh4Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripBCh5Arm0Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh5Arm0Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 1 && braccio == 0 && strato ==2){
	  switch (camera){
	  case 0:
          // _NStripBCh0Arm0Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh0Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripBCh1Arm0Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh1Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripBCh2Arm0Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh2Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripBCh3Arm0Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh3Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripBCh4Arm0Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh4Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripBCh5Arm0Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh5Arm0Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 2 && braccio == 0 && strato ==2){
	  switch (camera){
	  case 0:
          // _NStripBCh0Arm0Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh0Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripBCh1Arm0Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh1Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripBCh2Arm0Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh2Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripBCh3Arm0Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh3Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripBCh4Arm0Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh4Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripBCh5Arm0Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh5Arm0Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 3 && braccio == 0 && strato ==2){
	  switch (camera){
	  case 0:
          // _NStripBCh0Arm0Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh0Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripBCh1Arm0Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh1Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripBCh2Arm0Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh2Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripBCh3Arm0Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh3Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripBCh4Arm0Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh4Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripBCh5Arm0Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh5Arm0Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 4 && braccio == 0 && strato ==2){
	  switch (camera){
	  case 0:
          // _NStripBCh0Arm0Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh0Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripBCh1Arm0Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh1Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripBCh2Arm0Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh2Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripBCh3Arm0Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh3Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripBCh4Arm0Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh4Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripBCh5Arm0Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh5Arm0Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
	}



      ///////



        if(piano == 0 && braccio == 1 && strato ==2){
	  switch (camera){
	  case 0:
          // _NStripBCh0Arm1Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh0Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripBCh1Arm1Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh1Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripBCh2Arm1Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh2Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripBCh3Arm1Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh3Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripBCh4Arm1Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh4Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripBCh5Arm1Pl0++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh5Arm1Pl0++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 1 && braccio == 1 && strato ==2){
	  switch (camera){
	  case 0:
          // _NStripBCh0Arm1Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh0Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripBCh1Arm1Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh1Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripBCh2Arm1Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh2Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripBCh3Arm1Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh3Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripBCh4Arm1Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh4Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripBCh5Arm1Pl1++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh5Arm1Pl1++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 2 && braccio == 1 && strato ==2){
	  switch (camera){
	  case 0:
          // _NStripBCh0Arm1Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh0Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripBCh1Arm1Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh1Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripBCh2Arm1Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh2Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripBCh3Arm1Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh3Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripBCh4Arm1Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh4Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripBCh5Arm1Pl2++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh5Arm1Pl2++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 3 && braccio == 1 && strato ==2){
	  switch (camera){
	  case 0:
          // _NStripBCh0Arm1Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh0Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripBCh1Arm1Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh1Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripBCh2Arm1Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh2Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripBCh3Arm1Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh3Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripBCh4Arm1Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh4Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripBCh5Arm1Pl3++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh5Arm1Pl3++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  }
        }

        if(piano == 4 && braccio == 1 && strato ==2){
	  switch (camera){
	  case 0:
          // _NStripBCh0Arm1Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh0Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 1:
          // _NStripBCh1Arm1Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh1Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 2:
          // _NStripBCh2Arm1Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh2Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 3:
          // _NStripBCh3Arm1Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh3Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;
	  case 4:
          // _NStripBCh4Arm1Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh4Arm1Pl4++ ; //GetXaxis()++ ; //SetTitle("N");
	    break;

	  case 5:
          // _NStripBCh5Arm1Pl4++ ; //Fill(digiIt++ ; //strip() );
	    _NStripBCh5Arm1Pl4++ ; //GetXaxis()->SetTitle("N");
	    break;

	  }
	}






      




      }
    }

    NWireCh0Arm0Pl0->Fill(_NWireCh0Arm0Pl0);
    NStripACh0Arm0Pl0->Fill(_NStripACh0Arm0Pl0);
    NStripBCh0Arm0Pl0->Fill(_NStripBCh0Arm0Pl0);
    NWireCh0Arm0Pl1->Fill(_NWireCh0Arm0Pl1);
    NStripACh0Arm0Pl1->Fill(_NStripACh0Arm0Pl1);
    NStripBCh0Arm0Pl1->Fill(_NStripBCh0Arm0Pl1);
    NWireCh0Arm0Pl2->Fill(_NWireCh0Arm0Pl2);
    NStripACh0Arm0Pl2->Fill(_NStripACh0Arm0Pl2);
    NStripBCh0Arm0Pl2->Fill(_NStripBCh0Arm0Pl2);
    NWireCh0Arm0Pl3->Fill(_NWireCh0Arm0Pl3);
    NStripACh0Arm0Pl3->Fill(_NStripACh0Arm0Pl3);
    NStripBCh0Arm0Pl3->Fill(_NStripBCh0Arm0Pl3);
    NWireCh0Arm0Pl4->Fill(_NWireCh0Arm0Pl4);
    NStripACh0Arm0Pl4->Fill(_NStripACh0Arm0Pl4);
    NStripBCh0Arm0Pl4->Fill(_NStripBCh0Arm0Pl4);
    NWireCh0Arm1Pl0->Fill(_NWireCh0Arm1Pl0);
    NStripACh0Arm1Pl0->Fill(_NStripACh0Arm1Pl0);
    NStripBCh0Arm1Pl0->Fill(_NStripBCh0Arm1Pl0);
    NWireCh0Arm1Pl1->Fill(_NWireCh0Arm1Pl1);
    NStripACh0Arm1Pl1->Fill(_NStripACh0Arm1Pl1);
    NStripBCh0Arm1Pl1->Fill(_NStripBCh0Arm1Pl1);
    NWireCh0Arm1Pl2->Fill(_NWireCh0Arm1Pl2);
    NStripACh0Arm1Pl2->Fill(_NStripACh0Arm1Pl2);
    NStripBCh0Arm1Pl2->Fill(_NStripBCh0Arm1Pl2);
    NWireCh0Arm1Pl3->Fill(_NWireCh0Arm1Pl3);
    NStripACh0Arm1Pl3->Fill(_NStripACh0Arm1Pl3);
    NStripBCh0Arm1Pl3->Fill(_NStripBCh0Arm1Pl3);
    NWireCh0Arm1Pl4->Fill(_NWireCh0Arm1Pl4);
    NStripACh0Arm1Pl4->Fill(_NStripACh0Arm1Pl4);
    NStripBCh0Arm1Pl4->Fill(_NStripBCh0Arm1Pl4);
    NWireCh1Arm0Pl0->Fill(_NWireCh1Arm0Pl0);
    NStripACh1Arm0Pl0->Fill(_NStripACh1Arm0Pl0);
    NStripBCh1Arm0Pl0->Fill(_NStripBCh1Arm0Pl0);
    NWireCh1Arm0Pl1->Fill(_NWireCh1Arm0Pl1);
    NStripACh1Arm0Pl1->Fill(_NStripACh1Arm0Pl1);
    NStripBCh1Arm0Pl1->Fill(_NStripBCh1Arm0Pl1);
    NWireCh1Arm0Pl2->Fill(_NWireCh1Arm0Pl2);
    NStripACh1Arm0Pl2->Fill(_NStripACh1Arm0Pl2);
    NStripBCh1Arm0Pl2->Fill(_NStripBCh1Arm0Pl2);
    NWireCh1Arm0Pl3->Fill(_NWireCh1Arm0Pl3);
    NStripACh1Arm0Pl3->Fill(_NStripACh1Arm0Pl3);
    NStripBCh1Arm0Pl3->Fill(_NStripBCh1Arm0Pl3);
    NWireCh1Arm0Pl4->Fill(_NWireCh1Arm0Pl4);
    NStripACh1Arm0Pl4->Fill(_NStripACh1Arm0Pl4);
    NStripBCh1Arm0Pl4->Fill(_NStripBCh1Arm0Pl4);
    NWireCh1Arm1Pl0->Fill(_NWireCh1Arm1Pl0);
    NStripACh1Arm1Pl0->Fill(_NStripACh1Arm1Pl0);
    NStripBCh1Arm1Pl0->Fill(_NStripBCh1Arm1Pl0);
    NWireCh1Arm1Pl1->Fill(_NWireCh1Arm1Pl1);
    NStripACh1Arm1Pl1->Fill(_NStripACh1Arm1Pl1);
    NStripBCh1Arm1Pl1->Fill(_NStripBCh1Arm1Pl1);
    NWireCh1Arm1Pl2->Fill(_NWireCh1Arm1Pl2);
    NStripACh1Arm1Pl2->Fill(_NStripACh1Arm1Pl2);
    NStripBCh1Arm1Pl2->Fill(_NStripBCh1Arm1Pl2);
    NWireCh1Arm1Pl3->Fill(_NWireCh1Arm1Pl3);
    NStripACh1Arm1Pl3->Fill(_NStripACh1Arm1Pl3);
    NStripBCh1Arm1Pl3->Fill(_NStripBCh1Arm1Pl3);
    NWireCh1Arm1Pl4->Fill(_NWireCh1Arm1Pl4);
    NStripACh1Arm1Pl4->Fill(_NStripACh1Arm1Pl4);
    NStripBCh1Arm1Pl4->Fill(_NStripBCh1Arm1Pl4);
    NWireCh2Arm0Pl0->Fill(_NWireCh2Arm0Pl0);
    NStripACh2Arm0Pl0->Fill(_NStripACh2Arm0Pl0);
    NStripBCh2Arm0Pl0->Fill(_NStripBCh2Arm0Pl0);
    NWireCh2Arm0Pl1->Fill(_NWireCh2Arm0Pl1);
    NStripACh2Arm0Pl1->Fill(_NStripACh2Arm0Pl1);
    NStripBCh2Arm0Pl1->Fill(_NStripBCh2Arm0Pl1);
    NWireCh2Arm0Pl2->Fill(_NWireCh2Arm0Pl2);
    NStripACh2Arm0Pl2->Fill(_NStripACh2Arm0Pl2);
    NStripBCh2Arm0Pl2->Fill(_NStripBCh2Arm0Pl2);
    NWireCh2Arm0Pl3->Fill(_NWireCh2Arm0Pl3);
    NStripACh2Arm0Pl3->Fill(_NStripACh2Arm0Pl3);
    NStripBCh2Arm0Pl3->Fill(_NStripBCh2Arm0Pl3);
    NWireCh2Arm0Pl4->Fill(_NWireCh2Arm0Pl4);
    NStripACh2Arm0Pl4->Fill(_NStripACh2Arm0Pl4);
    NStripBCh2Arm0Pl4->Fill(_NStripBCh2Arm0Pl4);
    NWireCh2Arm1Pl0->Fill(_NWireCh2Arm1Pl0);
    NStripACh2Arm1Pl0->Fill(_NStripACh2Arm1Pl0);
    NStripBCh2Arm1Pl0->Fill(_NStripBCh2Arm1Pl0);
    NWireCh2Arm1Pl1->Fill(_NWireCh2Arm1Pl1);
    NStripACh2Arm1Pl1->Fill(_NStripACh2Arm1Pl1);
    NStripBCh2Arm1Pl1->Fill(_NStripBCh2Arm1Pl1);
    NWireCh2Arm1Pl2->Fill(_NWireCh2Arm1Pl2);
    NStripACh2Arm1Pl2->Fill(_NStripACh2Arm1Pl2);
    NStripBCh2Arm1Pl2->Fill(_NStripBCh2Arm1Pl2);
    NWireCh2Arm1Pl3->Fill(_NWireCh2Arm1Pl3);
    NStripACh2Arm1Pl3->Fill(_NStripACh2Arm1Pl3);
    NStripBCh2Arm1Pl3->Fill(_NStripBCh2Arm1Pl3);
    NWireCh2Arm1Pl4->Fill(_NWireCh2Arm1Pl4);
    NStripACh2Arm1Pl4->Fill(_NStripACh2Arm1Pl4);
    NStripBCh2Arm1Pl4->Fill(_NStripBCh2Arm1Pl4);
    NWireCh3Arm0Pl0->Fill(_NWireCh3Arm0Pl0);
    NStripACh3Arm0Pl0->Fill(_NStripACh3Arm0Pl0);
    NStripBCh3Arm0Pl0->Fill(_NStripBCh3Arm0Pl0);
    NWireCh3Arm0Pl1->Fill(_NWireCh3Arm0Pl1);
    NStripACh3Arm0Pl1->Fill(_NStripACh3Arm0Pl1);
    NStripBCh3Arm0Pl1->Fill(_NStripBCh3Arm0Pl1);
    NWireCh3Arm0Pl2->Fill(_NWireCh3Arm0Pl2);
    NStripACh3Arm0Pl2->Fill(_NStripACh3Arm0Pl2);
    NStripBCh3Arm0Pl2->Fill(_NStripBCh3Arm0Pl2);
    NWireCh3Arm0Pl3->Fill(_NWireCh3Arm0Pl3);
    NStripACh3Arm0Pl3->Fill(_NStripACh3Arm0Pl3);
    NStripBCh3Arm0Pl3->Fill(_NStripBCh3Arm0Pl3);
    NWireCh3Arm0Pl4->Fill(_NWireCh3Arm0Pl4);
    NStripACh3Arm0Pl4->Fill(_NStripACh3Arm0Pl4);
    NStripBCh3Arm0Pl4->Fill(_NStripBCh3Arm0Pl4);
    NWireCh3Arm1Pl0->Fill(_NWireCh3Arm1Pl0);
    NStripACh3Arm1Pl0->Fill(_NStripACh3Arm1Pl0);
    NStripBCh3Arm1Pl0->Fill(_NStripBCh3Arm1Pl0);
    NWireCh3Arm1Pl1->Fill(_NWireCh3Arm1Pl1);
    NStripACh3Arm1Pl1->Fill(_NStripACh3Arm1Pl1);
    NStripBCh3Arm1Pl1->Fill(_NStripBCh3Arm1Pl1);
    NWireCh3Arm1Pl2->Fill(_NWireCh3Arm1Pl2);
    NStripACh3Arm1Pl2->Fill(_NStripACh3Arm1Pl2);
    NStripBCh3Arm1Pl2->Fill(_NStripBCh3Arm1Pl2);
    NWireCh3Arm1Pl3->Fill(_NWireCh3Arm1Pl3);
    NStripACh3Arm1Pl3->Fill(_NStripACh3Arm1Pl3);
    NStripBCh3Arm1Pl3->Fill(_NStripBCh3Arm1Pl3);
    NWireCh3Arm1Pl4->Fill(_NWireCh3Arm1Pl4);
    NStripACh3Arm1Pl4->Fill(_NStripACh3Arm1Pl4);
    NStripBCh3Arm1Pl4->Fill(_NStripBCh3Arm1Pl4);
    NWireCh4Arm0Pl0->Fill(_NWireCh4Arm0Pl0);
    NStripACh4Arm0Pl0->Fill(_NStripACh4Arm0Pl0);
    NStripBCh4Arm0Pl0->Fill(_NStripBCh4Arm0Pl0);
    NWireCh4Arm0Pl1->Fill(_NWireCh4Arm0Pl1);
    NStripACh4Arm0Pl1->Fill(_NStripACh4Arm0Pl1);
    NStripBCh4Arm0Pl1->Fill(_NStripBCh4Arm0Pl1);
    NWireCh4Arm0Pl2->Fill(_NWireCh4Arm0Pl2);
    NStripACh4Arm0Pl2->Fill(_NStripACh4Arm0Pl2);
    NStripBCh4Arm0Pl2->Fill(_NStripBCh4Arm0Pl2);
    NWireCh4Arm0Pl3->Fill(_NWireCh4Arm0Pl3);
    NStripACh4Arm0Pl3->Fill(_NStripACh4Arm0Pl3);
    NStripBCh4Arm0Pl3->Fill(_NStripBCh4Arm0Pl3);
    NWireCh4Arm0Pl4->Fill(_NWireCh4Arm0Pl4);
    NStripACh4Arm0Pl4->Fill(_NStripACh4Arm0Pl4);
    NStripBCh4Arm0Pl4->Fill(_NStripBCh4Arm0Pl4);
    NWireCh4Arm1Pl0->Fill(_NWireCh4Arm1Pl0);
    NStripACh4Arm1Pl0->Fill(_NStripACh4Arm1Pl0);
    NStripBCh4Arm1Pl0->Fill(_NStripBCh4Arm1Pl0);
    NWireCh4Arm1Pl1->Fill(_NWireCh4Arm1Pl1);
    NStripACh4Arm1Pl1->Fill(_NStripACh4Arm1Pl1);
    NStripBCh4Arm1Pl1->Fill(_NStripBCh4Arm1Pl1);
    NWireCh4Arm1Pl2->Fill(_NWireCh4Arm1Pl2);
    NStripACh4Arm1Pl2->Fill(_NStripACh4Arm1Pl2);
    NStripBCh4Arm1Pl2->Fill(_NStripBCh4Arm1Pl2);
    NWireCh4Arm1Pl3->Fill(_NWireCh4Arm1Pl3);
    NStripACh4Arm1Pl3->Fill(_NStripACh4Arm1Pl3);
    NStripBCh4Arm1Pl3->Fill(_NStripBCh4Arm1Pl3);
    NWireCh4Arm1Pl4->Fill(_NWireCh4Arm1Pl4);
    NStripACh4Arm1Pl4->Fill(_NStripACh4Arm1Pl4);
    NStripBCh4Arm1Pl4->Fill(_NStripBCh4Arm1Pl4);
    NWireCh5Arm0Pl0->Fill(_NWireCh5Arm0Pl0);
    NStripACh5Arm0Pl0->Fill(_NStripACh5Arm0Pl0);
    NStripBCh5Arm0Pl0->Fill(_NStripBCh5Arm0Pl0);
    NWireCh5Arm0Pl1->Fill(_NWireCh5Arm0Pl1);
    NStripACh5Arm0Pl1->Fill(_NStripACh5Arm0Pl1);
    NStripBCh5Arm0Pl1->Fill(_NStripBCh5Arm0Pl1);
    NWireCh5Arm0Pl2->Fill(_NWireCh5Arm0Pl2);
    NStripACh5Arm0Pl2->Fill(_NStripACh5Arm0Pl2);
    NStripBCh5Arm0Pl2->Fill(_NStripBCh5Arm0Pl2);
    NWireCh5Arm0Pl3->Fill(_NWireCh5Arm0Pl3);
    NStripACh5Arm0Pl3->Fill(_NStripACh5Arm0Pl3);
    NStripBCh5Arm0Pl3->Fill(_NStripBCh5Arm0Pl3);
    NWireCh5Arm0Pl4->Fill(_NWireCh5Arm0Pl4);
    NStripACh5Arm0Pl4->Fill(_NStripACh5Arm0Pl4);
    NStripBCh5Arm0Pl4->Fill(_NStripBCh5Arm0Pl4);
    NWireCh5Arm1Pl0->Fill(_NWireCh5Arm1Pl0);
    NStripACh5Arm1Pl0->Fill(_NStripACh5Arm1Pl0);
    NStripBCh5Arm1Pl0->Fill(_NStripBCh5Arm1Pl0);
    NWireCh5Arm1Pl1->Fill(_NWireCh5Arm1Pl1);
    NStripACh5Arm1Pl1->Fill(_NStripACh5Arm1Pl1);
    NStripBCh5Arm1Pl1->Fill(_NStripBCh5Arm1Pl1);
    NWireCh5Arm1Pl2->Fill(_NWireCh5Arm1Pl2);
    NStripACh5Arm1Pl2->Fill(_NStripACh5Arm1Pl2);
    NStripBCh5Arm1Pl2->Fill(_NStripBCh5Arm1Pl2);
    NWireCh5Arm1Pl3->Fill(_NWireCh5Arm1Pl3);
    NStripACh5Arm1Pl3->Fill(_NStripACh5Arm1Pl3);
    NStripBCh5Arm1Pl3->Fill(_NStripBCh5Arm1Pl3);
    NWireCh5Arm1Pl4->Fill(_NWireCh5Arm1Pl4);
    NStripACh5Arm1Pl4->Fill(_NStripACh5Arm1Pl4);
    NStripBCh5Arm1Pl4->Fill(_NStripBCh5Arm1Pl4);


  }
  
  
#endif

  

}




// ------------ method called once each job just before starting event loop  ------------
void
T1Validation::beginJob()
{
  initHistograms();
}

// ------------ method called once each job just after ending the event loop  ------------
void
T1Validation::endJob() {
  theFile = TFile::Open(outputFileName.c_str(), "recreate");
  if(!theFile || !theFile->IsWritable())
    {
      std::cout<<"Output file not opened correctly!!"<<std::endl;
    }
  writeHistograms();
  theFile->Close();
}

void
T1Validation::initHistograms() {

//  hClusterWidth=std::auto_ptr<TH1D>(new TH1D("Cluster Width","Cluster Width",200,0,200));
//  hClusterWidth->SetDirectory(0);
  RecoXY_00 = std::auto_ptr<TH2D>(new TH2D("XY_00","XY_00",500,-1200,1200,500,-1200,1200));
  RecoXY_01 = std::auto_ptr<TH2D>(new TH2D("XY_01","XY_01",500,-1200,1200,500,-1200,1200));
  RecoXY_02 = std::auto_ptr<TH2D>(new TH2D("XY_02","XY_02",500,-1200,1200,500,-1200,1200));
  RecoXY_03 = std::auto_ptr<TH2D>(new TH2D("XY_03","XY_03",500,-1200,1200,500,-1200,1200));
  RecoXY_04 = std::auto_ptr<TH2D>(new TH2D("XY_04","XY_04",500,-1200,1200,500,-1200,1200));
  RecoXY_10 = std::auto_ptr<TH2D>(new TH2D("XY_10","XY_10",500,-1200,1200,500,-1200,1200));
  RecoXY_11 = std::auto_ptr<TH2D>(new TH2D("XY_11","XY_11",500,-1200,1200,500,-1200,1200));
  RecoXY_12 = std::auto_ptr<TH2D>(new TH2D("XY_12","XY_12",500,-1200,1200,500,-1200,1200));
  RecoXY_13 = std::auto_ptr<TH2D>(new TH2D("XY_13","XY_13",500,-1200,1200,500,-1200,1200));
  RecoXY_14 = std::auto_ptr<TH2D>(new TH2D("XY_14","XY_14",500,-1200,1200,500,-1200,1200));
  RecoXY_00 -> SetDirectory(0); 
  RecoXY_01 -> SetDirectory(0); 
  RecoXY_02 -> SetDirectory(0); 
  RecoXY_03 -> SetDirectory(0); 
  RecoXY_04 -> SetDirectory(0); 
  RecoXY_10 -> SetDirectory(0); 
  RecoXY_11 -> SetDirectory(0); 
  RecoXY_12 -> SetDirectory(0); 
  RecoXY_13 -> SetDirectory(0); 
  RecoXY_14 -> SetDirectory(0); 

  lEta = std::auto_ptr<TH1D>(new TH1D("SimEta","SimEta",55,-8,8));
  lEta->SetDirectory(0);
  ilPT = std::auto_ptr<TH1D>(new TH1D("SimPT","SimPT",500,0,100));
  ilPT->SetDirectory(0);

  hClDiffBetweenLayers = std::auto_ptr<TH1D>(new TH1D("ClDiffBetweenLayers","ClDiffBetweenLayers",200,-10,10));
  hClDiffBetweenLayers ->SetDirectory(0);
  hClRMSSex = std::auto_ptr<TH1D>(new TH1D("ClRMS","ClRMS",200,0,10));
  hClRMSSex->SetDirectory(0);

  ilWireCh0Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("WireCh0Arm0Pl0","WireCh0 Arm0Pl0",250,0.5,250.5));
  ilWireCh0Arm0Pl0->SetDirectory(0);
  ilWireCh0Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("WireCh0Arm0Pl1","WireCh0 Arm0Pl1",250,0.5,250.5));
  ilWireCh0Arm0Pl1->SetDirectory(0);
  ilWireCh0Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("WireCh0Arm0Pl2","WireCh0 Arm0Pl2",250,0.5,250.5));
  ilWireCh0Arm0Pl2->SetDirectory(0);
  ilWireCh0Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("WireCh0Arm0Pl3","WireCh0 Arm0Pl3",250,0.5,250.5));
  ilWireCh0Arm0Pl3->SetDirectory(0);
  ilWireCh0Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("WireCh0Arm0Pl4","WireCh0 Arm0Pl4",250,0.5,250.5));
  ilWireCh0Arm0Pl4->SetDirectory(0);

  ilWireCh0Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("WireCh0Arm1Pl0","WireCh0 Arm1Pl0",250,0.5,250.5));
  ilWireCh0Arm1Pl0->SetDirectory(0);
  ilWireCh0Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("WireCh0Arm1Pl1","WireCh0 Arm1Pl1",250,0.5,250.5));
  ilWireCh0Arm1Pl1->SetDirectory(0);
  ilWireCh0Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("WireCh0Arm1Pl2","WireCh0 Arm1Pl2",250,0.5,250.5));
  ilWireCh0Arm1Pl2->SetDirectory(0);
  ilWireCh0Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("WireCh0Arm1Pl3","WireCh0 Arm1Pl3",250,0.5,250.5));
  ilWireCh0Arm1Pl3->SetDirectory(0);
  ilWireCh0Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("WireCh0Arm1Pl4","WireCh0 Arm1Pl4",250,0.5,250.5));
  ilWireCh0Arm1Pl4->SetDirectory(0);


  ilClWidthCh0Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh0Arm0Pl0","ClWidthCh0 Arm0Pl0",250,0.5,250.5));
  ilClWidthCh0Arm0Pl0->SetDirectory(0);
  ilClWidthCh0Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh0Arm0Pl1","ClWidthCh0 Arm0Pl1",250,0.5,250.5));
  ilClWidthCh0Arm0Pl1->SetDirectory(0);
  ilClWidthCh0Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh0Arm0Pl2","ClWidthCh0 Arm0Pl2",250,0.5,250.5));
  ilClWidthCh0Arm0Pl2->SetDirectory(0);
  ilClWidthCh0Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh0Arm0Pl3","ClWidthCh0 Arm0Pl3",250,0.5,250.5));
  ilClWidthCh0Arm0Pl3->SetDirectory(0);
  ilClWidthCh0Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh0Arm0Pl4","ClWidthCh0 Arm0Pl4",250,0.5,250.5));
  ilClWidthCh0Arm0Pl4->SetDirectory(0);

  ilClWidthCh0Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh0Arm1Pl0","ClWidthCh0 Arm1Pl0",250,0.5,250.5));
  ilClWidthCh0Arm1Pl0->SetDirectory(0);
  ilClWidthCh0Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh0Arm1Pl1","ClWidthCh0 Arm1Pl1",250,0.5,250.5));
  ilClWidthCh0Arm1Pl1->SetDirectory(0);
  ilClWidthCh0Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh0Arm1Pl2","ClWidthCh0 Arm1Pl2",250,0.5,250.5));
  ilClWidthCh0Arm1Pl2->SetDirectory(0);
  ilClWidthCh0Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh0Arm1Pl3","ClWidthCh0 Arm1Pl3",250,0.5,250.5));
  ilClWidthCh0Arm1Pl3->SetDirectory(0);
  ilClWidthCh0Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh0Arm1Pl4","ClWidthCh0 Arm1Pl4",250,0.5,250.5));
  ilClWidthCh0Arm1Pl4->SetDirectory(0);




  laStripACh0Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("StripACh0Arm0Pl0","StripACh0 Arm0Pl0",250,0.5,250.5));
  laStripACh0Arm0Pl0->SetDirectory(0);
  laStripACh0Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("StripACh0Arm0Pl1","StripACh0 Arm0Pl1",250,0.5,250.5));
  laStripACh0Arm0Pl1->SetDirectory(0);
  laStripACh0Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("StripACh0Arm0Pl2","StripACh0 Arm0Pl2",250,0.5,250.5));
  laStripACh0Arm0Pl2->SetDirectory(0);
  laStripACh0Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("StripACh0Arm0Pl3","StripACh0 Arm0Pl3",250,0.5,250.5));
  laStripACh0Arm0Pl3->SetDirectory(0);
  laStripACh0Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("StripACh0Arm0Pl4","StripACh0 Arm0Pl4",250,0.5,250.5));
  laStripACh0Arm0Pl4->SetDirectory(0);

  laStripACh0Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("StripACh0Arm1Pl0","StripACh0 Arm1Pl0",250,0.5,250.5));
  laStripACh0Arm1Pl0->SetDirectory(0);
  laStripACh0Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("StripACh0Arm1Pl1","StripACh0 Arm1Pl1",250,0.5,250.5));
  laStripACh0Arm1Pl1->SetDirectory(0);
  laStripACh0Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("StripACh0Arm1Pl2","StripACh0 Arm1Pl2",250,0.5,250.5));
  laStripACh0Arm1Pl2->SetDirectory(0);
  laStripACh0Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("StripACh0Arm1Pl3","StripACh0 Arm1Pl3",250,0.5,250.5));
  laStripACh0Arm1Pl3->SetDirectory(0);
  laStripACh0Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("StripACh0Arm1Pl4","StripACh0 Arm1Pl4",250,0.5,250.5));
  laStripACh0Arm1Pl4->SetDirectory(0);


  laStripBCh0Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("StripBCh0Arm0Pl0","StripBCh0 Arm0Pl0",250,0.5,250.5));
  laStripBCh0Arm0Pl0->SetDirectory(0);
  laStripBCh0Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("StripBCh0Arm0Pl1","StripBCh0 Arm0Pl1",250,0.5,250.5));
  laStripBCh0Arm0Pl1->SetDirectory(0);
  laStripBCh0Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("StripBCh0Arm0Pl2","StripBCh0 Arm0Pl2",250,0.5,250.5));
  laStripBCh0Arm0Pl2->SetDirectory(0);
  laStripBCh0Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("StripBCh0Arm0Pl3","StripBCh0 Arm0Pl3",250,0.5,250.5));
  laStripBCh0Arm0Pl3->SetDirectory(0);
  laStripBCh0Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("StripBCh0Arm0Pl4","StripBCh0 Arm0Pl4",250,0.5,250.5));
  laStripBCh0Arm0Pl4->SetDirectory(0);

  laStripBCh0Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("StripBCh0Arm1Pl0","StripBCh0 Arm1Pl0",250,0.5,250.5));
  laStripBCh0Arm1Pl0->SetDirectory(0);
  laStripBCh0Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("StripBCh0Arm1Pl1","StripBCh0 Arm1Pl1",250,0.5,250.5));
  laStripBCh0Arm1Pl1->SetDirectory(0);
  laStripBCh0Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("StripBCh0Arm1Pl2","StripBCh0 Arm1Pl2",250,0.5,250.5));
  laStripBCh0Arm1Pl2->SetDirectory(0);
  laStripBCh0Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("StripBCh0Arm1Pl3","StripBCh0 Arm1Pl3",250,0.5,250.5));
  laStripBCh0Arm1Pl3->SetDirectory(0);
  laStripBCh0Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("StripBCh0Arm1Pl4","StripBCh0 Arm1Pl4",250,0.5,250.5));
  laStripBCh0Arm1Pl4->SetDirectory(0);


  ilWireCh1Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("WireCh1Arm0Pl0","WireCh1 Arm0Pl0",250,0.5,250.5));
  ilWireCh1Arm0Pl0->SetDirectory(0);
  ilWireCh1Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("WireCh1Arm0Pl1","WireCh1 Arm0Pl1",250,0.5,250.5));
  ilWireCh1Arm0Pl1->SetDirectory(0);
  ilWireCh1Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("WireCh1Arm0Pl2","WireCh1 Arm0Pl2",250,0.5,250.5));
  ilWireCh1Arm0Pl2->SetDirectory(0);
  ilWireCh1Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("WireCh1Arm0Pl3","WireCh1 Arm0Pl3",250,0.5,250.5));
  ilWireCh1Arm0Pl3->SetDirectory(0);
  ilWireCh1Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("WireCh1Arm0Pl4","WireCh1 Arm0Pl4",250,0.5,250.5));
  ilWireCh1Arm0Pl4->SetDirectory(0);

  ilWireCh1Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("WireCh1Arm1Pl0","WireCh1 Arm1Pl0",250,0.5,250.5));
  ilWireCh1Arm1Pl0->SetDirectory(0);
  ilWireCh1Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("WireCh1Arm1Pl1","WireCh1 Arm1Pl1",250,0.5,250.5));
  ilWireCh1Arm1Pl1->SetDirectory(0);
  ilWireCh1Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("WireCh1Arm1Pl2","WireCh1 Arm1Pl2",250,0.5,250.5));
  ilWireCh1Arm1Pl2->SetDirectory(0);
  ilWireCh1Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("WireCh1Arm1Pl3","WireCh1 Arm1Pl3",250,0.5,250.5));
  ilWireCh1Arm1Pl3->SetDirectory(0);
  ilWireCh1Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("WireCh1Arm1Pl4","WireCh1 Arm1Pl4",250,0.5,250.5));
  ilWireCh1Arm1Pl4->SetDirectory(0);




  ilClWidthCh1Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh1Arm0Pl0","ClWidthCh1 Arm0Pl0",250,0.5,250.5));
  ilClWidthCh1Arm0Pl0->SetDirectory(0);
  ilClWidthCh1Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh1Arm0Pl1","ClWidthCh1 Arm0Pl1",250,0.5,250.5));
  ilClWidthCh1Arm0Pl1->SetDirectory(0);
  ilClWidthCh1Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh1Arm0Pl2","ClWidthCh1 Arm0Pl2",250,0.5,250.5));
  ilClWidthCh1Arm0Pl2->SetDirectory(0);
  ilClWidthCh1Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh1Arm0Pl3","ClWidthCh1 Arm0Pl3",250,0.5,250.5));
  ilClWidthCh1Arm0Pl3->SetDirectory(0);
  ilClWidthCh1Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh1Arm0Pl4","ClWidthCh1 Arm0Pl4",250,0.5,250.5));
  ilClWidthCh1Arm0Pl4->SetDirectory(0);

  ilClWidthCh1Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh1Arm1Pl0","ClWidthCh1 Arm1Pl0",250,0.5,250.5));
  ilClWidthCh1Arm1Pl0->SetDirectory(0);
  ilClWidthCh1Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh1Arm1Pl1","ClWidthCh1 Arm1Pl1",250,0.5,250.5));
  ilClWidthCh1Arm1Pl1->SetDirectory(0);
  ilClWidthCh1Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh1Arm1Pl2","ClWidthCh1 Arm1Pl2",250,0.5,250.5));
  ilClWidthCh1Arm1Pl2->SetDirectory(0);
  ilClWidthCh1Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh1Arm1Pl3","ClWidthCh1 Arm1Pl3",250,0.5,250.5));
  ilClWidthCh1Arm1Pl3->SetDirectory(0);
  ilClWidthCh1Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh1Arm1Pl4","ClWidthCh1 Arm1Pl4",250,0.5,250.5));
  ilClWidthCh1Arm1Pl4->SetDirectory(0);





  laStripACh1Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("StripACh1Arm0Pl0","StripACh1 Arm0Pl0",250,0.5,250.5));
  laStripACh1Arm0Pl0->SetDirectory(0);
  laStripACh1Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("StripACh1Arm0Pl1","StripACh1 Arm0Pl1",250,0.5,250.5));
  laStripACh1Arm0Pl1->SetDirectory(0);
  laStripACh1Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("StripACh1Arm0Pl2","StripACh1 Arm0Pl2",250,0.5,250.5));
  laStripACh1Arm0Pl2->SetDirectory(0);
  laStripACh1Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("StripACh1Arm0Pl3","StripACh1 Arm0Pl3",250,0.5,250.5));
  laStripACh1Arm0Pl3->SetDirectory(0);
  laStripACh1Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("StripACh1Arm0Pl4","StripACh1 Arm0Pl4",250,0.5,250.5));
  laStripACh1Arm0Pl4->SetDirectory(0);

  laStripACh1Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("StripACh1Arm1Pl0","StripACh1 Arm1Pl0",250,0.5,250.5));
  laStripACh1Arm1Pl0->SetDirectory(0);
  laStripACh1Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("StripACh1Arm1Pl1","StripACh1 Arm1Pl1",250,0.5,250.5));
  laStripACh1Arm1Pl1->SetDirectory(0);
  laStripACh1Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("StripACh1Arm1Pl2","StripACh1 Arm1Pl2",250,0.5,250.5));
  laStripACh1Arm1Pl2->SetDirectory(0);
  laStripACh1Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("StripACh1Arm1Pl3","StripACh1 Arm1Pl3",250,0.5,250.5));
  laStripACh1Arm1Pl3->SetDirectory(0);
  laStripACh1Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("StripACh1Arm1Pl4","StripACh1 Arm1Pl4",250,0.5,250.5));
  laStripACh1Arm1Pl4->SetDirectory(0);


  laStripBCh1Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("StripBCh1Arm0Pl0","StripBCh1 Arm0Pl0",250,0.5,250.5));
  laStripBCh1Arm0Pl0->SetDirectory(0);
  laStripBCh1Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("StripBCh1Arm0Pl1","StripBCh1 Arm0Pl1",250,0.5,250.5));
  laStripBCh1Arm0Pl1->SetDirectory(0);
  laStripBCh1Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("StripBCh1Arm0Pl2","StripBCh1 Arm0Pl2",250,0.5,250.5));
  laStripBCh1Arm0Pl2->SetDirectory(0);
  laStripBCh1Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("StripBCh1Arm0Pl3","StripBCh1 Arm0Pl3",250,0.5,250.5));
  laStripBCh1Arm0Pl3->SetDirectory(0);
  laStripBCh1Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("StripBCh1Arm0Pl4","StripBCh1 Arm0Pl4",250,0.5,250.5));
  laStripBCh1Arm0Pl4->SetDirectory(0);

  laStripBCh1Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("StripBCh1Arm1Pl0","StripBCh1 Arm1Pl0",250,0.5,250.5));
  laStripBCh1Arm1Pl0->SetDirectory(0);
  laStripBCh1Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("StripBCh1Arm1Pl1","StripBCh1 Arm1Pl1",250,0.5,250.5));
  laStripBCh1Arm1Pl1->SetDirectory(0);
  laStripBCh1Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("StripBCh1Arm1Pl2","StripBCh1 Arm1Pl2",250,0.5,250.5));
  laStripBCh1Arm1Pl2->SetDirectory(0);
  laStripBCh1Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("StripBCh1Arm1Pl3","StripBCh1 Arm1Pl3",250,0.5,250.5));
  laStripBCh1Arm1Pl3->SetDirectory(0);
  laStripBCh1Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("StripBCh1Arm1Pl4","StripBCh1 Arm1Pl4",250,0.5,250.5));
  laStripBCh1Arm1Pl4->SetDirectory(0);


  ilWireCh2Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("WireCh2Arm0Pl0","WireCh2 Arm0Pl0",250,0.5,250.5));
  ilWireCh2Arm0Pl0->SetDirectory(0);
  ilWireCh2Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("WireCh2Arm0Pl1","WireCh2 Arm0Pl1",250,0.5,250.5));
  ilWireCh2Arm0Pl1->SetDirectory(0);
  ilWireCh2Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("WireCh2Arm0Pl2","WireCh2 Arm0Pl2",250,0.5,250.5));
  ilWireCh2Arm0Pl2->SetDirectory(0);
  ilWireCh2Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("WireCh2Arm0Pl3","WireCh2 Arm0Pl3",250,0.5,250.5));
  ilWireCh2Arm0Pl3->SetDirectory(0);
  ilWireCh2Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("WireCh2Arm0Pl4","WireCh2 Arm0Pl4",250,0.5,250.5));
  ilWireCh2Arm0Pl4->SetDirectory(0);

  ilWireCh2Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("WireCh2Arm1Pl0","WireCh2 Arm1Pl0",250,0.5,250.5));
  ilWireCh2Arm1Pl0->SetDirectory(0);
  ilWireCh2Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("WireCh2Arm1Pl1","WireCh2 Arm1Pl1",250,0.5,250.5));
  ilWireCh2Arm1Pl1->SetDirectory(0);
  ilWireCh2Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("WireCh2Arm1Pl2","WireCh2 Arm1Pl2",250,0.5,250.5));
  ilWireCh2Arm1Pl2->SetDirectory(0);
  ilWireCh2Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("WireCh2Arm1Pl3","WireCh2 Arm1Pl3",250,0.5,250.5));
  ilWireCh2Arm1Pl3->SetDirectory(0);
  ilWireCh2Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("WireCh2Arm1Pl4","WireCh2 Arm1Pl4",250,0.5,250.5));
  ilWireCh2Arm1Pl4->SetDirectory(0);






  ilClWidthCh2Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh2Arm0Pl0","ClWidthCh2 Arm0Pl0",250,0.5,250.5));
  ilClWidthCh2Arm0Pl0->SetDirectory(0);
  ilClWidthCh2Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh2Arm0Pl1","ClWidthCh2 Arm0Pl1",250,0.5,250.5));
  ilClWidthCh2Arm0Pl1->SetDirectory(0);
  ilClWidthCh2Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh2Arm0Pl2","ClWidthCh2 Arm0Pl2",250,0.5,250.5));
  ilClWidthCh2Arm0Pl2->SetDirectory(0);
  ilClWidthCh2Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh2Arm0Pl3","ClWidthCh2 Arm0Pl3",250,0.5,250.5));
  ilClWidthCh2Arm0Pl3->SetDirectory(0);
  ilClWidthCh2Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh2Arm0Pl4","ClWidthCh2 Arm0Pl4",250,0.5,250.5));
  ilClWidthCh2Arm0Pl4->SetDirectory(0);

  ilClWidthCh2Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh2Arm1Pl0","ClWidthCh2 Arm1Pl0",250,0.5,250.5));
  ilClWidthCh2Arm1Pl0->SetDirectory(0);
  ilClWidthCh2Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh2Arm1Pl1","ClWidthCh2 Arm1Pl1",250,0.5,250.5));
  ilClWidthCh2Arm1Pl1->SetDirectory(0);
  ilClWidthCh2Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh2Arm1Pl2","ClWidthCh2 Arm1Pl2",250,0.5,250.5));
  ilClWidthCh2Arm1Pl2->SetDirectory(0);
  ilClWidthCh2Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh2Arm1Pl3","ClWidthCh2 Arm1Pl3",250,0.5,250.5));
  ilClWidthCh2Arm1Pl3->SetDirectory(0);
  ilClWidthCh2Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh2Arm1Pl4","ClWidthCh2 Arm1Pl4",250,0.5,250.5));
  ilClWidthCh2Arm1Pl4->SetDirectory(0);





  laStripACh2Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("StripACh2Arm0Pl0","StripACh2 Arm0Pl0",250,0.5,250.5));
  laStripACh2Arm0Pl0->SetDirectory(0);
  laStripACh2Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("StripACh2Arm0Pl1","StripACh2 Arm0Pl1",250,0.5,250.5));
  laStripACh2Arm0Pl1->SetDirectory(0);
  laStripACh2Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("StripACh2Arm0Pl2","StripACh2 Arm0Pl2",250,0.5,250.5));
  laStripACh2Arm0Pl2->SetDirectory(0);
  laStripACh2Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("StripACh2Arm0Pl3","StripACh2 Arm0Pl3",250,0.5,250.5));
  laStripACh2Arm0Pl3->SetDirectory(0);
  laStripACh2Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("StripACh2Arm0Pl4","StripACh2 Arm0Pl4",250,0.5,250.5));
  laStripACh2Arm0Pl4->SetDirectory(0);

  laStripACh2Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("StripACh2Arm1Pl0","StripACh2 Arm1Pl0",250,0.5,250.5));
  laStripACh2Arm1Pl0->SetDirectory(0);
  laStripACh2Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("StripACh2Arm1Pl1","StripACh2 Arm1Pl1",250,0.5,250.5));
  laStripACh2Arm1Pl1->SetDirectory(0);
  laStripACh2Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("StripACh2Arm1Pl2","StripACh2 Arm1Pl2",250,0.5,250.5));
  laStripACh2Arm1Pl2->SetDirectory(0);
  laStripACh2Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("StripACh2Arm1Pl3","StripACh2 Arm1Pl3",250,0.5,250.5));
  laStripACh2Arm1Pl3->SetDirectory(0);
  laStripACh2Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("StripACh2Arm1Pl4","StripACh2 Arm1Pl4",250,0.5,250.5));
  laStripACh2Arm1Pl4->SetDirectory(0);


  laStripBCh2Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("StripBCh2Arm0Pl0","StripBCh2 Arm0Pl0",250,0.5,250.5));
  laStripBCh2Arm0Pl0->SetDirectory(0);
  laStripBCh2Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("StripBCh2Arm0Pl1","StripBCh2 Arm0Pl1",250,0.5,250.5));
  laStripBCh2Arm0Pl1->SetDirectory(0);
  laStripBCh2Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("StripBCh2Arm0Pl2","StripBCh2 Arm0Pl2",250,0.5,250.5));
  laStripBCh2Arm0Pl2->SetDirectory(0);
  laStripBCh2Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("StripBCh2Arm0Pl3","StripBCh2 Arm0Pl3",250,0.5,250.5));
  laStripBCh2Arm0Pl3->SetDirectory(0);
  laStripBCh2Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("StripBCh2Arm0Pl4","StripBCh2 Arm0Pl4",250,0.5,250.5));
  laStripBCh2Arm0Pl4->SetDirectory(0);

  laStripBCh2Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("StripBCh2Arm1Pl0","StripBCh2 Arm1Pl0",250,0.5,250.5));
  laStripBCh2Arm1Pl0->SetDirectory(0);
  laStripBCh2Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("StripBCh2Arm1Pl1","StripBCh2 Arm1Pl1",250,0.5,250.5));
  laStripBCh2Arm1Pl1->SetDirectory(0);
  laStripBCh2Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("StripBCh2Arm1Pl2","StripBCh2 Arm1Pl2",250,0.5,250.5));
  laStripBCh2Arm1Pl2->SetDirectory(0);
  laStripBCh2Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("StripBCh2Arm1Pl3","StripBCh2 Arm1Pl3",250,0.5,250.5));
  laStripBCh2Arm1Pl3->SetDirectory(0);
  laStripBCh2Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("StripBCh2Arm1Pl4","StripBCh2 Arm1Pl4",250,0.5,250.5));
  laStripBCh2Arm1Pl4->SetDirectory(0);


  ilWireCh3Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("WireCh3Arm0Pl0","WireCh3 Arm0Pl0",250,0.5,250.5));
  ilWireCh3Arm0Pl0->SetDirectory(0);
  ilWireCh3Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("WireCh3Arm0Pl1","WireCh3 Arm0Pl1",250,0.5,250.5));
  ilWireCh3Arm0Pl1->SetDirectory(0);
  ilWireCh3Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("WireCh3Arm0Pl2","WireCh3 Arm0Pl2",250,0.5,250.5));
  ilWireCh3Arm0Pl2->SetDirectory(0);
  ilWireCh3Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("WireCh3Arm0Pl3","WireCh3 Arm0Pl3",250,0.5,250.5));
  ilWireCh3Arm0Pl3->SetDirectory(0);
  ilWireCh3Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("WireCh3Arm0Pl4","WireCh3 Arm0Pl4",250,0.5,250.5));
  ilWireCh3Arm0Pl4->SetDirectory(0);

  ilWireCh3Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("WireCh3Arm1Pl0","WireCh3 Arm1Pl0",250,0.5,250.5));
  ilWireCh3Arm1Pl0->SetDirectory(0);
  ilWireCh3Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("WireCh3Arm1Pl1","WireCh3 Arm1Pl1",250,0.5,250.5));
  ilWireCh3Arm1Pl1->SetDirectory(0);
  ilWireCh3Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("WireCh3Arm1Pl2","WireCh3 Arm1Pl2",250,0.5,250.5));
  ilWireCh3Arm1Pl2->SetDirectory(0);
  ilWireCh3Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("WireCh3Arm1Pl3","WireCh3 Arm1Pl3",250,0.5,250.5));
  ilWireCh3Arm1Pl3->SetDirectory(0);
  ilWireCh3Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("WireCh3Arm1Pl4","WireCh3 Arm1Pl4",250,0.5,250.5));
  ilWireCh3Arm1Pl4->SetDirectory(0);






  ilClWidthCh3Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh3Arm0Pl0","ClWidthCh3 Arm0Pl0",250,0.5,250.5));
  ilClWidthCh3Arm0Pl0->SetDirectory(0);
  ilClWidthCh3Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh3Arm0Pl1","ClWidthCh3 Arm0Pl1",250,0.5,250.5));
  ilClWidthCh3Arm0Pl1->SetDirectory(0);
  ilClWidthCh3Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh3Arm0Pl2","ClWidthCh3 Arm0Pl2",250,0.5,250.5));
  ilClWidthCh3Arm0Pl2->SetDirectory(0);
  ilClWidthCh3Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh3Arm0Pl3","ClWidthCh3 Arm0Pl3",250,0.5,250.5));
  ilClWidthCh3Arm0Pl3->SetDirectory(0);
  ilClWidthCh3Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh3Arm0Pl4","ClWidthCh3 Arm0Pl4",250,0.5,250.5));
  ilClWidthCh3Arm0Pl4->SetDirectory(0);

  ilClWidthCh3Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh3Arm1Pl0","ClWidthCh3 Arm1Pl0",250,0.5,250.5));
  ilClWidthCh3Arm1Pl0->SetDirectory(0);
  ilClWidthCh3Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh3Arm1Pl1","ClWidthCh3 Arm1Pl1",250,0.5,250.5));
  ilClWidthCh3Arm1Pl1->SetDirectory(0);
  ilClWidthCh3Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh3Arm1Pl2","ClWidthCh3 Arm1Pl2",250,0.5,250.5));
  ilClWidthCh3Arm1Pl2->SetDirectory(0);
  ilClWidthCh3Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh3Arm1Pl3","ClWidthCh3 Arm1Pl3",250,0.5,250.5));
  ilClWidthCh3Arm1Pl3->SetDirectory(0);
  ilClWidthCh3Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh3Arm1Pl4","ClWidthCh3 Arm1Pl4",250,0.5,250.5));
  ilClWidthCh3Arm1Pl4->SetDirectory(0);






  laStripACh3Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("StripACh3Arm0Pl0","StripACh3 Arm0Pl0",250,0.5,250.5));
  laStripACh3Arm0Pl0->SetDirectory(0);
  laStripACh3Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("StripACh3Arm0Pl1","StripACh3 Arm0Pl1",250,0.5,250.5));
  laStripACh3Arm0Pl1->SetDirectory(0);
  laStripACh3Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("StripACh3Arm0Pl2","StripACh3 Arm0Pl2",250,0.5,250.5));
  laStripACh3Arm0Pl2->SetDirectory(0);
  laStripACh3Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("StripACh3Arm0Pl3","StripACh3 Arm0Pl3",250,0.5,250.5));
  laStripACh3Arm0Pl3->SetDirectory(0);
  laStripACh3Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("StripACh3Arm0Pl4","StripACh3 Arm0Pl4",250,0.5,250.5));
  laStripACh3Arm0Pl4->SetDirectory(0);

  laStripACh3Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("StripACh3Arm1Pl0","StripACh3 Arm1Pl0",250,0.5,250.5));
  laStripACh3Arm1Pl0->SetDirectory(0);
  laStripACh3Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("StripACh3Arm1Pl1","StripACh3 Arm1Pl1",250,0.5,250.5));
  laStripACh3Arm1Pl1->SetDirectory(0);
  laStripACh3Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("StripACh3Arm1Pl2","StripACh3 Arm1Pl2",250,0.5,250.5));
  laStripACh3Arm1Pl2->SetDirectory(0);
  laStripACh3Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("StripACh3Arm1Pl3","StripACh3 Arm1Pl3",250,0.5,250.5));
  laStripACh3Arm1Pl3->SetDirectory(0);
  laStripACh3Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("StripACh3Arm1Pl4","StripACh3 Arm1Pl4",250,0.5,250.5));
  laStripACh3Arm1Pl4->SetDirectory(0);


  laStripBCh3Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("StripBCh3Arm0Pl0","StripBCh3 Arm0Pl0",250,0.5,250.5));
  laStripBCh3Arm0Pl0->SetDirectory(0);
  laStripBCh3Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("StripBCh3Arm0Pl1","StripBCh3 Arm0Pl1",250,0.5,250.5));
  laStripBCh3Arm0Pl1->SetDirectory(0);
  laStripBCh3Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("StripBCh3Arm0Pl2","StripBCh3 Arm0Pl2",250,0.5,250.5));
  laStripBCh3Arm0Pl2->SetDirectory(0);
  laStripBCh3Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("StripBCh3Arm0Pl3","StripBCh3 Arm0Pl3",250,0.5,250.5));
  laStripBCh3Arm0Pl3->SetDirectory(0);
  laStripBCh3Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("StripBCh3Arm0Pl4","StripBCh3 Arm0Pl4",250,0.5,250.5));
  laStripBCh3Arm0Pl4->SetDirectory(0);

  laStripBCh3Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("StripBCh3Arm1Pl0","StripBCh3 Arm1Pl0",250,0.5,250.5));
  laStripBCh3Arm1Pl0->SetDirectory(0);
  laStripBCh3Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("StripBCh3Arm1Pl1","StripBCh3 Arm1Pl1",250,0.5,250.5));
  laStripBCh3Arm1Pl1->SetDirectory(0);
  laStripBCh3Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("StripBCh3Arm1Pl2","StripBCh3 Arm1Pl2",250,0.5,250.5));
  laStripBCh3Arm1Pl2->SetDirectory(0);
  laStripBCh3Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("StripBCh3Arm1Pl3","StripBCh3 Arm1Pl3",250,0.5,250.5));
  laStripBCh3Arm1Pl3->SetDirectory(0);
  laStripBCh3Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("StripBCh3Arm1Pl4","StripBCh3 Arm1Pl4",250,0.5,250.5));
  laStripBCh3Arm1Pl4->SetDirectory(0);


  ilWireCh4Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("WireCh4Arm0Pl0","WireCh4 Arm0Pl0",250,0.5,250.5));
  ilWireCh4Arm0Pl0->SetDirectory(0);
  ilWireCh4Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("WireCh4Arm0Pl1","WireCh4 Arm0Pl1",250,0.5,250.5));
  ilWireCh4Arm0Pl1->SetDirectory(0);
  ilWireCh4Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("WireCh4Arm0Pl2","WireCh4 Arm0Pl2",250,0.5,250.5));
  ilWireCh4Arm0Pl2->SetDirectory(0);
  ilWireCh4Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("WireCh4Arm0Pl3","WireCh4 Arm0Pl3",250,0.5,250.5));
  ilWireCh4Arm0Pl3->SetDirectory(0);
  ilWireCh4Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("WireCh4Arm0Pl4","WireCh4 Arm0Pl4",250,0.5,250.5));
  ilWireCh4Arm0Pl4->SetDirectory(0);

  ilWireCh4Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("WireCh4Arm1Pl0","WireCh4 Arm1Pl0",250,0.5,250.5));
  ilWireCh4Arm1Pl0->SetDirectory(0);
  ilWireCh4Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("WireCh4Arm1Pl1","WireCh4 Arm1Pl1",250,0.5,250.5));
  ilWireCh4Arm1Pl1->SetDirectory(0);
  ilWireCh4Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("WireCh4Arm1Pl2","WireCh4 Arm1Pl2",250,0.5,250.5));
  ilWireCh4Arm1Pl2->SetDirectory(0);
  ilWireCh4Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("WireCh4Arm1Pl3","WireCh4 Arm1Pl3",250,0.5,250.5));
  ilWireCh4Arm1Pl3->SetDirectory(0);
  ilWireCh4Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("WireCh4Arm1Pl4","WireCh4 Arm1Pl4",250,0.5,250.5));
  ilWireCh4Arm1Pl4->SetDirectory(0);




  ilClWidthCh4Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh4Arm0Pl0","ClWidthCh4 Arm0Pl0",250,0.5,250.5));
  ilClWidthCh4Arm0Pl0->SetDirectory(0);
  ilClWidthCh4Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh4Arm0Pl1","ClWidthCh4 Arm0Pl1",250,0.5,250.5));
  ilClWidthCh4Arm0Pl1->SetDirectory(0);
  ilClWidthCh4Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh4Arm0Pl2","ClWidthCh4 Arm0Pl2",250,0.5,250.5));
  ilClWidthCh4Arm0Pl2->SetDirectory(0);
  ilClWidthCh4Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh4Arm0Pl3","ClWidthCh4 Arm0Pl3",250,0.5,250.5));
  ilClWidthCh4Arm0Pl3->SetDirectory(0);
  ilClWidthCh4Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh4Arm0Pl4","ClWidthCh4 Arm0Pl4",250,0.5,250.5));
  ilClWidthCh4Arm0Pl4->SetDirectory(0);

  ilClWidthCh4Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh4Arm1Pl0","ClWidthCh4 Arm1Pl0",250,0.5,250.5));
  ilClWidthCh4Arm1Pl0->SetDirectory(0);
  ilClWidthCh4Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh4Arm1Pl1","ClWidthCh4 Arm1Pl1",250,0.5,250.5));
  ilClWidthCh4Arm1Pl1->SetDirectory(0);
  ilClWidthCh4Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh4Arm1Pl2","ClWidthCh4 Arm1Pl2",250,0.5,250.5));
  ilClWidthCh4Arm1Pl2->SetDirectory(0);
  ilClWidthCh4Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh4Arm1Pl3","ClWidthCh4 Arm1Pl3",250,0.5,250.5));
  ilClWidthCh4Arm1Pl3->SetDirectory(0);
  ilClWidthCh4Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh4Arm1Pl4","ClWidthCh4 Arm1Pl4",250,0.5,250.5));
  ilClWidthCh4Arm1Pl4->SetDirectory(0);



  laStripACh4Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("StripACh4Arm0Pl0","StripACh4 Arm0Pl0",250,0.5,250.5));
  laStripACh4Arm0Pl0->SetDirectory(0);
  laStripACh4Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("StripACh4Arm0Pl1","StripACh4 Arm0Pl1",250,0.5,250.5));
  laStripACh4Arm0Pl1->SetDirectory(0);
  laStripACh4Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("StripACh4Arm0Pl2","StripACh4 Arm0Pl2",250,0.5,250.5));
  laStripACh4Arm0Pl2->SetDirectory(0);
  laStripACh4Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("StripACh4Arm0Pl3","StripACh4 Arm0Pl3",250,0.5,250.5));
  laStripACh4Arm0Pl3->SetDirectory(0);
  laStripACh4Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("StripACh4Arm0Pl4","StripACh4 Arm0Pl4",250,0.5,250.5));
  laStripACh4Arm0Pl4->SetDirectory(0);

  laStripACh4Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("StripACh4Arm1Pl0","StripACh4 Arm1Pl0",250,0.5,250.5));
  laStripACh4Arm1Pl0->SetDirectory(0);
  laStripACh4Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("StripACh4Arm1Pl1","StripACh4 Arm1Pl1",250,0.5,250.5));
  laStripACh4Arm1Pl1->SetDirectory(0);
  laStripACh4Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("StripACh4Arm1Pl2","StripACh4 Arm1Pl2",250,0.5,250.5));
  laStripACh4Arm1Pl2->SetDirectory(0);
  laStripACh4Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("StripACh4Arm1Pl3","StripACh4 Arm1Pl3",250,0.5,250.5));
  laStripACh4Arm1Pl3->SetDirectory(0);
  laStripACh4Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("StripACh4Arm1Pl4","StripACh4 Arm1Pl4",250,0.5,250.5));
  laStripACh4Arm1Pl4->SetDirectory(0);


  laStripBCh4Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("StripBCh4Arm0Pl0","StripBCh4 Arm0Pl0",250,0.5,250.5));
  laStripBCh4Arm0Pl0->SetDirectory(0);
  laStripBCh4Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("StripBCh4Arm0Pl1","StripBCh4 Arm0Pl1",250,0.5,250.5));
  laStripBCh4Arm0Pl1->SetDirectory(0);
  laStripBCh4Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("StripBCh4Arm0Pl2","StripBCh4 Arm0Pl2",250,0.5,250.5));
  laStripBCh4Arm0Pl2->SetDirectory(0);
  laStripBCh4Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("StripBCh4Arm0Pl3","StripBCh4 Arm0Pl3",250,0.5,250.5));
  laStripBCh4Arm0Pl3->SetDirectory(0);
  laStripBCh4Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("StripBCh4Arm0Pl4","StripBCh4 Arm0Pl4",250,0.5,250.5));
  laStripBCh4Arm0Pl4->SetDirectory(0);

  laStripBCh4Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("StripBCh4Arm1Pl0","StripBCh4 Arm1Pl0",250,0.5,250.5));
  laStripBCh4Arm1Pl0->SetDirectory(0);
  laStripBCh4Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("StripBCh4Arm1Pl1","StripBCh4 Arm1Pl1",250,0.5,250.5));
  laStripBCh4Arm1Pl1->SetDirectory(0);
  laStripBCh4Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("StripBCh4Arm1Pl2","StripBCh4 Arm1Pl2",250,0.5,250.5));
  laStripBCh4Arm1Pl2->SetDirectory(0);
  laStripBCh4Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("StripBCh4Arm1Pl3","StripBCh4 Arm1Pl3",250,0.5,250.5));
  laStripBCh4Arm1Pl3->SetDirectory(0);
  laStripBCh4Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("StripBCh4Arm1Pl4","StripBCh4 Arm1Pl4",250,0.5,250.5));
  laStripBCh4Arm1Pl4->SetDirectory(0);


  ilWireCh5Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("WireCh5Arm0Pl0","WireCh5 Arm0Pl0",250,0.5,250.5));
  ilWireCh5Arm0Pl0->SetDirectory(0);
  ilWireCh5Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("WireCh5Arm0Pl1","WireCh5 Arm0Pl1",250,0.5,250.5));
  ilWireCh5Arm0Pl1->SetDirectory(0);
  ilWireCh5Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("WireCh5Arm0Pl2","WireCh5 Arm0Pl2",250,0.5,250.5));
  ilWireCh5Arm0Pl2->SetDirectory(0);
  ilWireCh5Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("WireCh5Arm0Pl3","WireCh5 Arm0Pl3",250,0.5,250.5));
  ilWireCh5Arm0Pl3->SetDirectory(0);
  ilWireCh5Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("WireCh5Arm0Pl4","WireCh5 Arm0Pl4",250,0.5,250.5));
  ilWireCh5Arm0Pl4->SetDirectory(0);

  ilWireCh5Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("WireCh5Arm1Pl0","WireCh5 Arm1Pl0",250,0.5,250.5));
  ilWireCh5Arm1Pl0->SetDirectory(0);
  ilWireCh5Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("WireCh5Arm1Pl1","WireCh5 Arm1Pl1",250,0.5,250.5));
  ilWireCh5Arm1Pl1->SetDirectory(0);
  ilWireCh5Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("WireCh5Arm1Pl2","WireCh5 Arm1Pl2",250,0.5,250.5));
  ilWireCh5Arm1Pl2->SetDirectory(0);
  ilWireCh5Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("WireCh5Arm1Pl3","WireCh5 Arm1Pl3",250,0.5,250.5));
  ilWireCh5Arm1Pl3->SetDirectory(0);
  ilWireCh5Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("WireCh5Arm1Pl4","WireCh5 Arm1Pl4",250,0.5,250.5));
  ilWireCh5Arm1Pl4->SetDirectory(0);




  ilClWidthCh5Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh5Arm0Pl0","ClWidthCh5 Arm0Pl0",250,0.5,250.5));
  ilClWidthCh5Arm0Pl0->SetDirectory(0);
  ilClWidthCh5Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh5Arm0Pl1","ClWidthCh5 Arm0Pl1",250,0.5,250.5));
  ilClWidthCh5Arm0Pl1->SetDirectory(0);
  ilClWidthCh5Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh5Arm0Pl2","ClWidthCh5 Arm0Pl2",250,0.5,250.5));
  ilClWidthCh5Arm0Pl2->SetDirectory(0);
  ilClWidthCh5Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh5Arm0Pl3","ClWidthCh5 Arm0Pl3",250,0.5,250.5));
  ilClWidthCh5Arm0Pl3->SetDirectory(0);
  ilClWidthCh5Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh5Arm0Pl4","ClWidthCh5 Arm0Pl4",250,0.5,250.5));
  ilClWidthCh5Arm0Pl4->SetDirectory(0);

  ilClWidthCh5Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh5Arm1Pl0","ClWidthCh5 Arm1Pl0",250,0.5,250.5));
  ilClWidthCh5Arm1Pl0->SetDirectory(0);
  ilClWidthCh5Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh5Arm1Pl1","ClWidthCh5 Arm1Pl1",250,0.5,250.5));
  ilClWidthCh5Arm1Pl1->SetDirectory(0);
  ilClWidthCh5Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh5Arm1Pl2","ClWidthCh5 Arm1Pl2",250,0.5,250.5));
  ilClWidthCh5Arm1Pl2->SetDirectory(0);
  ilClWidthCh5Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh5Arm1Pl3","ClWidthCh5 Arm1Pl3",250,0.5,250.5));
  ilClWidthCh5Arm1Pl3->SetDirectory(0);
  ilClWidthCh5Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("ClWidthCh5Arm1Pl4","ClWidthCh5 Arm1Pl4",250,0.5,250.5));
  ilClWidthCh5Arm1Pl4->SetDirectory(0);



  laStripACh5Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("StripACh5Arm0Pl0","StripACh5 Arm0Pl0",250,0.5,250.5));
  laStripACh5Arm0Pl0->SetDirectory(0);
  laStripACh5Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("StripACh5Arm0Pl1","StripACh5 Arm0Pl1",250,0.5,250.5));
  laStripACh5Arm0Pl1->SetDirectory(0);
  laStripACh5Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("StripACh5Arm0Pl2","StripACh5 Arm0Pl2",250,0.5,250.5));
  laStripACh5Arm0Pl2->SetDirectory(0);
  laStripACh5Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("StripACh5Arm0Pl3","StripACh5 Arm0Pl3",250,0.5,250.5));
  laStripACh5Arm0Pl3->SetDirectory(0);
  laStripACh5Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("StripACh5Arm0Pl4","StripACh5 Arm0Pl4",250,0.5,250.5));
  laStripACh5Arm0Pl4->SetDirectory(0);

  laStripACh5Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("StripACh5Arm1Pl0","StripACh5 Arm1Pl0",250,0.5,250.5));
  laStripACh5Arm1Pl0->SetDirectory(0);
  laStripACh5Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("StripACh5Arm1Pl1","StripACh5 Arm1Pl1",250,0.5,250.5));
  laStripACh5Arm1Pl1->SetDirectory(0);
  laStripACh5Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("StripACh5Arm1Pl2","StripACh5 Arm1Pl2",250,0.5,250.5));
  laStripACh5Arm1Pl2->SetDirectory(0);
  laStripACh5Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("StripACh5Arm1Pl3","StripACh5 Arm1Pl3",250,0.5,250.5));
  laStripACh5Arm1Pl3->SetDirectory(0);
  laStripACh5Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("StripACh5Arm1Pl4","StripACh5 Arm1Pl4",250,0.5,250.5));
  laStripACh5Arm1Pl4->SetDirectory(0);


  laStripBCh5Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("StripBCh5Arm0Pl0","StripBCh5 Arm0Pl0",250,0.5,250.5));
  laStripBCh5Arm0Pl0->SetDirectory(0);
  laStripBCh5Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("StripBCh5Arm0Pl1","StripBCh5 Arm0Pl1",250,0.5,250.5));
  laStripBCh5Arm0Pl1->SetDirectory(0);
  laStripBCh5Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("StripBCh5Arm0Pl2","StripBCh5 Arm0Pl2",250,0.5,250.5));
  laStripBCh5Arm0Pl2->SetDirectory(0);
  laStripBCh5Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("StripBCh5Arm0Pl3","StripBCh5 Arm0Pl3",250,0.5,250.5));
  laStripBCh5Arm0Pl3->SetDirectory(0);
  laStripBCh5Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("StripBCh5Arm0Pl4","StripBCh5 Arm0Pl4",250,0.5,250.5));
  laStripBCh5Arm0Pl4->SetDirectory(0);

  laStripBCh5Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("StripBCh5Arm1Pl0","StripBCh5 Arm1Pl0",250,0.5,250.5));
  laStripBCh5Arm1Pl0->SetDirectory(0);
  laStripBCh5Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("StripBCh5Arm1Pl1","StripBCh5 Arm1Pl1",250,0.5,250.5));
  laStripBCh5Arm1Pl1->SetDirectory(0);
  laStripBCh5Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("StripBCh5Arm1Pl2","StripBCh5 Arm1Pl2",250,0.5,250.5));
  laStripBCh5Arm1Pl2->SetDirectory(0);
  laStripBCh5Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("StripBCh5Arm1Pl3","StripBCh5 Arm1Pl3",250,0.5,250.5));
  laStripBCh5Arm1Pl3->SetDirectory(0);
  laStripBCh5Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("StripBCh5Arm1Pl4","StripBCh5 Arm1Pl4",250,0.5,250.5));
  laStripBCh5Arm1Pl4->SetDirectory(0);
//////////////////////////////////////////////////

  NWireCh0Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NWireCh0Arm0Pl0","NWireCh0 Arm0Pl0",250,-0.5,249.5));
  NWireCh0Arm0Pl0->SetDirectory(0);
  NWireCh0Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NWireCh0Arm0Pl1","NWireCh0 Arm0Pl1",250,-0.5,249.5));
  NWireCh0Arm0Pl1->SetDirectory(0);
  NWireCh0Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NWireCh0Arm0Pl2","NWireCh0 Arm0Pl2",250,-0.5,249.5));
  NWireCh0Arm0Pl2->SetDirectory(0);
  NWireCh0Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NWireCh0Arm0Pl3","NWireCh0 Arm0Pl3",250,-0.5,249.5));
  NWireCh0Arm0Pl3->SetDirectory(0);
  NWireCh0Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NWireCh0Arm0Pl4","NWireCh0 Arm0Pl4",250,-0.5,249.5));
  NWireCh0Arm0Pl4->SetDirectory(0);

  NWireCh0Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NWireCh0Arm1Pl0","NWireCh0 Arm1Pl0",250,-0.5,249.5));
  NWireCh0Arm1Pl0->SetDirectory(0);
  NWireCh0Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NWireCh0Arm1Pl1","NWireCh0 Arm1Pl1",250,-0.5,249.5));
  NWireCh0Arm1Pl1->SetDirectory(0);
  NWireCh0Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NWireCh0Arm1Pl2","NWireCh0 Arm1Pl2",250,-0.5,249.5));
  NWireCh0Arm1Pl2->SetDirectory(0);
  NWireCh0Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NWireCh0Arm1Pl3","NWireCh0 Arm1Pl3",250,-0.5,249.5));
  NWireCh0Arm1Pl3->SetDirectory(0);
  NWireCh0Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NWireCh0Arm1Pl4","NWireCh0 Arm1Pl4",250,-0.5,249.5));
  NWireCh0Arm1Pl4->SetDirectory(0);

  NStripACh0Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripACh0Arm0Pl0","NStripACh0 Arm0Pl0",250,-0.5,249.5));
  NStripACh0Arm0Pl0->SetDirectory(0);
  NStripACh0Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripACh0Arm0Pl1","NStripACh0 Arm0Pl1",250,-0.5,249.5));
  NStripACh0Arm0Pl1->SetDirectory(0);
  NStripACh0Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripACh0Arm0Pl2","NStripACh0 Arm0Pl2",250,-0.5,249.5));
  NStripACh0Arm0Pl2->SetDirectory(0);
  NStripACh0Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripACh0Arm0Pl3","NStripACh0 Arm0Pl3",250,-0.5,249.5));
  NStripACh0Arm0Pl3->SetDirectory(0);
  NStripACh0Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripACh0Arm0Pl4","NStripACh0 Arm0Pl4",250,-0.5,249.5));
  NStripACh0Arm0Pl4->SetDirectory(0);

  NStripACh0Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripACh0Arm1Pl0","NStripACh0 Arm1Pl0",250,-0.5,249.5));
  NStripACh0Arm1Pl0->SetDirectory(0);
  NStripACh0Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripACh0Arm1Pl1","NStripACh0 Arm1Pl1",250,-0.5,249.5));
  NStripACh0Arm1Pl1->SetDirectory(0);
  NStripACh0Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripACh0Arm1Pl2","NStripACh0 Arm1Pl2",250,-0.5,249.5));
  NStripACh0Arm1Pl2->SetDirectory(0);
  NStripACh0Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripACh0Arm1Pl3","NStripACh0 Arm1Pl3",250,-0.5,249.5));
  NStripACh0Arm1Pl3->SetDirectory(0);
  NStripACh0Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripACh0Arm1Pl4","NStripACh0 Arm1Pl4",250,-0.5,249.5));
  NStripACh0Arm1Pl4->SetDirectory(0);


  NStripBCh0Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripBCh0Arm0Pl0","NStripBCh0 Arm0Pl0",250,-0.5,249.5));
  NStripBCh0Arm0Pl0->SetDirectory(0);
  NStripBCh0Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripBCh0Arm0Pl1","NStripBCh0 Arm0Pl1",250,-0.5,249.5));
  NStripBCh0Arm0Pl1->SetDirectory(0);
  NStripBCh0Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripBCh0Arm0Pl2","NStripBCh0 Arm0Pl2",250,-0.5,249.5));
  NStripBCh0Arm0Pl2->SetDirectory(0);
  NStripBCh0Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripBCh0Arm0Pl3","NStripBCh0 Arm0Pl3",250,-0.5,249.5));
  NStripBCh0Arm0Pl3->SetDirectory(0);
  NStripBCh0Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripBCh0Arm0Pl4","NStripBCh0 Arm0Pl4",250,-0.5,249.5));
  NStripBCh0Arm0Pl4->SetDirectory(0);

  NStripBCh0Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripBCh0Arm1Pl0","NStripBCh0 Arm1Pl0",250,-0.5,249.5));
  NStripBCh0Arm1Pl0->SetDirectory(0);
  NStripBCh0Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripBCh0Arm1Pl1","NStripBCh0 Arm1Pl1",250,-0.5,249.5));
  NStripBCh0Arm1Pl1->SetDirectory(0);
  NStripBCh0Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripBCh0Arm1Pl2","NStripBCh0 Arm1Pl2",250,-0.5,249.5));
  NStripBCh0Arm1Pl2->SetDirectory(0);
  NStripBCh0Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripBCh0Arm1Pl3","NStripBCh0 Arm1Pl3",250,-0.5,249.5));
  NStripBCh0Arm1Pl3->SetDirectory(0);
  NStripBCh0Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripBCh0Arm1Pl4","NStripBCh0 Arm1Pl4",250,-0.5,249.5));
  NStripBCh0Arm1Pl4->SetDirectory(0);


  NWireCh1Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NWireCh1Arm0Pl0","NWireCh1 Arm0Pl0",250,-0.5,249.5));
  NWireCh1Arm0Pl0->SetDirectory(0);
  NWireCh1Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NWireCh1Arm0Pl1","NWireCh1 Arm0Pl1",250,-0.5,249.5));
  NWireCh1Arm0Pl1->SetDirectory(0);
  NWireCh1Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NWireCh1Arm0Pl2","NWireCh1 Arm0Pl2",250,-0.5,249.5));
  NWireCh1Arm0Pl2->SetDirectory(0);
  NWireCh1Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NWireCh1Arm0Pl3","NWireCh1 Arm0Pl3",250,-0.5,249.5));
  NWireCh1Arm0Pl3->SetDirectory(0);
  NWireCh1Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NWireCh1Arm0Pl4","NWireCh1 Arm0Pl4",250,-0.5,249.5));
  NWireCh1Arm0Pl4->SetDirectory(0);

  NWireCh1Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NWireCh1Arm1Pl0","NWireCh1 Arm1Pl0",250,-0.5,249.5));
  NWireCh1Arm1Pl0->SetDirectory(0);
  NWireCh1Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NWireCh1Arm1Pl1","NWireCh1 Arm1Pl1",250,-0.5,249.5));
  NWireCh1Arm1Pl1->SetDirectory(0);
  NWireCh1Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NWireCh1Arm1Pl2","NWireCh1 Arm1Pl2",250,-0.5,249.5));
  NWireCh1Arm1Pl2->SetDirectory(0);
  NWireCh1Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NWireCh1Arm1Pl3","NWireCh1 Arm1Pl3",250,-0.5,249.5));
  NWireCh1Arm1Pl3->SetDirectory(0);
  NWireCh1Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NWireCh1Arm1Pl4","NWireCh1 Arm1Pl4",250,-0.5,249.5));
  NWireCh1Arm1Pl4->SetDirectory(0);

  NStripACh1Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripACh1Arm0Pl0","NStripACh1 Arm0Pl0",250,-0.5,249.5));
  NStripACh1Arm0Pl0->SetDirectory(0);
  NStripACh1Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripACh1Arm0Pl1","NStripACh1 Arm0Pl1",250,-0.5,249.5));
  NStripACh1Arm0Pl1->SetDirectory(0);
  NStripACh1Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripACh1Arm0Pl2","NStripACh1 Arm0Pl2",250,-0.5,249.5));
  NStripACh1Arm0Pl2->SetDirectory(0);
  NStripACh1Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripACh1Arm0Pl3","NStripACh1 Arm0Pl3",250,-0.5,249.5));
  NStripACh1Arm0Pl3->SetDirectory(0);
  NStripACh1Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripACh1Arm0Pl4","NStripACh1 Arm0Pl4",250,-0.5,249.5));
  NStripACh1Arm0Pl4->SetDirectory(0);

  NStripACh1Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripACh1Arm1Pl0","NStripACh1 Arm1Pl0",250,-0.5,249.5));
  NStripACh1Arm1Pl0->SetDirectory(0);
  NStripACh1Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripACh1Arm1Pl1","NStripACh1 Arm1Pl1",250,-0.5,249.5));
  NStripACh1Arm1Pl1->SetDirectory(0);
  NStripACh1Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripACh1Arm1Pl2","NStripACh1 Arm1Pl2",250,-0.5,249.5));
  NStripACh1Arm1Pl2->SetDirectory(0);
  NStripACh1Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripACh1Arm1Pl3","NStripACh1 Arm1Pl3",250,-0.5,249.5));
  NStripACh1Arm1Pl3->SetDirectory(0);
  NStripACh1Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripACh1Arm1Pl4","NStripACh1 Arm1Pl4",250,-0.5,249.5));
  NStripACh1Arm1Pl4->SetDirectory(0);


  NStripBCh1Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripBCh1Arm0Pl0","NStripBCh1 Arm0Pl0",250,-0.5,249.5));
  NStripBCh1Arm0Pl0->SetDirectory(0);
  NStripBCh1Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripBCh1Arm0Pl1","NStripBCh1 Arm0Pl1",250,-0.5,249.5));
  NStripBCh1Arm0Pl1->SetDirectory(0);
  NStripBCh1Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripBCh1Arm0Pl2","NStripBCh1 Arm0Pl2",250,-0.5,249.5));
  NStripBCh1Arm0Pl2->SetDirectory(0);
  NStripBCh1Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripBCh1Arm0Pl3","NStripBCh1 Arm0Pl3",250,-0.5,249.5));
  NStripBCh1Arm0Pl3->SetDirectory(0);
  NStripBCh1Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripBCh1Arm0Pl4","NStripBCh1 Arm0Pl4",250,-0.5,249.5));
  NStripBCh1Arm0Pl4->SetDirectory(0);

  NStripBCh1Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripBCh1Arm1Pl0","NStripBCh1 Arm1Pl0",250,-0.5,249.5));
  NStripBCh1Arm1Pl0->SetDirectory(0);
  NStripBCh1Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripBCh1Arm1Pl1","NStripBCh1 Arm1Pl1",250,-0.5,249.5));
  NStripBCh1Arm1Pl1->SetDirectory(0);
  NStripBCh1Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripBCh1Arm1Pl2","NStripBCh1 Arm1Pl2",250,-0.5,249.5));
  NStripBCh1Arm1Pl2->SetDirectory(0);
  NStripBCh1Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripBCh1Arm1Pl3","NStripBCh1 Arm1Pl3",250,-0.5,249.5));
  NStripBCh1Arm1Pl3->SetDirectory(0);
  NStripBCh1Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripBCh1Arm1Pl4","NStripBCh1 Arm1Pl4",250,-0.5,249.5));
  NStripBCh1Arm1Pl4->SetDirectory(0);


  NWireCh2Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NWireCh2Arm0Pl0","NWireCh2 Arm0Pl0",250,-0.5,249.5));
  NWireCh2Arm0Pl0->SetDirectory(0);
  NWireCh2Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NWireCh2Arm0Pl1","NWireCh2 Arm0Pl1",250,-0.5,249.5));
  NWireCh2Arm0Pl1->SetDirectory(0);
  NWireCh2Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NWireCh2Arm0Pl2","NWireCh2 Arm0Pl2",250,-0.5,249.5));
  NWireCh2Arm0Pl2->SetDirectory(0);
  NWireCh2Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NWireCh2Arm0Pl3","NWireCh2 Arm0Pl3",250,-0.5,249.5));
  NWireCh2Arm0Pl3->SetDirectory(0);
  NWireCh2Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NWireCh2Arm0Pl4","NWireCh2 Arm0Pl4",250,-0.5,249.5));
  NWireCh2Arm0Pl4->SetDirectory(0);

  NWireCh2Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NWireCh2Arm1Pl0","NWireCh2 Arm1Pl0",250,-0.5,249.5));
  NWireCh2Arm1Pl0->SetDirectory(0);
  NWireCh2Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NWireCh2Arm1Pl1","NWireCh2 Arm1Pl1",250,-0.5,249.5));
  NWireCh2Arm1Pl1->SetDirectory(0);
  NWireCh2Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NWireCh2Arm1Pl2","NWireCh2Arm1Pl2",250,-0.5,249.5));
  NWireCh2Arm1Pl2->SetDirectory(0);
  NWireCh2Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NWireCh2Arm1Pl3","NWireCh2Arm1Pl3",250,-0.5,249.5));
  NWireCh2Arm1Pl3->SetDirectory(0);
  NWireCh2Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NWireCh2Arm1Pl4","NWireCh2Arm1Pl4",250,-0.5,249.5));
  NWireCh2Arm1Pl4->SetDirectory(0);

  NStripACh2Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripACh2Arm0Pl0","NStripACh2Arm0Pl0",250,-0.5,249.5));
  NStripACh2Arm0Pl0->SetDirectory(0);
  NStripACh2Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripACh2Arm0Pl1","NStripACh2Arm0Pl1",250,-0.5,249.5));
  NStripACh2Arm0Pl1->SetDirectory(0);
  NStripACh2Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripACh2Arm0Pl2","NStripACh2Arm0Pl2",250,-0.5,249.5));
  NStripACh2Arm0Pl2->SetDirectory(0);
  NStripACh2Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripACh2Arm0Pl3","NStripACh2Arm0Pl3",250,-0.5,249.5));
  NStripACh2Arm0Pl3->SetDirectory(0);
  NStripACh2Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripACh2Arm0Pl4","NStripACh2Arm0Pl4",250,-0.5,249.5));
  NStripACh2Arm0Pl4->SetDirectory(0);

  NStripACh2Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripACh2Arm1Pl0","NStripACh2Arm1Pl0",250,-0.5,249.5));
  NStripACh2Arm1Pl0->SetDirectory(0);
  NStripACh2Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripACh2Arm1Pl1","NStripACh2Arm1Pl1",250,-0.5,249.5));
  NStripACh2Arm1Pl1->SetDirectory(0);
  NStripACh2Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripACh2Arm1Pl2","NStripACh2Arm1Pl2",250,-0.5,249.5));
  NStripACh2Arm1Pl2->SetDirectory(0);
  NStripACh2Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripACh2Arm1Pl3","NStripACh2Arm1Pl3",250,-0.5,249.5));
  NStripACh2Arm1Pl3->SetDirectory(0);
  NStripACh2Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripACh2Arm1Pl4","NStripACh2Arm1Pl4",250,-0.5,249.5));
  NStripACh2Arm1Pl4->SetDirectory(0);


  NStripBCh2Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripBCh2Arm0Pl0","NStripBCh2Arm0Pl0",250,-0.5,249.5));
  NStripBCh2Arm0Pl0->SetDirectory(0);
  NStripBCh2Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripBCh2Arm0Pl1","NStripBCh2Arm0Pl1",250,-0.5,249.5));
  NStripBCh2Arm0Pl1->SetDirectory(0);
  NStripBCh2Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripBCh2Arm0Pl2","NStripBCh2Arm0Pl2",250,-0.5,249.5));
  NStripBCh2Arm0Pl2->SetDirectory(0);
  NStripBCh2Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripBCh2Arm0Pl3","NStripBCh2Arm0Pl3",250,-0.5,249.5));
  NStripBCh2Arm0Pl3->SetDirectory(0);
  NStripBCh2Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripBCh2Arm0Pl4","NStripBCh2Arm0Pl4",250,-0.5,249.5));
  NStripBCh2Arm0Pl4->SetDirectory(0);

  NStripBCh2Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripBCh2Arm1Pl0","NStripBCh2Arm1Pl0",250,-0.5,249.5));
  NStripBCh2Arm1Pl0->SetDirectory(0);
  NStripBCh2Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripBCh2Arm1Pl1","NStripBCh2Arm1Pl1",250,-0.5,249.5));
  NStripBCh2Arm1Pl1->SetDirectory(0);
  NStripBCh2Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripBCh2Arm1Pl2","NStripBCh2Arm1Pl2",250,-0.5,249.5));
  NStripBCh2Arm1Pl2->SetDirectory(0);
  NStripBCh2Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripBCh2Arm1Pl3","NStripBCh2Arm1Pl3",250,-0.5,249.5));
  NStripBCh2Arm1Pl3->SetDirectory(0);
  NStripBCh2Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripBCh2Arm1Pl4","NStripBCh2Arm1Pl4",250,-0.5,249.5));
  NStripBCh2Arm1Pl4->SetDirectory(0);


  NWireCh3Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NWireCh3Arm0Pl0","NWireCh3Arm0Pl0",250,-0.5,249.5));
  NWireCh3Arm0Pl0->SetDirectory(0);
  NWireCh3Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NWireCh3Arm0Pl1","NWireCh3Arm0Pl1",250,-0.5,249.5));
  NWireCh3Arm0Pl1->SetDirectory(0);
  NWireCh3Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NWireCh3Arm0Pl2","NWireCh3Arm0Pl2",250,-0.5,249.5));
  NWireCh3Arm0Pl2->SetDirectory(0);
  NWireCh3Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NWireCh3Arm0Pl3","NWireCh3Arm0Pl3",250,-0.5,249.5));
  NWireCh3Arm0Pl3->SetDirectory(0);
  NWireCh3Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NWireCh3Arm0Pl4","NWireCh3Arm0Pl4",250,-0.5,249.5));
  NWireCh3Arm0Pl4->SetDirectory(0);

  NWireCh3Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NWireCh3Arm1Pl0","NWireCh3Arm1Pl0",250,-0.5,249.5));
  NWireCh3Arm1Pl0->SetDirectory(0);
  NWireCh3Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NWireCh3Arm1Pl1","NWireCh3Arm1Pl1",250,-0.5,249.5));
  NWireCh3Arm1Pl1->SetDirectory(0);
  NWireCh3Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NWireCh3Arm1Pl2","NWireCh3Arm1Pl2",250,-0.5,249.5));
  NWireCh3Arm1Pl2->SetDirectory(0);
  NWireCh3Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NWireCh3Arm1Pl3","NWireCh3Arm1Pl3",250,-0.5,249.5));
  NWireCh3Arm1Pl3->SetDirectory(0);
  NWireCh3Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NWireCh3Arm1Pl4","NWireCh3Arm1Pl4",250,-0.5,249.5));
  NWireCh3Arm1Pl4->SetDirectory(0);

  NStripACh3Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripACh3Arm0Pl0","NStripACh3Arm0Pl0",250,-0.5,249.5));
  NStripACh3Arm0Pl0->SetDirectory(0);
  NStripACh3Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripACh3Arm0Pl1","NStripACh3Arm0Pl1",250,-0.5,249.5));
  NStripACh3Arm0Pl1->SetDirectory(0);
  NStripACh3Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripACh3Arm0Pl2","NStripACh3Arm0Pl2",250,-0.5,249.5));
  NStripACh3Arm0Pl2->SetDirectory(0);
  NStripACh3Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripACh3Arm0Pl3","NStripACh3Arm0Pl3",250,-0.5,249.5));
  NStripACh3Arm0Pl3->SetDirectory(0);
  NStripACh3Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripACh3Arm0Pl4","NStripACh3Arm0Pl4",250,-0.5,249.5));
  NStripACh3Arm0Pl4->SetDirectory(0);

  NStripACh3Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripACh3Arm1Pl0","NStripACh3Arm1Pl0",250,-0.5,249.5));
  NStripACh3Arm1Pl0->SetDirectory(0);
  NStripACh3Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripACh3Arm1Pl1","NStripACh3Arm1Pl1",250,-0.5,249.5));
  NStripACh3Arm1Pl1->SetDirectory(0);
  NStripACh3Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripACh3Arm1Pl2","NStripACh3Arm1Pl2",250,-0.5,249.5));
  NStripACh3Arm1Pl2->SetDirectory(0);
  NStripACh3Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripACh3Arm1Pl3","NStripACh3Arm1Pl3",250,-0.5,249.5));
  NStripACh3Arm1Pl3->SetDirectory(0);
  NStripACh3Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripACh3Arm1Pl4","NStripACh3Arm1Pl4",250,-0.5,249.5));
  NStripACh3Arm1Pl4->SetDirectory(0);


  NStripBCh3Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripBCh3Arm0Pl0","NStripBCh3Arm0Pl0",250,-0.5,249.5));
  NStripBCh3Arm0Pl0->SetDirectory(0);
  NStripBCh3Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripBCh3Arm0Pl1","NStripBCh3Arm0Pl1",250,-0.5,249.5));
  NStripBCh3Arm0Pl1->SetDirectory(0);
  NStripBCh3Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripBCh3Arm0Pl2","NStripBCh3Arm0Pl2",250,-0.5,249.5));
  NStripBCh3Arm0Pl2->SetDirectory(0);
  NStripBCh3Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripBCh3Arm0Pl3","NStripBCh3Arm0Pl3",250,-0.5,249.5));
  NStripBCh3Arm0Pl3->SetDirectory(0);
  NStripBCh3Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripBCh3Arm0Pl4","NStripBCh3Arm0Pl4",250,-0.5,249.5));
  NStripBCh3Arm0Pl4->SetDirectory(0);

  NStripBCh3Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripBCh3Arm1Pl0","NStripBCh3Arm1Pl0",250,-0.5,249.5));
  NStripBCh3Arm1Pl0->SetDirectory(0);
  NStripBCh3Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripBCh3Arm1Pl1","NStripBCh3Arm1Pl1",250,-0.5,249.5));
  NStripBCh3Arm1Pl1->SetDirectory(0);
  NStripBCh3Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripBCh3Arm1Pl2","NStripBCh3Arm1Pl2",250,-0.5,249.5));
  NStripBCh3Arm1Pl2->SetDirectory(0);
  NStripBCh3Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripBCh3Arm1Pl3","NStripBCh3Arm1Pl3",250,-0.5,249.5));
  NStripBCh3Arm1Pl3->SetDirectory(0);
  NStripBCh3Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripBCh3Arm1Pl4","NStripBCh3Arm1Pl4",250,-0.5,249.5));
  NStripBCh3Arm1Pl4->SetDirectory(0);


  NWireCh4Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NWireCh4Arm0Pl0","NWireCh4Arm0Pl0",250,-0.5,249.5));
  NWireCh4Arm0Pl0->SetDirectory(0);
  NWireCh4Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NWireCh4Arm0Pl1","NWireCh4Arm0Pl1",250,-0.5,249.5));
  NWireCh4Arm0Pl1->SetDirectory(0);
  NWireCh4Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NWireCh4Arm0Pl2","NWireCh4Arm0Pl2",250,-0.5,249.5));
  NWireCh4Arm0Pl2->SetDirectory(0);
  NWireCh4Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NWireCh4Arm0Pl3","NWireCh4Arm0Pl3",250,-0.5,249.5));
  NWireCh4Arm0Pl3->SetDirectory(0);
  NWireCh4Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NWireCh4Arm0Pl4","NWireCh4Arm0Pl4",250,-0.5,249.5));
  NWireCh4Arm0Pl4->SetDirectory(0);

  NWireCh4Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NWireCh4Arm1Pl0","NWireCh4Arm1Pl0",250,-0.5,249.5));
  NWireCh4Arm1Pl0->SetDirectory(0);
  NWireCh4Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NWireCh4Arm1Pl1","NWireCh4Arm1Pl1",250,-0.5,249.5));
  NWireCh4Arm1Pl1->SetDirectory(0);
  NWireCh4Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NWireCh4Arm1Pl2","NWireCh4Arm1Pl2",250,-0.5,249.5));
  NWireCh4Arm1Pl2->SetDirectory(0);
  NWireCh4Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NWireCh4Arm1Pl3","NWireCh4Arm1Pl3",250,-0.5,249.5));
  NWireCh4Arm1Pl3->SetDirectory(0);
  NWireCh4Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NWireCh4Arm1Pl4","NWireCh4Arm1Pl4",250,-0.5,249.5));
  NWireCh4Arm1Pl4->SetDirectory(0);

  NStripACh4Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripACh4Arm0Pl0","NStripACh4Arm0Pl0",250,-0.5,249.5));
  NStripACh4Arm0Pl0->SetDirectory(0);
  NStripACh4Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripACh4Arm0Pl1","NStripACh4Arm0Pl1",250,-0.5,249.5));
  NStripACh4Arm0Pl1->SetDirectory(0);
  NStripACh4Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripACh4Arm0Pl2","NStripACh4Arm0Pl2",250,-0.5,249.5));
  NStripACh4Arm0Pl2->SetDirectory(0);
  NStripACh4Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripACh4Arm0Pl3","NStripACh4Arm0Pl3",250,-0.5,249.5));
  NStripACh4Arm0Pl3->SetDirectory(0);
  NStripACh4Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripACh4Arm0Pl4","NStripACh4Arm0Pl4",250,-0.5,249.5));
  NStripACh4Arm0Pl4->SetDirectory(0);

  NStripACh4Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripACh4Arm1Pl0","NStripACh4Arm1Pl0",250,-0.5,249.5));
  NStripACh4Arm1Pl0->SetDirectory(0);
  NStripACh4Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripACh4Arm1Pl1","NStripACh4Arm1Pl1",250,-0.5,249.5));
  NStripACh4Arm1Pl1->SetDirectory(0);
  NStripACh4Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripACh4Arm1Pl2","NStripACh4Arm1Pl2",250,-0.5,249.5));
  NStripACh4Arm1Pl2->SetDirectory(0);
  NStripACh4Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripACh4Arm1Pl3","NStripACh4Arm1Pl3",250,-0.5,249.5));
  NStripACh4Arm1Pl3->SetDirectory(0);
  NStripACh4Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripACh4Arm1Pl4","NStripACh4Arm1Pl4",250,-0.5,249.5));
  NStripACh4Arm1Pl4->SetDirectory(0);


  NStripBCh4Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripBCh4Arm0Pl0","NStripBCh4Arm0Pl0",250,-0.5,249.5));
  NStripBCh4Arm0Pl0->SetDirectory(0);
  NStripBCh4Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripBCh4Arm0Pl1","NStripBCh4Arm0Pl1",250,-0.5,249.5));
  NStripBCh4Arm0Pl1->SetDirectory(0);
  NStripBCh4Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripBCh4Arm0Pl2","NStripBCh4Arm0Pl2",250,-0.5,249.5));
  NStripBCh4Arm0Pl2->SetDirectory(0);
  NStripBCh4Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripBCh4Arm0Pl3","NStripBCh4Arm0Pl3",250,-0.5,249.5));
  NStripBCh4Arm0Pl3->SetDirectory(0);
  NStripBCh4Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripBCh4Arm0Pl4","NStripBCh4Arm0Pl4",250,-0.5,249.5));
  NStripBCh4Arm0Pl4->SetDirectory(0);

  NStripBCh4Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripBCh4Arm1Pl0","NStripBCh4Arm1Pl0",250,-0.5,249.5));
  NStripBCh4Arm1Pl0->SetDirectory(0);
  NStripBCh4Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripBCh4Arm1Pl1","NStripBCh4Arm1Pl1",250,-0.5,249.5));
  NStripBCh4Arm1Pl1->SetDirectory(0);
  NStripBCh4Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripBCh4Arm1Pl2","NStripBCh4Arm1Pl2",250,-0.5,249.5));
  NStripBCh4Arm1Pl2->SetDirectory(0);
  NStripBCh4Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripBCh4Arm1Pl3","NStripBCh4Arm1Pl3",250,-0.5,249.5));
  NStripBCh4Arm1Pl3->SetDirectory(0);
  NStripBCh4Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripBCh4Arm1Pl4","NStripBCh4Arm1Pl4",250,-0.5,249.5));
  NStripBCh4Arm1Pl4->SetDirectory(0);


  NWireCh5Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NWireCh5Arm0Pl0","NWireCh5Arm0Pl0",250,-0.5,249.5));
  NWireCh5Arm0Pl0->SetDirectory(0);
  NWireCh5Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NWireCh5Arm0Pl1","NWireCh5Arm0Pl1",250,-0.5,249.5));
  NWireCh5Arm0Pl1->SetDirectory(0);
  NWireCh5Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NWireCh5Arm0Pl2","NWireCh5Arm0Pl2",250,-0.5,249.5));
  NWireCh5Arm0Pl2->SetDirectory(0);
  NWireCh5Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NWireCh5Arm0Pl3","NWireCh5Arm0Pl3",250,-0.5,249.5));
  NWireCh5Arm0Pl3->SetDirectory(0);
  NWireCh5Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NWireCh5Arm0Pl4","NWireCh5Arm0Pl4",250,-0.5,249.5));
  NWireCh5Arm0Pl4->SetDirectory(0);

  NWireCh5Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NWireCh5Arm1Pl0","NWireCh5Arm1Pl0",250,-0.5,249.5));
  NWireCh5Arm1Pl0->SetDirectory(0);
  NWireCh5Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NWireCh5Arm1Pl1","NWireCh5Arm1Pl1",250,-0.5,249.5));
  NWireCh5Arm1Pl1->SetDirectory(0);
  NWireCh5Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NWireCh5Arm1Pl2","NWireCh5Arm1Pl2",250,-0.5,249.5));
  NWireCh5Arm1Pl2->SetDirectory(0);
  NWireCh5Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NWireCh5Arm1Pl3","NWireCh5Arm1Pl3",250,-0.5,249.5));
  NWireCh5Arm1Pl3->SetDirectory(0);
  NWireCh5Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NWireCh5Arm1Pl4","NWireCh5Arm1Pl4",250,-0.5,249.5));
  NWireCh5Arm1Pl4->SetDirectory(0);

  NStripACh5Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripACh5Arm0Pl0","NStripACh5Arm0Pl0",250,-0.5,249.5));
  NStripACh5Arm0Pl0->SetDirectory(0);
  NStripACh5Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripACh5Arm0Pl1","NStripACh5Arm0Pl1",250,-0.5,249.5));
  NStripACh5Arm0Pl1->SetDirectory(0);
  NStripACh5Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripACh5Arm0Pl2","NStripACh5Arm0Pl2",250,-0.5,249.5));
  NStripACh5Arm0Pl2->SetDirectory(0);
  NStripACh5Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripACh5Arm0Pl3","NStripACh5Arm0Pl3",250,-0.5,249.5));
  NStripACh5Arm0Pl3->SetDirectory(0);
  NStripACh5Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripACh5Arm0Pl4","NStripACh5Arm0Pl4",250,-0.5,249.5));
  NStripACh5Arm0Pl4->SetDirectory(0);

  NStripACh5Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripACh5Arm1Pl0","NStripACh5Arm1Pl0",250,-0.5,249.5));
  NStripACh5Arm1Pl0->SetDirectory(0);
  NStripACh5Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripACh5Arm1Pl1","NStripACh5Arm1Pl1",250,-0.5,249.5));
  NStripACh5Arm1Pl1->SetDirectory(0);
  NStripACh5Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripACh5Arm1Pl2","NStripACh5Arm1Pl2",250,-0.5,249.5));
  NStripACh5Arm1Pl2->SetDirectory(0);
  NStripACh5Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripACh5Arm1Pl3","NStripACh5Arm1Pl3",250,-0.5,249.5));
  NStripACh5Arm1Pl3->SetDirectory(0);
  NStripACh5Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripACh5Arm1Pl4","NStripACh5Arm1Pl4",250,-0.5,249.5));
  NStripACh5Arm1Pl4->SetDirectory(0);


  NStripBCh5Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripBCh5Arm0Pl0","NStripBCh5Arm0Pl0",250,-0.5,249.5));
  NStripBCh5Arm0Pl0->SetDirectory(0);
  NStripBCh5Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripBCh5Arm0Pl1","NStripBCh5Arm0Pl1",250,-0.5,249.5));
  NStripBCh5Arm0Pl1->SetDirectory(0);
  NStripBCh5Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripBCh5Arm0Pl2","NStripBCh5Arm0Pl2",250,-0.5,249.5));
  NStripBCh5Arm0Pl2->SetDirectory(0);
  NStripBCh5Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripBCh5Arm0Pl3","NStripBCh5Arm0Pl3",250,-0.5,249.5));
  NStripBCh5Arm0Pl3->SetDirectory(0);
  NStripBCh5Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripBCh5Arm0Pl4","NStripBCh5Arm0Pl4",250,-0.5,249.5));
  NStripBCh5Arm0Pl4->SetDirectory(0);

  NStripBCh5Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("NStripBCh5Arm1Pl0","NStripBCh5Arm1Pl0",250,-0.5,249.5));
  NStripBCh5Arm1Pl0->SetDirectory(0);
  NStripBCh5Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("NStripBCh5Arm1Pl1","NStripBCh5Arm1Pl1",250,-0.5,249.5));
  NStripBCh5Arm1Pl1->SetDirectory(0);
  NStripBCh5Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("NStripBCh5Arm1Pl2","NStripBCh5Arm1Pl2",250,-0.5,249.5));
  NStripBCh5Arm1Pl2->SetDirectory(0);
  NStripBCh5Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("NStripBCh5Arm1Pl3","NStripBCh5Arm1Pl3",250,-0.5,249.5));
  NStripBCh5Arm1Pl3->SetDirectory(0);
  NStripBCh5Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("NStripBCh5Arm1Pl4","NStripBCh5Arm1 Pl4",250,-0.5,249.5));
  NStripBCh5Arm1Pl4->SetDirectory(0);
///////////////////


  RecoX_Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Pl4", "RecoXArm1Pl4", 300, -2000.,2000.));
  RecoX_Arm1Pl4->SetDirectory(0);
  RecoX_Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Pl3", "RecoXArm1Pl3", 300, -2000.,2000.));
  RecoX_Arm1Pl3->SetDirectory(0);
  RecoX_Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Pl2", "RecoXArm1Pl2", 300, -2000.,2000.));
  RecoX_Arm1Pl2->SetDirectory(0);
  RecoX_Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Pl1", "RecoXArm1Pl1", 300, -2000.,2000.));
  RecoX_Arm1Pl1->SetDirectory(0);
  RecoX_Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Pl0", "RecoXArm1Pl0", 300, -2000.,2000.));
  RecoX_Arm1Pl0->SetDirectory(0);

  RecoX_Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Pl4", "RecoXArm0Pl4", 300, -2000.,2000.));
  RecoX_Arm0Pl4->SetDirectory(0);
  RecoX_Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Pl3", "RecoXArm0Pl3", 300, -2000.,2000.));
  RecoX_Arm0Pl3->SetDirectory(0);
  RecoX_Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Pl2", "RecoXArm0Pl2", 300, -2000.,2000.));
  RecoX_Arm0Pl2->SetDirectory(0);
  RecoX_Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Pl1", "RecoXArm0Pl1", 300, -2000.,2000.));
  RecoX_Arm0Pl1->SetDirectory(0);
  RecoX_Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Pl0", "RecoXArm0Pl0", 300, -2000.,2000.));
  RecoX_Arm0Pl0->SetDirectory(0);


  RecoY_Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Pl4", "RecoYArm1Pl4", 300, -2000.,2000.));
  RecoY_Arm1Pl4->SetDirectory(0);
  RecoY_Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Pl3", "RecoYArm1Pl3", 300, -2000.,2000.));
  RecoY_Arm1Pl3->SetDirectory(0);
  RecoY_Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Pl2", "RecoYArm1Pl2", 300, -2000.,2000.));
  RecoY_Arm1Pl2->SetDirectory(0);
  RecoY_Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Pl1", "RecoYArm1Pl1", 300, -2000.,2000.));
  RecoY_Arm1Pl1->SetDirectory(0);
  RecoY_Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Pl0", "RecoYArm1Pl0", 300, -2000.,2000.));
  RecoY_Arm1Pl0->SetDirectory(0);

  RecoY_Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Pl4", "RecoYArm0Pl4", 300, -2000.,2000.));
  RecoY_Arm0Pl4->SetDirectory(0);
  RecoY_Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Pl3", "RecoYArm0Pl3", 300, -2000.,2000.));
  RecoY_Arm0Pl3->SetDirectory(0);
  RecoY_Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Pl2", "RecoYArm0Pl2", 300, -2000.,2000.));
  RecoY_Arm0Pl2->SetDirectory(0);
  RecoY_Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Pl1", "RecoYArm0Pl1", 300, -2000.,2000.));
  RecoY_Arm0Pl1->SetDirectory(0);
  RecoY_Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Pl0", "RecoYArm0Pl0", 300, -2000.,2000.));
  RecoY_Arm0Pl0->SetDirectory(0);
///----------------------------------------------

  DeltaX_Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("DeltaXArm1Pl4", "DeltaXArm1Pl4", 300, -100.,100.));
  DeltaX_Arm1Pl4->SetDirectory(0);
  DeltaX_Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("DeltaXArm1Pl3", "DeltaXArm1Pl3", 300, -100.,100.));
  DeltaX_Arm1Pl3->SetDirectory(0);
  DeltaX_Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("DeltaXArm1Pl2", "DeltaXArm1Pl2", 300, -100.,100.));
  DeltaX_Arm1Pl2->SetDirectory(0);
  DeltaX_Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("DeltaXArm1Pl1", "DeltaXArm1Pl1", 300, -100.,100.));
  DeltaX_Arm1Pl1->SetDirectory(0);
  DeltaX_Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("DeltaXArm1Pl0", "DeltaXArm1Pl0", 300, -100.,100.));
  DeltaX_Arm1Pl0->SetDirectory(0);

  DeltaX_Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("DeltaXArm0Pl4", "DeltaXArm0Pl4", 300, -100.,100.));
  DeltaX_Arm0Pl4->SetDirectory(0);
  DeltaX_Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("DeltaXArm0Pl3", "DeltaXArm0Pl3", 300, -100.,100.));
  DeltaX_Arm0Pl3->SetDirectory(0);
  DeltaX_Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("DeltaXArm0Pl2", "DeltaXArm0Pl2", 300, -100.,100.));
  DeltaX_Arm0Pl2->SetDirectory(0);
  DeltaX_Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("DeltaXArm0Pl1", "DeltaXArm0Pl1", 300, -100.,100.));
  DeltaX_Arm0Pl1->SetDirectory(0);
  DeltaX_Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("DeltaXArm0Pl0", "DeltaXArm0Pl0", 300, -100.,100.));
  DeltaX_Arm0Pl0->SetDirectory(0);


  DeltaY_Arm1Pl4 = std::auto_ptr<TH1D>(new TH1D("DeltaYArm1Pl4", "DeltaYArm1Pl4", 300, -100.,100.));
  DeltaY_Arm1Pl4->SetDirectory(0);
  DeltaY_Arm1Pl3 = std::auto_ptr<TH1D>(new TH1D("DeltaYArm1Pl3", "DeltaYArm1Pl3", 300, -100.,100.));
  DeltaY_Arm1Pl3->SetDirectory(0);
  DeltaY_Arm1Pl2 = std::auto_ptr<TH1D>(new TH1D("DeltaYArm1Pl2", "DeltaYArm1Pl2", 300, -100.,100.));
  DeltaY_Arm1Pl2->SetDirectory(0);
  DeltaY_Arm1Pl1 = std::auto_ptr<TH1D>(new TH1D("DeltaYArm1Pl1", "DeltaYArm1Pl1", 300, -100.,100.));
  DeltaY_Arm1Pl1->SetDirectory(0);
  DeltaY_Arm1Pl0 = std::auto_ptr<TH1D>(new TH1D("DeltaYArm1Pl0", "DeltaYArm1Pl0", 300, -100.,100.));
  DeltaY_Arm1Pl0->SetDirectory(0);

  DeltaY_Arm0Pl4 = std::auto_ptr<TH1D>(new TH1D("DeltaYArm0Pl4", "DeltaYArm0Pl4", 300, -100.,100.));
  DeltaY_Arm0Pl4->SetDirectory(0);
  DeltaY_Arm0Pl3 = std::auto_ptr<TH1D>(new TH1D("DeltaYArm0Pl3", "DeltaYArm0Pl3", 300, -100.,100.));
  DeltaY_Arm0Pl3->SetDirectory(0);
  DeltaY_Arm0Pl2 = std::auto_ptr<TH1D>(new TH1D("DeltaYArm0Pl2", "DeltaYArm0Pl2", 300, -100.,100.));
  DeltaY_Arm0Pl2->SetDirectory(0);
  DeltaY_Arm0Pl1 = std::auto_ptr<TH1D>(new TH1D("DeltaYArm0Pl1", "DeltaYArm0Pl1", 300, -100.,100.));
  DeltaY_Arm0Pl1->SetDirectory(0);
  DeltaY_Arm0Pl0 = std::auto_ptr<TH1D>(new TH1D("DeltaYArm0Pl0", "DeltaYArm0Pl0", 300, -100.,100.));
  DeltaY_Arm0Pl0->SetDirectory(0);


RecoLX_Arm0Plane0CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane0CSC0","RecoXArm0Plane0CSC0",300,-1500,1500));
RecoLX_Arm0Plane0CSC0->SetDirectory(0);
RecoLY_Arm0Plane0CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane0CSC0","RecoYArm0Plane0CSC0",300,-1500,1500));
RecoLY_Arm0Plane0CSC0->SetDirectory(0);
NReco_Arm0Plane0CSC0 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane0CSC0","NRecoArm0Plane0CSC0",100,-0.5,99.5));
NReco_Arm0Plane0CSC0->SetDirectory(0);
RecoLX_Arm0Plane0CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane0CSC1","RecoXArm0Plane0CSC1",300,-1500,1500));
RecoLX_Arm0Plane0CSC1->SetDirectory(0);
RecoLY_Arm0Plane0CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane0CSC1","RecoYArm0Plane0CSC1",300,-1500,1500));
RecoLY_Arm0Plane0CSC1->SetDirectory(0);
NReco_Arm0Plane0CSC1 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane0CSC1","NRecoArm0Plane0CSC1",100,-0.5,99.5));
NReco_Arm0Plane0CSC1->SetDirectory(0);
RecoLX_Arm0Plane0CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane0CSC2","RecoXArm0Plane0CSC2",300,-1500,1500));
RecoLX_Arm0Plane0CSC2->SetDirectory(0);
RecoLY_Arm0Plane0CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane0CSC2","RecoYArm0Plane0CSC2",300,-1500,1500));
RecoLY_Arm0Plane0CSC2->SetDirectory(0);
NReco_Arm0Plane0CSC2 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane0CSC2","NRecoArm0Plane0CSC2",100,-0.5,99.5));
NReco_Arm0Plane0CSC2->SetDirectory(0);
RecoLX_Arm0Plane0CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane0CSC3","RecoXArm0Plane0CSC3",300,-1500,1500));
RecoLX_Arm0Plane0CSC3->SetDirectory(0);
RecoLY_Arm0Plane0CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane0CSC3","RecoYArm0Plane0CSC3",300,-1500,1500));
RecoLY_Arm0Plane0CSC3->SetDirectory(0);
NReco_Arm0Plane0CSC3 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane0CSC3","NRecoArm0Plane0CSC3",100,-0.5,99.5));
NReco_Arm0Plane0CSC3->SetDirectory(0);
RecoLX_Arm0Plane0CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane0CSC4","RecoXArm0Plane0CSC4",300,-1500,1500));
RecoLX_Arm0Plane0CSC4->SetDirectory(0);
RecoLY_Arm0Plane0CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane0CSC4","RecoYArm0Plane0CSC4",300,-1500,1500));
RecoLY_Arm0Plane0CSC4->SetDirectory(0);
NReco_Arm0Plane0CSC4 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane0CSC4","NRecoArm0Plane0CSC4",100,-0.5,99.5));
NReco_Arm0Plane0CSC4->SetDirectory(0);
RecoLX_Arm0Plane0CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane0CSC5","RecoXArm0Plane0CSC5",300,-1500,1500));
RecoLX_Arm0Plane0CSC5->SetDirectory(0);
RecoLY_Arm0Plane0CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane0CSC5","RecoYArm0Plane0CSC5",300,-1500,1500));
RecoLY_Arm0Plane0CSC5->SetDirectory(0);
NReco_Arm0Plane0CSC5 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane0CSC5","NRecoArm0Plane0CSC5",100,-0.5,99.5));
NReco_Arm0Plane0CSC5->SetDirectory(0);
RecoLX_Arm0Plane1CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane1CSC0","RecoXArm0Plane1CSC0",300,-1500,1500));
RecoLX_Arm0Plane1CSC0->SetDirectory(0);
RecoLY_Arm0Plane1CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane1CSC0","RecoYArm0Plane1CSC0",300,-1500,1500));
RecoLY_Arm0Plane1CSC0->SetDirectory(0);
NReco_Arm0Plane1CSC0 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane1CSC0","NRecoArm0Plane1CSC0",100,-0.5,99.5));
NReco_Arm0Plane1CSC0->SetDirectory(0);
RecoLX_Arm0Plane1CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane1CSC1","RecoXArm0Plane1CSC1",300,-1500,1500));
RecoLX_Arm0Plane1CSC1->SetDirectory(0);
RecoLY_Arm0Plane1CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane1CSC1","RecoYArm0Plane1CSC1",300,-1500,1500));
RecoLY_Arm0Plane1CSC1->SetDirectory(0);
NReco_Arm0Plane1CSC1 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane1CSC1","NRecoArm0Plane1CSC1",100,-0.5,99.5));
NReco_Arm0Plane1CSC1->SetDirectory(0);
RecoLX_Arm0Plane1CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane1CSC2","RecoXArm0Plane1CSC2",300,-1500,1500));
RecoLX_Arm0Plane1CSC2->SetDirectory(0);
RecoLY_Arm0Plane1CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane1CSC2","RecoYArm0Plane1CSC2",300,-1500,1500));
RecoLY_Arm0Plane1CSC2->SetDirectory(0);
NReco_Arm0Plane1CSC2 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane1CSC2","NRecoArm0Plane1CSC2",100,-0.5,99.5));
NReco_Arm0Plane1CSC2->SetDirectory(0);
RecoLX_Arm0Plane1CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane1CSC3","RecoXArm0Plane1CSC3",300,-1500,1500));
RecoLX_Arm0Plane1CSC3->SetDirectory(0);
RecoLY_Arm0Plane1CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane1CSC3","RecoYArm0Plane1CSC3",300,-1500,1500));
RecoLY_Arm0Plane1CSC3->SetDirectory(0);
NReco_Arm0Plane1CSC3 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane1CSC3","NRecoArm0Plane1CSC3",100,-0.5,99.5));
NReco_Arm0Plane1CSC3->SetDirectory(0);
RecoLX_Arm0Plane1CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane1CSC4","RecoXArm0Plane1CSC4",300,-1500,1500));
RecoLX_Arm0Plane1CSC4->SetDirectory(0);
RecoLY_Arm0Plane1CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane1CSC4","RecoYArm0Plane1CSC4",300,-1500,1500));
RecoLY_Arm0Plane1CSC4->SetDirectory(0);
NReco_Arm0Plane1CSC4 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane1CSC4","NRecoArm0Plane1CSC4",100,-0.5,99.5));
NReco_Arm0Plane1CSC4->SetDirectory(0);
RecoLX_Arm0Plane1CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane1CSC5","RecoXArm0Plane1CSC5",300,-1500,1500));
RecoLX_Arm0Plane1CSC5->SetDirectory(0);
RecoLY_Arm0Plane1CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane1CSC5","RecoYArm0Plane1CSC5",300,-1500,1500));
RecoLY_Arm0Plane1CSC5->SetDirectory(0);
NReco_Arm0Plane1CSC5 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane1CSC5","NRecoArm0Plane1CSC5",100,-0.5,99.5));
NReco_Arm0Plane1CSC5->SetDirectory(0);
RecoLX_Arm0Plane2CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane2CSC0","RecoXArm0Plane2CSC0",300,-1500,1500));
RecoLX_Arm0Plane2CSC0->SetDirectory(0);
RecoLY_Arm0Plane2CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane2CSC0","RecoYArm0Plane2CSC0",300,-1500,1500));
RecoLY_Arm0Plane2CSC0->SetDirectory(0);
NReco_Arm0Plane2CSC0 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane2CSC0","NRecoArm0Plane2CSC0",100,-0.5,99.5));
NReco_Arm0Plane2CSC0->SetDirectory(0);
RecoLX_Arm0Plane2CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane2CSC1","RecoXArm0Plane2CSC1",300,-1500,1500));
RecoLX_Arm0Plane2CSC1->SetDirectory(0);
RecoLY_Arm0Plane2CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane2CSC1","RecoYArm0Plane2CSC1",300,-1500,1500));
RecoLY_Arm0Plane2CSC1->SetDirectory(0);
NReco_Arm0Plane2CSC1 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane2CSC1","NRecoArm0Plane2CSC1",100,-0.5,99.5));
NReco_Arm0Plane2CSC1->SetDirectory(0);
RecoLX_Arm0Plane2CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane2CSC2","RecoXArm0Plane2CSC2",300,-1500,1500));
RecoLX_Arm0Plane2CSC2->SetDirectory(0);
RecoLY_Arm0Plane2CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane2CSC2","RecoYArm0Plane2CSC2",300,-1500,1500));
RecoLY_Arm0Plane2CSC2->SetDirectory(0);
NReco_Arm0Plane2CSC2 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane2CSC2","NRecoArm0Plane2CSC2",100,-0.5,99.5));
NReco_Arm0Plane2CSC2->SetDirectory(0);
RecoLX_Arm0Plane2CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane2CSC3","RecoXArm0Plane2CSC3",300,-1500,1500));
RecoLX_Arm0Plane2CSC3->SetDirectory(0);
RecoLY_Arm0Plane2CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane2CSC3","RecoYArm0Plane2CSC3",300,-1500,1500));
RecoLY_Arm0Plane2CSC3->SetDirectory(0);
NReco_Arm0Plane2CSC3 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane2CSC3","NRecoArm0Plane2CSC3",100,-0.5,99.5));
NReco_Arm0Plane2CSC3->SetDirectory(0);
RecoLX_Arm0Plane2CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane2CSC4","RecoXArm0Plane2CSC4",300,-1500,1500));
RecoLX_Arm0Plane2CSC4->SetDirectory(0);
RecoLY_Arm0Plane2CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane2CSC4","RecoYArm0Plane2CSC4",300,-1500,1500));
RecoLY_Arm0Plane2CSC4->SetDirectory(0);
NReco_Arm0Plane2CSC4 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane2CSC4","NRecoArm0Plane2CSC4",100,-0.5,99.5));
NReco_Arm0Plane2CSC4->SetDirectory(0);
RecoLX_Arm0Plane2CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane2CSC5","RecoXArm0Plane2CSC5",300,-1500,1500));
RecoLX_Arm0Plane2CSC5->SetDirectory(0);
RecoLY_Arm0Plane2CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane2CSC5","RecoYArm0Plane2CSC5",300,-1500,1500));
RecoLY_Arm0Plane2CSC5->SetDirectory(0);
NReco_Arm0Plane2CSC5 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane2CSC5","NRecoArm0Plane2CSC5",100,-0.5,99.5));
NReco_Arm0Plane2CSC5->SetDirectory(0);
RecoLX_Arm0Plane3CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane3CSC0","RecoXArm0Plane3CSC0",300,-1500,1500));
RecoLX_Arm0Plane3CSC0->SetDirectory(0);
RecoLY_Arm0Plane3CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane3CSC0","RecoYArm0Plane3CSC0",300,-1500,1500));
RecoLY_Arm0Plane3CSC0->SetDirectory(0);
NReco_Arm0Plane3CSC0 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane3CSC0","NRecoArm0Plane3CSC0",100,-0.5,99.5));
NReco_Arm0Plane3CSC0->SetDirectory(0);
RecoLX_Arm0Plane3CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane3CSC1","RecoXArm0Plane3CSC1",300,-1500,1500));
RecoLX_Arm0Plane3CSC1->SetDirectory(0);
RecoLY_Arm0Plane3CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane3CSC1","RecoYArm0Plane3CSC1",300,-1500,1500));
RecoLY_Arm0Plane3CSC1->SetDirectory(0);
NReco_Arm0Plane3CSC1 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane3CSC1","NRecoArm0Plane3CSC1",100,-0.5,99.5));
NReco_Arm0Plane3CSC1->SetDirectory(0);
RecoLX_Arm0Plane3CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane3CSC2","RecoXArm0Plane3CSC2",300,-1500,1500));
RecoLX_Arm0Plane3CSC2->SetDirectory(0);
RecoLY_Arm0Plane3CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane3CSC2","RecoYArm0Plane3CSC2",300,-1500,1500));
RecoLY_Arm0Plane3CSC2->SetDirectory(0);
NReco_Arm0Plane3CSC2 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane3CSC2","NRecoArm0Plane3CSC2",100,-0.5,99.5));
NReco_Arm0Plane3CSC2->SetDirectory(0);
RecoLX_Arm0Plane3CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane3CSC3","RecoXArm0Plane3CSC3",300,-1500,1500));
RecoLX_Arm0Plane3CSC3->SetDirectory(0);
RecoLY_Arm0Plane3CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane3CSC3","RecoYArm0Plane3CSC3",300,-1500,1500));
RecoLY_Arm0Plane3CSC3->SetDirectory(0);
NReco_Arm0Plane3CSC3 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane3CSC3","NRecoArm0Plane3CSC3",100,-0.5,99.5));
NReco_Arm0Plane3CSC3->SetDirectory(0);
RecoLX_Arm0Plane3CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane3CSC4","RecoXArm0Plane3CSC4",300,-1500,1500));
RecoLX_Arm0Plane3CSC4->SetDirectory(0);
RecoLY_Arm0Plane3CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane3CSC4","RecoYArm0Plane3CSC4",300,-1500,1500));
RecoLY_Arm0Plane3CSC4->SetDirectory(0);
NReco_Arm0Plane3CSC4 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane3CSC4","NRecoArm0Plane3CSC4",100,-0.5,99.5));
NReco_Arm0Plane3CSC4->SetDirectory(0);
RecoLX_Arm0Plane3CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane3CSC5","RecoXArm0Plane3CSC5",300,-1500,1500));
RecoLX_Arm0Plane3CSC5->SetDirectory(0);
RecoLY_Arm0Plane3CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane3CSC5","RecoYArm0Plane3CSC5",300,-1500,1500));
RecoLY_Arm0Plane3CSC5->SetDirectory(0);
NReco_Arm0Plane3CSC5 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane3CSC5","NRecoArm0Plane3CSC5",100,-0.5,99.5));
NReco_Arm0Plane3CSC5->SetDirectory(0);
RecoLX_Arm0Plane4CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane4CSC0","RecoXArm0Plane4CSC0",300,-1500,1500));
RecoLX_Arm0Plane4CSC0->SetDirectory(0);
RecoLY_Arm0Plane4CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane4CSC0","RecoYArm0Plane4CSC0",300,-1500,1500));
RecoLY_Arm0Plane4CSC0->SetDirectory(0);
NReco_Arm0Plane4CSC0 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane4CSC0","NRecoArm0Plane4CSC0",100,-0.5,99.5));
NReco_Arm0Plane4CSC0->SetDirectory(0);
RecoLX_Arm0Plane4CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane4CSC1","RecoXArm0Plane4CSC1",300,-1500,1500));
RecoLX_Arm0Plane4CSC1->SetDirectory(0);
RecoLY_Arm0Plane4CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane4CSC1","RecoYArm0Plane4CSC1",300,-1500,1500));
RecoLY_Arm0Plane4CSC1->SetDirectory(0);
NReco_Arm0Plane4CSC1 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane4CSC1","NRecoArm0Plane4CSC1",100,-0.5,99.5));
NReco_Arm0Plane4CSC1->SetDirectory(0);
RecoLX_Arm0Plane4CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane4CSC2","RecoXArm0Plane4CSC2",300,-1500,1500));
RecoLX_Arm0Plane4CSC2->SetDirectory(0);
RecoLY_Arm0Plane4CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane4CSC2","RecoYArm0Plane4CSC2",300,-1500,1500));
RecoLY_Arm0Plane4CSC2->SetDirectory(0);
NReco_Arm0Plane4CSC2 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane4CSC2","NRecoArm0Plane4CSC2",100,-0.5,99.5));
NReco_Arm0Plane4CSC2->SetDirectory(0);
RecoLX_Arm0Plane4CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane4CSC3","RecoXArm0Plane4CSC3",300,-1500,1500));
RecoLX_Arm0Plane4CSC3->SetDirectory(0);
RecoLY_Arm0Plane4CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane4CSC3","RecoYArm0Plane4CSC3",300,-1500,1500));
RecoLY_Arm0Plane4CSC3->SetDirectory(0);
NReco_Arm0Plane4CSC3 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane4CSC3","NRecoArm0Plane4CSC3",100,-0.5,99.5));
NReco_Arm0Plane4CSC3->SetDirectory(0);
RecoLX_Arm0Plane4CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane4CSC4","RecoXArm0Plane4CSC4",300,-1500,1500));
RecoLX_Arm0Plane4CSC4->SetDirectory(0);
RecoLY_Arm0Plane4CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane4CSC4","RecoYArm0Plane4CSC4",300,-1500,1500));
RecoLY_Arm0Plane4CSC4->SetDirectory(0);
NReco_Arm0Plane4CSC4 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane4CSC4","NRecoArm0Plane4CSC4",100,-0.5,99.5));
NReco_Arm0Plane4CSC4->SetDirectory(0);
RecoLX_Arm0Plane4CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoXArm0Plane4CSC5","RecoXArm0Plane4CSC5",300,-1500,1500));
RecoLX_Arm0Plane4CSC5->SetDirectory(0);
RecoLY_Arm0Plane4CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoYArm0Plane4CSC5","RecoYArm0Plane4CSC5",300,-1500,1500));
RecoLY_Arm0Plane4CSC5->SetDirectory(0);
NReco_Arm0Plane4CSC5 = std::auto_ptr<TH1D>(new TH1D("NRecoArm0Plane4CSC5","NRecoArm0Plane4CSC5",100,-0.5,99.5));
NReco_Arm0Plane4CSC5->SetDirectory(0);
RecoLX_Arm1Plane0CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane0CSC0","RecoXArm1Plane0CSC0",300,-1500,1500));
RecoLX_Arm1Plane0CSC0->SetDirectory(0);
RecoLY_Arm1Plane0CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane0CSC0","RecoYArm1Plane0CSC0",300,-1500,1500));
RecoLY_Arm1Plane0CSC0->SetDirectory(0);
NReco_Arm1Plane0CSC0 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane0CSC0","NRecoArm1Plane0CSC0",100,-0.5,99.5));
NReco_Arm1Plane0CSC0->SetDirectory(0);
RecoLX_Arm1Plane0CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane0CSC1","RecoXArm1Plane0CSC1",300,-1500,1500));
RecoLX_Arm1Plane0CSC1->SetDirectory(0);
RecoLY_Arm1Plane0CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane0CSC1","RecoYArm1Plane0CSC1",300,-1500,1500));
RecoLY_Arm1Plane0CSC1->SetDirectory(0);
NReco_Arm1Plane0CSC1 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane0CSC1","NRecoArm1Plane0CSC1",100,-0.5,99.5));
NReco_Arm1Plane0CSC1->SetDirectory(0);
RecoLX_Arm1Plane0CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane0CSC2","RecoXArm1Plane0CSC2",300,-1500,1500));
RecoLX_Arm1Plane0CSC2->SetDirectory(0);
RecoLY_Arm1Plane0CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane0CSC2","RecoYArm1Plane0CSC2",300,-1500,1500));
RecoLY_Arm1Plane0CSC2->SetDirectory(0);
NReco_Arm1Plane0CSC2 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane0CSC2","NRecoArm1Plane0CSC2",100,-0.5,99.5));
NReco_Arm1Plane0CSC2->SetDirectory(0);
RecoLX_Arm1Plane0CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane0CSC3","RecoXArm1Plane0CSC3",300,-1500,1500));
RecoLX_Arm1Plane0CSC3->SetDirectory(0);
RecoLY_Arm1Plane0CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane0CSC3","RecoYArm1Plane0CSC3",300,-1500,1500));
RecoLY_Arm1Plane0CSC3->SetDirectory(0);
NReco_Arm1Plane0CSC3 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane0CSC3","NRecoArm1Plane0CSC3",100,-0.5,99.5));
NReco_Arm1Plane0CSC3->SetDirectory(0);
RecoLX_Arm1Plane0CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane0CSC4","RecoXArm1Plane0CSC4",300,-1500,1500));
RecoLX_Arm1Plane0CSC4->SetDirectory(0);
RecoLY_Arm1Plane0CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane0CSC4","RecoYArm1Plane0CSC4",300,-1500,1500));
RecoLY_Arm1Plane0CSC4->SetDirectory(0);
NReco_Arm1Plane0CSC4 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane0CSC4","NRecoArm1Plane0CSC4",100,-0.5,99.5));
NReco_Arm1Plane0CSC4->SetDirectory(0);
RecoLX_Arm1Plane0CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane0CSC5","RecoXArm1Plane0CSC5",300,-1500,1500));
RecoLX_Arm1Plane0CSC5->SetDirectory(0);
RecoLY_Arm1Plane0CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane0CSC5","RecoYArm1Plane0CSC5",300,-1500,1500));
RecoLY_Arm1Plane0CSC5->SetDirectory(0);
NReco_Arm1Plane0CSC5 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane0CSC5","NRecoArm1Plane0CSC5",100,-0.5,99.5));
NReco_Arm1Plane0CSC5->SetDirectory(0);
RecoLX_Arm1Plane1CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane1CSC0","RecoXArm1Plane1CSC0",300,-1500,1500));
RecoLX_Arm1Plane1CSC0->SetDirectory(0);
RecoLY_Arm1Plane1CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane1CSC0","RecoYArm1Plane1CSC0",300,-1500,1500));
RecoLY_Arm1Plane1CSC0->SetDirectory(0);
NReco_Arm1Plane1CSC0 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane1CSC0","NRecoArm1Plane1CSC0",100,-0.5,99.5));
NReco_Arm1Plane1CSC0->SetDirectory(0);
RecoLX_Arm1Plane1CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane1CSC1","RecoXArm1Plane1CSC1",300,-1500,1500));
RecoLX_Arm1Plane1CSC1->SetDirectory(0);
RecoLY_Arm1Plane1CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane1CSC1","RecoYArm1Plane1CSC1",300,-1500,1500));
RecoLY_Arm1Plane1CSC1->SetDirectory(0);
NReco_Arm1Plane1CSC1 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane1CSC1","NRecoArm1Plane1CSC1",100,-0.5,99.5));
NReco_Arm1Plane1CSC1->SetDirectory(0);
RecoLX_Arm1Plane1CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane1CSC2","RecoXArm1Plane1CSC2",300,-1500,1500));
RecoLX_Arm1Plane1CSC2->SetDirectory(0);
RecoLY_Arm1Plane1CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane1CSC2","RecoYArm1Plane1CSC2",300,-1500,1500));
RecoLY_Arm1Plane1CSC2->SetDirectory(0);
NReco_Arm1Plane1CSC2 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane1CSC2","NRecoArm1Plane1CSC2",100,-0.5,99.5));
NReco_Arm1Plane1CSC2->SetDirectory(0);
RecoLX_Arm1Plane1CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane1CSC3","RecoXArm1Plane1CSC3",300,-1500,1500));
RecoLX_Arm1Plane1CSC3->SetDirectory(0);
RecoLY_Arm1Plane1CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane1CSC3","RecoYArm1Plane1CSC3",300,-1500,1500));
RecoLY_Arm1Plane1CSC3->SetDirectory(0);
NReco_Arm1Plane1CSC3 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane1CSC3","NRecoArm1Plane1CSC3",100,-0.5,99.5));
NReco_Arm1Plane1CSC3->SetDirectory(0);
RecoLX_Arm1Plane1CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane1CSC4","RecoXArm1Plane1CSC4",300,-1500,1500));
RecoLX_Arm1Plane1CSC4->SetDirectory(0);
RecoLY_Arm1Plane1CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane1CSC4","RecoYArm1Plane1CSC4",300,-1500,1500));
RecoLY_Arm1Plane1CSC4->SetDirectory(0);
NReco_Arm1Plane1CSC4 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane1CSC4","NRecoArm1Plane1CSC4",100,-0.5,99.5));
NReco_Arm1Plane1CSC4->SetDirectory(0);
RecoLX_Arm1Plane1CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane1CSC5","RecoXArm1Plane1CSC5",300,-1500,1500));
RecoLX_Arm1Plane1CSC5->SetDirectory(0);
RecoLY_Arm1Plane1CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane1CSC5","RecoYArm1Plane1CSC5",300,-1500,1500));
RecoLY_Arm1Plane1CSC5->SetDirectory(0);
NReco_Arm1Plane1CSC5 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane1CSC5","NRecoArm1Plane1CSC5",100,-0.5,99.5));
NReco_Arm1Plane1CSC5->SetDirectory(0);
RecoLX_Arm1Plane2CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane2CSC0","RecoXArm1Plane2CSC0",300,-1500,1500));
RecoLX_Arm1Plane2CSC0->SetDirectory(0);
RecoLY_Arm1Plane2CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane2CSC0","RecoYArm1Plane2CSC0",300,-1500,1500));
RecoLY_Arm1Plane2CSC0->SetDirectory(0);
NReco_Arm1Plane2CSC0 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane2CSC0","NRecoArm1Plane2CSC0",100,-0.5,99.5));
NReco_Arm1Plane2CSC0->SetDirectory(0);
RecoLX_Arm1Plane2CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane2CSC1","RecoXArm1Plane2CSC1",300,-1500,1500));
RecoLX_Arm1Plane2CSC1->SetDirectory(0);
RecoLY_Arm1Plane2CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane2CSC1","RecoYArm1Plane2CSC1",300,-1500,1500));
RecoLY_Arm1Plane2CSC1->SetDirectory(0);
NReco_Arm1Plane2CSC1 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane2CSC1","NRecoArm1Plane2CSC1",100,-0.5,99.5));
NReco_Arm1Plane2CSC1->SetDirectory(0);
RecoLX_Arm1Plane2CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane2CSC2","RecoXArm1Plane2CSC2",300,-1500,1500));
RecoLX_Arm1Plane2CSC2->SetDirectory(0);
RecoLY_Arm1Plane2CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane2CSC2","RecoYArm1Plane2CSC2",300,-1500,1500));
RecoLY_Arm1Plane2CSC2->SetDirectory(0);
NReco_Arm1Plane2CSC2 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane2CSC2","NRecoArm1Plane2CSC2",100,-0.5,99.5));
NReco_Arm1Plane2CSC2->SetDirectory(0);
RecoLX_Arm1Plane2CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane2CSC3","RecoXArm1Plane2CSC3",300,-1500,1500));
RecoLX_Arm1Plane2CSC3->SetDirectory(0);
RecoLY_Arm1Plane2CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane2CSC3","RecoYArm1Plane2CSC3",300,-1500,1500));
RecoLY_Arm1Plane2CSC3->SetDirectory(0);
NReco_Arm1Plane2CSC3 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane2CSC3","NRecoArm1Plane2CSC3",100,-0.5,99.5));
NReco_Arm1Plane2CSC3->SetDirectory(0);
RecoLX_Arm1Plane2CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane2CSC4","RecoXArm1Plane2CSC4",300,-1500,1500));
RecoLX_Arm1Plane2CSC4->SetDirectory(0);
RecoLY_Arm1Plane2CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane2CSC4","RecoYArm1Plane2CSC4",300,-1500,1500));
RecoLY_Arm1Plane2CSC4->SetDirectory(0);
NReco_Arm1Plane2CSC4 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane2CSC4","NRecoArm1Plane2CSC4",100,-0.5,99.5));
NReco_Arm1Plane2CSC4->SetDirectory(0);
RecoLX_Arm1Plane2CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane2CSC5","RecoXArm1Plane2CSC5",300,-1500,1500));
RecoLX_Arm1Plane2CSC5->SetDirectory(0);
RecoLY_Arm1Plane2CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane2CSC5","RecoYArm1Plane2CSC5",300,-1500,1500));
RecoLY_Arm1Plane2CSC5->SetDirectory(0);
NReco_Arm1Plane2CSC5 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane2CSC5","NRecoArm1Plane2CSC5",100,-0.5,99.5));
NReco_Arm1Plane2CSC5->SetDirectory(0);
RecoLX_Arm1Plane3CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane3CSC0","RecoXArm1Plane3CSC0",300,-1500,1500));
RecoLX_Arm1Plane3CSC0->SetDirectory(0);
RecoLY_Arm1Plane3CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane3CSC0","RecoYArm1Plane3CSC0",300,-1500,1500));
RecoLY_Arm1Plane3CSC0->SetDirectory(0);
NReco_Arm1Plane3CSC0 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane3CSC0","NRecoArm1Plane3CSC0",100,-0.5,99.5));
NReco_Arm1Plane3CSC0->SetDirectory(0);
RecoLX_Arm1Plane3CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane3CSC1","RecoXArm1Plane3CSC1",300,-1500,1500));
RecoLX_Arm1Plane3CSC1->SetDirectory(0);
RecoLY_Arm1Plane3CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane3CSC1","RecoYArm1Plane3CSC1",300,-1500,1500));
RecoLY_Arm1Plane3CSC1->SetDirectory(0);
NReco_Arm1Plane3CSC1 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane3CSC1","NRecoArm1Plane3CSC1",100,-0.5,99.5));
NReco_Arm1Plane3CSC1->SetDirectory(0);
RecoLX_Arm1Plane3CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane3CSC2","RecoXArm1Plane3CSC2",300,-1500,1500));
RecoLX_Arm1Plane3CSC2->SetDirectory(0);
RecoLY_Arm1Plane3CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane3CSC2","RecoYArm1Plane3CSC2",300,-1500,1500));
RecoLY_Arm1Plane3CSC2->SetDirectory(0);
NReco_Arm1Plane3CSC2 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane3CSC2","NRecoArm1Plane3CSC2",100,-0.5,99.5));
NReco_Arm1Plane3CSC2->SetDirectory(0);
RecoLX_Arm1Plane3CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane3CSC3","RecoXArm1Plane3CSC3",300,-1500,1500));
RecoLX_Arm1Plane3CSC3->SetDirectory(0);
RecoLY_Arm1Plane3CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane3CSC3","RecoYArm1Plane3CSC3",300,-1500,1500));
RecoLY_Arm1Plane3CSC3->SetDirectory(0);
NReco_Arm1Plane3CSC3 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane3CSC3","NRecoArm1Plane3CSC3",100,-0.5,99.5));
NReco_Arm1Plane3CSC3->SetDirectory(0);
RecoLX_Arm1Plane3CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane3CSC4","RecoXArm1Plane3CSC4",300,-1500,1500));
RecoLX_Arm1Plane3CSC4->SetDirectory(0);
RecoLY_Arm1Plane3CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane3CSC4","RecoYArm1Plane3CSC4",300,-1500,1500));
RecoLY_Arm1Plane3CSC4->SetDirectory(0);
NReco_Arm1Plane3CSC4 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane3CSC4","NRecoArm1Plane3CSC4",100,-0.5,99.5));
NReco_Arm1Plane3CSC4->SetDirectory(0);
RecoLX_Arm1Plane3CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane3CSC5","RecoXArm1Plane3CSC5",300,-1500,1500));
RecoLX_Arm1Plane3CSC5->SetDirectory(0);
RecoLY_Arm1Plane3CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane3CSC5","RecoYArm1Plane3CSC5",300,-1500,1500));
RecoLY_Arm1Plane3CSC5->SetDirectory(0);
NReco_Arm1Plane3CSC5 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane3CSC5","NRecoArm1Plane3CSC5",100,-0.5,99.5));
NReco_Arm1Plane3CSC5->SetDirectory(0);
RecoLX_Arm1Plane4CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane4CSC0","RecoXArm1Plane4CSC0",300,-1500,1500));
RecoLX_Arm1Plane4CSC0->SetDirectory(0);
RecoLY_Arm1Plane4CSC0 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane4CSC0","RecoYArm1Plane4CSC0",300,-1500,1500));
RecoLY_Arm1Plane4CSC0->SetDirectory(0);
NReco_Arm1Plane4CSC0 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane4CSC0","NRecoArm1Plane4CSC0",100,-0.5,99.5));
NReco_Arm1Plane4CSC0->SetDirectory(0);
RecoLX_Arm1Plane4CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane4CSC1","RecoXArm1Plane4CSC1",300,-1500,1500));
RecoLX_Arm1Plane4CSC1->SetDirectory(0);
RecoLY_Arm1Plane4CSC1 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane4CSC1","RecoYArm1Plane4CSC1",300,-1500,1500));
RecoLY_Arm1Plane4CSC1->SetDirectory(0);
NReco_Arm1Plane4CSC1 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane4CSC1","NRecoArm1Plane4CSC1",100,-0.5,99.5));
NReco_Arm1Plane4CSC1->SetDirectory(0);
RecoLX_Arm1Plane4CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane4CSC2","RecoXArm1Plane4CSC2",300,-1500,1500));
RecoLX_Arm1Plane4CSC2->SetDirectory(0);
RecoLY_Arm1Plane4CSC2 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane4CSC2","RecoYArm1Plane4CSC2",300,-1500,1500));
RecoLY_Arm1Plane4CSC2->SetDirectory(0);
NReco_Arm1Plane4CSC2 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane4CSC2","NRecoArm1Plane4CSC2",100,-0.5,99.5));
NReco_Arm1Plane4CSC2->SetDirectory(0);
RecoLX_Arm1Plane4CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane4CSC3","RecoXArm1Plane4CSC3",300,-1500,1500));
RecoLX_Arm1Plane4CSC3->SetDirectory(0);
RecoLY_Arm1Plane4CSC3 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane4CSC3","RecoYArm1Plane4CSC3",300,-1500,1500));
RecoLY_Arm1Plane4CSC3->SetDirectory(0);
 NReco_Arm1Plane4CSC3 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane4CSC3","NRecoArm1Plane4CSC3",100,-0.5,99.5));
NReco_Arm1Plane4CSC3->SetDirectory(0);
RecoLX_Arm1Plane4CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane4CSC4","RecoXArm1Plane4CSC4",300,-1500,1500));
RecoLX_Arm1Plane4CSC4->SetDirectory(0);
RecoLY_Arm1Plane4CSC4 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane4CSC4","RecoYArm1Plane4CSC4",300,-1500,1500));
RecoLY_Arm1Plane4CSC4->SetDirectory(0);
NReco_Arm1Plane4CSC4 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane4CSC4","NRecoArm1Plane4CSC4",100,-0.5,99.5));
NReco_Arm1Plane4CSC4->SetDirectory(0);
RecoLX_Arm1Plane4CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoXArm1Plane4CSC5","RecoXArm1Plane4CSC5",300,-1500,1500));
RecoLX_Arm1Plane4CSC5->SetDirectory(0);
RecoLY_Arm1Plane4CSC5 = std::auto_ptr<TH1D>(new TH1D("RecoYArm1Plane4CSC5","RecoYArm1Plane4CSC5",300,-1500,1500));
RecoLY_Arm1Plane4CSC5->SetDirectory(0);
NReco_Arm1Plane4CSC5 = std::auto_ptr<TH1D>(new TH1D("NRecoArm1Plane4CSC5","NRecoArm1Plane4CSC5",100,-0.5,99.5));
NReco_Arm1Plane4CSC5->SetDirectory(0);








  hAllEtaRec = std::auto_ptr<TH1D>(new TH1D("TrackEtaRec","TrackEtaRec",480,-15,15));
  hAllEtaRec->SetDirectory(0);
  hAllPhiRec = std::auto_ptr<TH1D>(new TH1D("TrackPhiRec","TrackPhiRec",50,0,6.3));
  hAllPhiRec->SetDirectory(0);

  hChiSquaredOverN = std::auto_ptr<TH1D>(new TH1D("ChiSquaredOverN","ChiSquaredOverN",100,0,30));
  hChiSquaredOverN->SetDirectory(0);

  hRoadSize = std::auto_ptr<TH1D>(new TH1D("RoadSize","RoadSize",100,0,100));
  hRoadSize->SetDirectory(0);

  hRoadHitSigmaX = std::auto_ptr<TH1D>(new TH1D("RoadHitSigmaX","RoadHitSigmaX",100,0,400));
  hRoadHitSigmaX->SetDirectory(0);

  hRoadHitSigmaY = std::auto_ptr<TH1D>(new TH1D("RoadHitSigmaY","RoadHitSigmaY",100,0,400));
  hRoadHitSigmaY->SetDirectory(0);

  hRoadNumber = std::auto_ptr<TH1D>(new TH1D("RoadNumber","RoadNumber",100,0,100));
  hRoadNumber->SetDirectory(0);

  hChiSquaredXOverN = std::auto_ptr<TH1D>(new TH1D("ChiSquaredXOverN","ChiSquaredXOverN",100,0,30));
  hChiSquaredXOverN->SetDirectory(0);

  hChiSquaredYOverN = std::auto_ptr<TH1D>(new TH1D("ChiSquaredYOverN","ChiSquaredYOverN",100,0,30));
  hChiSquaredYOverN->SetDirectory(0);

  hChiSquaredProb =  std::auto_ptr<TH1D>(new TH1D("ChiSquaredProb","ChiSquaredProb",100,0,1));
  hChiSquaredProb->SetDirectory(0);

  hChiSquaredXProb =  std::auto_ptr<TH1D>(new TH1D("ChiSquaredXProb","ChiSquaredXProb",100,0,1));
  hChiSquaredXProb->SetDirectory(0);

  hChiSquaredYProb =  std::auto_ptr<TH1D>(new TH1D("ChiSquaredYProb","ChiSquaredYProb",100,0,1));
  hChiSquaredYProb->SetDirectory(0);

  ilMomento = std::auto_ptr<TH1D>(new TH1D("SimMomentum","SimMomentum",100,0.,100));
  ilMomento->SetDirectory(0);
  ilMomentoLog = std::auto_ptr<TH1D>(new TH1D("SimMomentumLog","SimMomentumLog",100,-10.,4));
  ilMomentoLog->SetDirectory(0);
  lEnergia = std::auto_ptr<TH1D>(new TH1D("SimEnergy loss","SimEnergy loss",100,0.,1e-05));
  lEnergia->SetDirectory(0);
  ilTipo = std::auto_ptr<TH1D>(new TH1D("SimType","SimType",2000, -1000, 1000));
  ilTipo->SetDirectory(0);
  ilDetector = std::auto_ptr<TH1D>(new TH1D("SimUnitID", "SimUnitID", 100, 1910e06, 1940e06));
  ilDetector->SetDirectory(0);
}

void
T1Validation::writeHistograms() {
  lEta->Write();
  ilPT->Write();
//  hClusterWidth->Write();
  hClDiffBetweenLayers ->Write();

  hClRMSSex->Write();
  ilWireCh0Arm0Pl0->Write();
  ilWireCh0Arm0Pl1->Write();
  ilWireCh0Arm0Pl2->Write();
  ilWireCh0Arm0Pl3->Write();
  ilWireCh0Arm0Pl4->Write();
  ilWireCh0Arm1Pl0->Write();
  ilWireCh0Arm1Pl1->Write();
  ilWireCh0Arm1Pl2->Write();
  ilWireCh0Arm1Pl3->Write();
  ilWireCh0Arm1Pl4->Write();
  ilClWidthCh0Arm0Pl0->Write();
  ilClWidthCh0Arm0Pl1->Write();
  ilClWidthCh0Arm0Pl2->Write();
  ilClWidthCh0Arm0Pl3->Write();
  ilClWidthCh0Arm0Pl4->Write();
  ilClWidthCh0Arm1Pl0->Write();
  ilClWidthCh0Arm1Pl1->Write();
  ilClWidthCh0Arm1Pl2->Write();
  ilClWidthCh0Arm1Pl3->Write();
  ilClWidthCh0Arm1Pl4->Write();
  laStripACh0Arm0Pl0->Write();
  laStripACh0Arm0Pl1->Write();
  laStripACh0Arm0Pl2->Write();
  laStripACh0Arm0Pl3->Write();
  laStripACh0Arm0Pl4->Write();
  laStripACh0Arm1Pl0->Write();
  laStripACh0Arm1Pl1->Write();
  laStripACh0Arm1Pl2->Write();
  laStripACh0Arm1Pl3->Write();
  laStripACh0Arm1Pl4->Write();
  laStripBCh0Arm0Pl0->Write();
  laStripBCh0Arm0Pl1->Write();
  laStripBCh0Arm0Pl2->Write();
  laStripBCh0Arm0Pl3->Write();
  laStripBCh0Arm0Pl4->Write();
  laStripBCh0Arm1Pl0->Write();
  laStripBCh0Arm1Pl1->Write();
  laStripBCh0Arm1Pl2->Write();
  laStripBCh0Arm1Pl3->Write();
  laStripBCh0Arm1Pl4->Write();
  ilWireCh1Arm0Pl0->Write();
  ilWireCh1Arm0Pl1->Write();
  ilWireCh1Arm0Pl2->Write();
  ilWireCh1Arm0Pl3->Write();
  ilWireCh1Arm0Pl4->Write();
  ilWireCh1Arm1Pl0->Write();
  ilWireCh1Arm1Pl1->Write();
  ilWireCh1Arm1Pl2->Write();
  ilWireCh1Arm1Pl3->Write();
  ilWireCh1Arm1Pl4->Write();
  ilClWidthCh1Arm0Pl0->Write();
  ilClWidthCh1Arm0Pl1->Write();
  ilClWidthCh1Arm0Pl2->Write();
  ilClWidthCh1Arm0Pl3->Write();
  ilClWidthCh1Arm0Pl4->Write();
  ilClWidthCh1Arm1Pl0->Write();
  ilClWidthCh1Arm1Pl1->Write();
  ilClWidthCh1Arm1Pl2->Write();
  ilClWidthCh1Arm1Pl3->Write();
  ilClWidthCh1Arm1Pl4->Write();
  laStripACh1Arm0Pl0->Write();
  laStripACh1Arm0Pl1->Write();
  laStripACh1Arm0Pl2->Write();
  laStripACh1Arm0Pl3->Write();
  laStripACh1Arm0Pl4->Write();
  laStripACh1Arm1Pl0->Write();
  laStripACh1Arm1Pl1->Write();
  laStripACh1Arm1Pl2->Write();
  laStripACh1Arm1Pl3->Write();
  laStripACh1Arm1Pl4->Write();
  laStripBCh1Arm0Pl0->Write();
  laStripBCh1Arm0Pl1->Write();
  laStripBCh1Arm0Pl2->Write();
  laStripBCh1Arm0Pl3->Write();
  laStripBCh1Arm0Pl4->Write();
  laStripBCh1Arm1Pl0->Write();
  laStripBCh1Arm1Pl1->Write();
  laStripBCh1Arm1Pl2->Write();
  laStripBCh1Arm1Pl3->Write();
  laStripBCh1Arm1Pl4->Write();
  ilWireCh2Arm0Pl0->Write();
  ilWireCh2Arm0Pl1->Write();
  ilWireCh2Arm0Pl2->Write();
  ilWireCh2Arm0Pl3->Write();
  ilWireCh2Arm0Pl4->Write();
  ilWireCh2Arm1Pl0->Write();
  ilWireCh2Arm1Pl1->Write();
  ilWireCh2Arm1Pl2->Write();
  ilWireCh2Arm1Pl3->Write();
  ilWireCh2Arm1Pl4->Write();
  ilClWidthCh2Arm0Pl0->Write();
  ilClWidthCh2Arm0Pl1->Write();
  ilClWidthCh2Arm0Pl2->Write();
  ilClWidthCh2Arm0Pl3->Write();
  ilClWidthCh2Arm0Pl4->Write();
  ilClWidthCh2Arm1Pl0->Write();
  ilClWidthCh2Arm1Pl1->Write();
  ilClWidthCh2Arm1Pl2->Write();
  ilClWidthCh2Arm1Pl3->Write();
  ilClWidthCh2Arm1Pl4->Write();
  laStripACh2Arm0Pl0->Write();
  laStripACh2Arm0Pl1->Write();
  laStripACh2Arm0Pl2->Write();
  laStripACh2Arm0Pl3->Write();
  laStripACh2Arm0Pl4->Write();
  laStripACh2Arm1Pl0->Write();
  laStripACh2Arm1Pl1->Write();
  laStripACh2Arm1Pl2->Write();
  laStripACh2Arm1Pl3->Write();
  laStripACh2Arm1Pl4->Write();
  laStripBCh2Arm0Pl0->Write();
  laStripBCh2Arm0Pl1->Write();
  laStripBCh2Arm0Pl2->Write();
  laStripBCh2Arm0Pl3->Write();
  laStripBCh2Arm0Pl4->Write();
  laStripBCh2Arm1Pl0->Write();
  laStripBCh2Arm1Pl1->Write();
  laStripBCh2Arm1Pl2->Write();
  laStripBCh2Arm1Pl3->Write();
  laStripBCh2Arm1Pl4->Write();
  ilWireCh3Arm0Pl0->Write();
  ilWireCh3Arm0Pl1->Write();
  ilWireCh3Arm0Pl2->Write();
  ilWireCh3Arm0Pl3->Write();
  ilWireCh3Arm0Pl4->Write();
  ilWireCh3Arm1Pl0->Write();
  ilWireCh3Arm1Pl1->Write();
  ilWireCh3Arm1Pl2->Write();
  ilWireCh3Arm1Pl3->Write();
  ilWireCh3Arm1Pl4->Write();
  ilClWidthCh3Arm0Pl0->Write();
  ilClWidthCh3Arm0Pl1->Write();
  ilClWidthCh3Arm0Pl2->Write();
  ilClWidthCh3Arm0Pl3->Write();
  ilClWidthCh3Arm0Pl4->Write();
  ilClWidthCh3Arm1Pl0->Write();
  ilClWidthCh3Arm1Pl1->Write();
  ilClWidthCh3Arm1Pl2->Write();
  ilClWidthCh3Arm1Pl3->Write();
  ilClWidthCh3Arm1Pl4->Write();
  laStripACh3Arm0Pl0->Write();
  laStripACh3Arm0Pl1->Write();
  laStripACh3Arm0Pl2->Write();
  laStripACh3Arm0Pl3->Write();
  laStripACh3Arm0Pl4->Write();
  laStripACh3Arm1Pl0->Write();
  laStripACh3Arm1Pl1->Write();
  laStripACh3Arm1Pl2->Write();
  laStripACh3Arm1Pl3->Write();
  laStripACh3Arm1Pl4->Write();
  laStripBCh3Arm0Pl0->Write();
  laStripBCh3Arm0Pl1->Write();
  laStripBCh3Arm0Pl2->Write();
  laStripBCh3Arm0Pl3->Write();
  laStripBCh3Arm0Pl4->Write();
  laStripBCh3Arm1Pl0->Write();
  laStripBCh3Arm1Pl1->Write();
  laStripBCh3Arm1Pl2->Write();
  laStripBCh3Arm1Pl3->Write();
  laStripBCh3Arm1Pl4->Write();
  ilWireCh4Arm0Pl0->Write();
  ilWireCh4Arm0Pl1->Write();
  ilWireCh4Arm0Pl2->Write();
  ilWireCh4Arm0Pl3->Write();
  ilWireCh4Arm0Pl4->Write();
  ilWireCh4Arm1Pl0->Write();
  ilWireCh4Arm1Pl1->Write();
  ilWireCh4Arm1Pl2->Write();
  ilWireCh4Arm1Pl3->Write();
  ilWireCh4Arm1Pl4->Write();
  ilClWidthCh4Arm0Pl0->Write();
  ilClWidthCh4Arm0Pl1->Write();
  ilClWidthCh4Arm0Pl2->Write();
  ilClWidthCh4Arm0Pl3->Write();
  ilClWidthCh4Arm0Pl4->Write();
  ilClWidthCh4Arm1Pl0->Write();
  ilClWidthCh4Arm1Pl1->Write();
  ilClWidthCh4Arm1Pl2->Write();
  ilClWidthCh4Arm1Pl3->Write();
  ilClWidthCh4Arm1Pl4->Write();
  laStripACh4Arm0Pl0->Write();
  laStripACh4Arm0Pl1->Write();
  laStripACh4Arm0Pl2->Write();
  laStripACh4Arm0Pl3->Write();
  laStripACh4Arm0Pl4->Write();
  laStripACh4Arm1Pl0->Write();
  laStripACh4Arm1Pl1->Write();
  laStripACh4Arm1Pl2->Write();
  laStripACh4Arm1Pl3->Write();
  laStripACh4Arm1Pl4->Write();
  laStripBCh4Arm0Pl0->Write();
  laStripBCh4Arm0Pl1->Write();
  laStripBCh4Arm0Pl2->Write();
  laStripBCh4Arm0Pl3->Write();
  laStripBCh4Arm0Pl4->Write();
  laStripBCh4Arm1Pl0->Write();
  laStripBCh4Arm1Pl1->Write();
  laStripBCh4Arm1Pl2->Write();
  laStripBCh4Arm1Pl3->Write();
  laStripBCh4Arm1Pl4->Write();
  ilWireCh5Arm0Pl0->Write();
  ilWireCh5Arm0Pl1->Write();
  ilWireCh5Arm0Pl2->Write();
  ilWireCh5Arm0Pl3->Write();
  ilWireCh5Arm0Pl4->Write();
  ilWireCh5Arm1Pl0->Write();
  ilWireCh5Arm1Pl1->Write();
  ilWireCh5Arm1Pl2->Write();
  ilWireCh5Arm1Pl3->Write();
  ilWireCh5Arm1Pl4->Write();
  ilClWidthCh5Arm0Pl0->Write();
  ilClWidthCh5Arm0Pl1->Write();
  ilClWidthCh5Arm0Pl2->Write();
  ilClWidthCh5Arm0Pl3->Write();
  ilClWidthCh5Arm0Pl4->Write();
  ilClWidthCh5Arm1Pl0->Write();
  ilClWidthCh5Arm1Pl1->Write();
  ilClWidthCh5Arm1Pl2->Write();
  ilClWidthCh5Arm1Pl3->Write();
  ilClWidthCh5Arm1Pl4->Write();
  laStripACh5Arm0Pl0->Write();
  laStripACh5Arm0Pl1->Write();
  laStripACh5Arm0Pl2->Write();
  laStripACh5Arm0Pl3->Write();
  laStripACh5Arm0Pl4->Write();
  laStripACh5Arm1Pl0->Write();
  laStripACh5Arm1Pl1->Write();
  laStripACh5Arm1Pl2->Write();
  laStripACh5Arm1Pl3->Write();
  laStripACh5Arm1Pl4->Write();
  laStripBCh5Arm0Pl0->Write();
  laStripBCh5Arm0Pl1->Write();
  laStripBCh5Arm0Pl2->Write();
  laStripBCh5Arm0Pl3->Write();
  laStripBCh5Arm0Pl4->Write();
  laStripBCh5Arm1Pl0->Write();
  laStripBCh5Arm1Pl1->Write();
  laStripBCh5Arm1Pl2->Write();
  laStripBCh5Arm1Pl3->Write();
  laStripBCh5Arm1Pl4->Write();
  RecoX_Arm1Pl4->Write();
  RecoX_Arm1Pl3->Write();
  RecoX_Arm1Pl2->Write();
  RecoX_Arm1Pl1->Write();
  RecoX_Arm1Pl0->Write();
  RecoX_Arm0Pl4->Write();
  RecoX_Arm0Pl3->Write();
  RecoX_Arm0Pl2->Write();
  RecoX_Arm0Pl1->Write();
  RecoX_Arm0Pl0->Write();
  RecoY_Arm1Pl4->Write();
  RecoY_Arm1Pl3->Write();
  RecoY_Arm1Pl2->Write();
  RecoY_Arm1Pl1->Write();
  RecoY_Arm1Pl0->Write();
  RecoY_Arm0Pl4->Write();
  RecoY_Arm0Pl3->Write();
  RecoY_Arm0Pl2->Write();
  RecoY_Arm0Pl1->Write();
  RecoY_Arm0Pl0->Write();


  DeltaX_Arm1Pl4->Write();
  DeltaX_Arm1Pl3->Write();
  DeltaX_Arm1Pl2->Write();
  DeltaX_Arm1Pl1->Write();
  DeltaX_Arm1Pl0->Write();
  DeltaX_Arm0Pl4->Write();
  DeltaX_Arm0Pl3->Write();
  DeltaX_Arm0Pl2->Write();
  DeltaX_Arm0Pl1->Write();
  DeltaX_Arm0Pl0->Write();
  DeltaY_Arm1Pl4->Write();
  DeltaY_Arm1Pl3->Write();
  DeltaY_Arm1Pl2->Write();
  DeltaY_Arm1Pl1->Write();
  DeltaY_Arm1Pl0->Write();
  DeltaY_Arm0Pl4->Write();
  DeltaY_Arm0Pl3->Write();
  DeltaY_Arm0Pl2->Write();
  DeltaY_Arm0Pl1->Write();
  DeltaY_Arm0Pl0->Write();



  hAllEtaRec->Write();
  hAllPhiRec->Write();
  hRoadSize->Write();
  hRoadHitSigmaX->Write();
  hRoadHitSigmaY->Write();
  hRoadNumber->Write();
  hChiSquaredOverN->Write();
  hChiSquaredProb->Write();
  hChiSquaredXOverN->Write();
  hChiSquaredXProb->Write();
  hChiSquaredYOverN->Write();
  hChiSquaredYProb->Write();
  ilMomento->Write();
  ilMomentoLog->Write();
  lEnergia->Write();
  ilTipo->Write();
  ilDetector->Write();

  NWireCh0Arm0Pl0->Write();
  NWireCh0Arm0Pl1->Write();
  NWireCh0Arm0Pl2->Write();
  NWireCh0Arm0Pl3->Write();
  NWireCh0Arm0Pl4->Write();
  NWireCh0Arm1Pl0->Write();
  NWireCh0Arm1Pl1->Write();
  NWireCh0Arm1Pl2->Write();
  NWireCh0Arm1Pl3->Write();
  NWireCh0Arm1Pl4->Write();
  NStripACh0Arm0Pl0->Write();
  NStripACh0Arm0Pl1->Write();
  NStripACh0Arm0Pl2->Write();
  NStripACh0Arm0Pl3->Write();
  NStripACh0Arm0Pl4->Write();
  NStripACh0Arm1Pl0->Write();
  NStripACh0Arm1Pl1->Write();
  NStripACh0Arm1Pl2->Write();
  NStripACh0Arm1Pl3->Write();
  NStripACh0Arm1Pl4->Write();
  NStripBCh0Arm0Pl0->Write();
  NStripBCh0Arm0Pl1->Write();
  NStripBCh0Arm0Pl2->Write();
  NStripBCh0Arm0Pl3->Write();
  NStripBCh0Arm0Pl4->Write();
  NStripBCh0Arm1Pl0->Write();
  NStripBCh0Arm1Pl1->Write();
  NStripBCh0Arm1Pl2->Write();
  NStripBCh0Arm1Pl3->Write();
  NStripBCh0Arm1Pl4->Write();
  NWireCh1Arm0Pl0->Write();
  NWireCh1Arm0Pl1->Write();
  NWireCh1Arm0Pl2->Write();
  NWireCh1Arm0Pl3->Write();
  NWireCh1Arm0Pl4->Write();
  NWireCh1Arm1Pl0->Write();
  NWireCh1Arm1Pl1->Write();
  NWireCh1Arm1Pl2->Write();
  NWireCh1Arm1Pl3->Write();
  NWireCh1Arm1Pl4->Write();
  NStripACh1Arm0Pl0->Write();
  NStripACh1Arm0Pl1->Write();
  NStripACh1Arm0Pl2->Write();
  NStripACh1Arm0Pl3->Write();
  NStripACh1Arm0Pl4->Write();
  NStripACh1Arm1Pl0->Write();
  NStripACh1Arm1Pl1->Write();
  NStripACh1Arm1Pl2->Write();
  NStripACh1Arm1Pl3->Write();
  NStripACh1Arm1Pl4->Write();
  NStripBCh1Arm0Pl0->Write();
  NStripBCh1Arm0Pl1->Write();
  NStripBCh1Arm0Pl2->Write();
  NStripBCh1Arm0Pl3->Write();
  NStripBCh1Arm0Pl4->Write();
  NStripBCh1Arm1Pl0->Write();
  NStripBCh1Arm1Pl1->Write();
  NStripBCh1Arm1Pl2->Write();
  NStripBCh1Arm1Pl3->Write();
  NStripBCh1Arm1Pl4->Write();
  NWireCh2Arm0Pl0->Write();
  NWireCh2Arm0Pl1->Write();
  NWireCh2Arm0Pl2->Write();
  NWireCh2Arm0Pl3->Write();
  NWireCh2Arm0Pl4->Write();
  NWireCh2Arm1Pl0->Write();
  NWireCh2Arm1Pl1->Write();
  NWireCh2Arm1Pl2->Write();
  NWireCh2Arm1Pl3->Write();
  NWireCh2Arm1Pl4->Write();
  NStripACh2Arm0Pl0->Write();
  NStripACh2Arm0Pl1->Write();
  NStripACh2Arm0Pl2->Write();
  NStripACh2Arm0Pl3->Write();
  NStripACh2Arm0Pl4->Write();
  NStripACh2Arm1Pl0->Write();
  NStripACh2Arm1Pl1->Write();
  NStripACh2Arm1Pl2->Write();
  NStripACh2Arm1Pl3->Write();
  NStripACh2Arm1Pl4->Write();
  NStripBCh2Arm0Pl0->Write();
  NStripBCh2Arm0Pl1->Write();
  NStripBCh2Arm0Pl2->Write();
  NStripBCh2Arm0Pl3->Write();
  NStripBCh2Arm0Pl4->Write();
  NStripBCh2Arm1Pl0->Write();
  NStripBCh2Arm1Pl1->Write();
  NStripBCh2Arm1Pl2->Write();
  NStripBCh2Arm1Pl3->Write();
  NStripBCh2Arm1Pl4->Write();
  NWireCh3Arm0Pl0->Write();
  NWireCh3Arm0Pl1->Write();
  NWireCh3Arm0Pl2->Write();
  NWireCh3Arm0Pl3->Write();
  NWireCh3Arm0Pl4->Write();
  NWireCh3Arm1Pl0->Write();
  NWireCh3Arm1Pl1->Write();
  NWireCh3Arm1Pl2->Write();
  NWireCh3Arm1Pl3->Write();
  NWireCh3Arm1Pl4->Write();
  NStripACh3Arm0Pl0->Write();
  NStripACh3Arm0Pl1->Write();
  NStripACh3Arm0Pl2->Write();
  NStripACh3Arm0Pl3->Write();
  NStripACh3Arm0Pl4->Write();
  NStripACh3Arm1Pl0->Write();
  NStripACh3Arm1Pl1->Write();
  NStripACh3Arm1Pl2->Write();
  NStripACh3Arm1Pl3->Write();
  NStripACh3Arm1Pl4->Write();
  NStripBCh3Arm0Pl0->Write();
  NStripBCh3Arm0Pl1->Write();
  NStripBCh3Arm0Pl2->Write();
  NStripBCh3Arm0Pl3->Write();
  NStripBCh3Arm0Pl4->Write();
  NStripBCh3Arm1Pl0->Write();
  NStripBCh3Arm1Pl1->Write();
  NStripBCh3Arm1Pl2->Write();
  NStripBCh3Arm1Pl3->Write();
  NStripBCh3Arm1Pl4->Write();
  NWireCh4Arm0Pl0->Write();
  NWireCh4Arm0Pl1->Write();
  NWireCh4Arm0Pl2->Write();
  NWireCh4Arm0Pl3->Write();
  NWireCh4Arm0Pl4->Write();
  NWireCh4Arm1Pl0->Write();
  NWireCh4Arm1Pl1->Write();
  NWireCh4Arm1Pl2->Write();
  NWireCh4Arm1Pl3->Write();
  NWireCh4Arm1Pl4->Write();
  NStripACh4Arm0Pl0->Write();
  NStripACh4Arm0Pl1->Write();
  NStripACh4Arm0Pl2->Write();
  NStripACh4Arm0Pl3->Write();
  NStripACh4Arm0Pl4->Write();
  NStripACh4Arm1Pl0->Write();
  NStripACh4Arm1Pl1->Write();
  NStripACh4Arm1Pl2->Write();
  NStripACh4Arm1Pl3->Write();
  NStripACh4Arm1Pl4->Write();
  NStripBCh4Arm0Pl0->Write();
  NStripBCh4Arm0Pl1->Write();
  NStripBCh4Arm0Pl2->Write();
  NStripBCh4Arm0Pl3->Write();
  NStripBCh4Arm0Pl4->Write();
  NStripBCh4Arm1Pl0->Write();
  NStripBCh4Arm1Pl1->Write();
  NStripBCh4Arm1Pl2->Write();
  NStripBCh4Arm1Pl3->Write();
  NStripBCh4Arm1Pl4->Write();
  NWireCh5Arm0Pl0->Write();
  NWireCh5Arm0Pl1->Write();
  NWireCh5Arm0Pl2->Write();
  NWireCh5Arm0Pl3->Write();
  NWireCh5Arm0Pl4->Write();
  NWireCh5Arm1Pl0->Write();
  NWireCh5Arm1Pl1->Write();
  NWireCh5Arm1Pl2->Write();
  NWireCh5Arm1Pl3->Write();
  NWireCh5Arm1Pl4->Write();
  NStripACh5Arm0Pl0->Write();
  NStripACh5Arm0Pl1->Write();
  NStripACh5Arm0Pl2->Write();
  NStripACh5Arm0Pl3->Write();
  NStripACh5Arm0Pl4->Write();
  NStripACh5Arm1Pl0->Write();
  NStripACh5Arm1Pl1->Write();
  NStripACh5Arm1Pl2->Write();
  NStripACh5Arm1Pl3->Write();
  NStripACh5Arm1Pl4->Write();
  NStripBCh5Arm0Pl0->Write();
  NStripBCh5Arm0Pl1->Write();
  NStripBCh5Arm0Pl2->Write();
  NStripBCh5Arm0Pl3->Write();
  NStripBCh5Arm0Pl4->Write();
  NStripBCh5Arm1Pl0->Write();
  NStripBCh5Arm1Pl1->Write();
  NStripBCh5Arm1Pl2->Write();
  NStripBCh5Arm1Pl3->Write();
  NStripBCh5Arm1Pl4->Write();

RecoLX_Arm0Plane0CSC0->Write();
RecoLY_Arm0Plane0CSC0->Write();
NReco_Arm0Plane0CSC0->Write();
RecoLX_Arm0Plane0CSC1->Write();
RecoLY_Arm0Plane0CSC1->Write();
NReco_Arm0Plane0CSC1->Write();
RecoLX_Arm0Plane0CSC2->Write();
RecoLY_Arm0Plane0CSC2->Write();
NReco_Arm0Plane0CSC2->Write();
RecoLX_Arm0Plane0CSC3->Write();
RecoLY_Arm0Plane0CSC3->Write();
NReco_Arm0Plane0CSC3->Write();
RecoLX_Arm0Plane0CSC4->Write();
RecoLY_Arm0Plane0CSC4->Write();
NReco_Arm0Plane0CSC4->Write();
RecoLX_Arm0Plane0CSC5->Write();
RecoLY_Arm0Plane0CSC5->Write();
NReco_Arm0Plane0CSC5->Write();
RecoLX_Arm0Plane1CSC0->Write();
RecoLY_Arm0Plane1CSC0->Write();
NReco_Arm0Plane1CSC0->Write();
RecoLX_Arm0Plane1CSC1->Write();
RecoLY_Arm0Plane1CSC1->Write();
NReco_Arm0Plane1CSC1->Write();
RecoLX_Arm0Plane1CSC2->Write();
RecoLY_Arm0Plane1CSC2->Write();
NReco_Arm0Plane1CSC2->Write();
RecoLX_Arm0Plane1CSC3->Write();
RecoLY_Arm0Plane1CSC3->Write();
NReco_Arm0Plane1CSC3->Write();
RecoLX_Arm0Plane1CSC4->Write();
RecoLY_Arm0Plane1CSC4->Write();
NReco_Arm0Plane1CSC4->Write();
RecoLX_Arm0Plane1CSC5->Write();
RecoLY_Arm0Plane1CSC5->Write();
NReco_Arm0Plane1CSC5->Write();
RecoLX_Arm0Plane2CSC0->Write();
RecoLY_Arm0Plane2CSC0->Write();
NReco_Arm0Plane2CSC0->Write();
RecoLX_Arm0Plane2CSC1->Write();
RecoLY_Arm0Plane2CSC1->Write();
NReco_Arm0Plane2CSC1->Write();
RecoLX_Arm0Plane2CSC2->Write();
RecoLY_Arm0Plane2CSC2->Write();
NReco_Arm0Plane2CSC2->Write();
RecoLX_Arm0Plane2CSC3->Write();
RecoLY_Arm0Plane2CSC3->Write();
NReco_Arm0Plane2CSC3->Write();
RecoLX_Arm0Plane2CSC4->Write();
RecoLY_Arm0Plane2CSC4->Write();
NReco_Arm0Plane2CSC4->Write();
RecoLX_Arm0Plane2CSC5->Write();
RecoLY_Arm0Plane2CSC5->Write();
NReco_Arm0Plane2CSC5->Write();
RecoLX_Arm0Plane3CSC0->Write();
RecoLY_Arm0Plane3CSC0->Write();
NReco_Arm0Plane3CSC0->Write();
RecoLX_Arm0Plane3CSC1->Write();
RecoLY_Arm0Plane3CSC1->Write();
NReco_Arm0Plane3CSC1->Write();
RecoLX_Arm0Plane3CSC2->Write();
RecoLY_Arm0Plane3CSC2->Write();
NReco_Arm0Plane3CSC2->Write();
RecoLX_Arm0Plane3CSC3->Write();
RecoLY_Arm0Plane3CSC3->Write();
NReco_Arm0Plane3CSC3->Write();
RecoLX_Arm0Plane3CSC4->Write();
RecoLY_Arm0Plane3CSC4->Write();
NReco_Arm0Plane3CSC4->Write();
RecoLX_Arm0Plane3CSC5->Write();
RecoLY_Arm0Plane3CSC5->Write();
NReco_Arm0Plane3CSC5->Write();
RecoLX_Arm0Plane4CSC0->Write();
RecoLY_Arm0Plane4CSC0->Write();
NReco_Arm0Plane4CSC0->Write();
RecoLX_Arm0Plane4CSC1->Write();
RecoLY_Arm0Plane4CSC1->Write();
NReco_Arm0Plane4CSC1->Write();
RecoLX_Arm0Plane4CSC2->Write();
RecoLY_Arm0Plane4CSC2->Write();
NReco_Arm0Plane4CSC2->Write();
RecoLX_Arm0Plane4CSC3->Write();
RecoLY_Arm0Plane4CSC3->Write();
NReco_Arm0Plane4CSC3->Write();
RecoLX_Arm0Plane4CSC4->Write();
RecoLY_Arm0Plane4CSC4->Write();
NReco_Arm0Plane4CSC4->Write();
RecoLX_Arm0Plane4CSC5->Write();
RecoLY_Arm0Plane4CSC5->Write();
NReco_Arm0Plane4CSC5->Write();
RecoLX_Arm1Plane0CSC0->Write();
RecoLY_Arm1Plane0CSC0->Write();
NReco_Arm1Plane0CSC0->Write();
RecoLX_Arm1Plane0CSC1->Write();
RecoLY_Arm1Plane0CSC1->Write();
NReco_Arm1Plane0CSC1->Write();
RecoLX_Arm1Plane0CSC2->Write();
RecoLY_Arm1Plane0CSC2->Write();
NReco_Arm1Plane0CSC2->Write();
RecoLX_Arm1Plane0CSC3->Write();
RecoLY_Arm1Plane0CSC3->Write();
NReco_Arm1Plane0CSC3->Write();
RecoLX_Arm1Plane0CSC4->Write();
RecoLY_Arm1Plane0CSC4->Write();
NReco_Arm1Plane0CSC4->Write();
RecoLX_Arm1Plane0CSC5->Write();
RecoLY_Arm1Plane0CSC5->Write();
NReco_Arm1Plane0CSC5->Write();
RecoLX_Arm1Plane1CSC0->Write();
RecoLY_Arm1Plane1CSC0->Write();
NReco_Arm1Plane1CSC0->Write();
RecoLX_Arm1Plane1CSC1->Write();
RecoLY_Arm1Plane1CSC1->Write();
NReco_Arm1Plane1CSC1->Write();
RecoLX_Arm1Plane1CSC2->Write();
RecoLY_Arm1Plane1CSC2->Write();
NReco_Arm1Plane1CSC2->Write();
RecoLX_Arm1Plane1CSC3->Write();
RecoLY_Arm1Plane1CSC3->Write();
NReco_Arm1Plane1CSC3->Write();
RecoLX_Arm1Plane1CSC4->Write();
RecoLY_Arm1Plane1CSC4->Write();
NReco_Arm1Plane1CSC4->Write();
RecoLX_Arm1Plane1CSC5->Write();
RecoLY_Arm1Plane1CSC5->Write();
NReco_Arm1Plane1CSC5->Write();
RecoLX_Arm1Plane2CSC0->Write();
RecoLY_Arm1Plane2CSC0->Write();
NReco_Arm1Plane2CSC0->Write();
RecoLX_Arm1Plane2CSC1->Write();
RecoLY_Arm1Plane2CSC1->Write();
NReco_Arm1Plane2CSC1->Write();
RecoLX_Arm1Plane2CSC2->Write();
RecoLY_Arm1Plane2CSC2->Write();
NReco_Arm1Plane2CSC2->Write();
RecoLX_Arm1Plane2CSC3->Write();
RecoLY_Arm1Plane2CSC3->Write();
NReco_Arm1Plane2CSC3->Write();
RecoLX_Arm1Plane2CSC4->Write();
RecoLY_Arm1Plane2CSC4->Write();
NReco_Arm1Plane2CSC4->Write();
RecoLX_Arm1Plane2CSC5->Write();
RecoLY_Arm1Plane2CSC5->Write();
NReco_Arm1Plane2CSC5->Write();
RecoLX_Arm1Plane3CSC0->Write();
RecoLY_Arm1Plane3CSC0->Write();
NReco_Arm1Plane3CSC0->Write();
RecoLX_Arm1Plane3CSC1->Write();
RecoLY_Arm1Plane3CSC1->Write();
NReco_Arm1Plane3CSC1->Write();
RecoLX_Arm1Plane3CSC2->Write();
RecoLY_Arm1Plane3CSC2->Write();
NReco_Arm1Plane3CSC2->Write();
RecoLX_Arm1Plane3CSC3->Write();
RecoLY_Arm1Plane3CSC3->Write();
NReco_Arm1Plane3CSC3->Write();
RecoLX_Arm1Plane3CSC4->Write();
RecoLY_Arm1Plane3CSC4->Write();
NReco_Arm1Plane3CSC4->Write();
RecoLX_Arm1Plane3CSC5->Write();
RecoLY_Arm1Plane3CSC5->Write();
NReco_Arm1Plane3CSC5->Write();
RecoLX_Arm1Plane4CSC0->Write();
RecoLY_Arm1Plane4CSC0->Write();
NReco_Arm1Plane4CSC0->Write();
RecoLX_Arm1Plane4CSC1->Write();
RecoLY_Arm1Plane4CSC1->Write();
NReco_Arm1Plane4CSC1->Write();
RecoLX_Arm1Plane4CSC2->Write();
RecoLY_Arm1Plane4CSC2->Write();
NReco_Arm1Plane4CSC2->Write();
RecoLX_Arm1Plane4CSC3->Write();
RecoLY_Arm1Plane4CSC3->Write();
NReco_Arm1Plane4CSC3->Write();
RecoLX_Arm1Plane4CSC4->Write();
RecoLY_Arm1Plane4CSC4->Write();
NReco_Arm1Plane4CSC4->Write();
RecoLX_Arm1Plane4CSC5->Write();
RecoLY_Arm1Plane4CSC5->Write();
NReco_Arm1Plane4CSC5->Write();




  RecoXY_00 -> Write();
  RecoXY_01 -> Write();
  RecoXY_02 -> Write();
  RecoXY_03 -> Write();
  RecoXY_04 -> Write();
  RecoXY_10 -> Write();
  RecoXY_11 -> Write();
  RecoXY_12 -> Write();
  RecoXY_13 -> Write();
  RecoXY_14 -> Write();

}

DEFINE_FWK_MODULE(T1Validation);
