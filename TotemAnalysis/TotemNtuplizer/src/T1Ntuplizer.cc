
#include "TotemAnalysis/TotemNtuplizer/interface/T1Ntuplizer.h"

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include "CLHEP/Vector/LorentzVector.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TTree.h"
#include "TMath.h"


T1Ntuplizer::T1Ntuplizer(const edm::ParameterSet& iConfig) : Ntuplizer(iConfig)
{

/*
  CluLabel = iConfig.getParameter<std::string>("CluLabel");
  HitLabel = iConfig.getParameter<std::string>("HitLabel");
  RoadLabel = iConfig.getParameter<std::string>("RoadLabel");
  TrackLabel= iConfig.getParameter<std::string>("TrackLabel");
  std::vector<int> qused;qused.push_back(0);qused.push_back(1);qused.push_back(2);qused.push_back(3);
  
  T1CutsUtil.SetCuts(4.5,7.5,4,11,2.,0.01,0.01,qused,0.001,true);
*/
  std::string thegeometryfile = "Geometry/TotemGeometry/data/T1_data_geometry.dat";
  layer = std::auto_ptr<T1Geometry>(new T1Geometry(thegeometryfile));
  trackLabel= iConfig.getParameter<std::string>("T1TrackLabel");
  t1DigiWireCollectionLabel = iConfig.getParameter<edm::InputTag>("T1DigiWireCollectionLabel");
  t1DigiVfatCollectionLabel = iConfig.getParameter<edm::InputTag>("T1DigiVfatCollectionLabel");
  t1RecHit2DCollectionLabel = iConfig.getParameter<edm::InputTag>("T1RecHit2DCollectionLabel");

}



//
// member functions
//
// ------------ method called to for each event  ------------
void T1Ntuplizer::FillEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  /* LOADING OF ALL THE RECORDS FROM THE EVENT */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

    edm::Handle<T1DigiWireCollection> myDigiColl;
    iEvent.getByLabel(t1DigiWireCollectionLabel, myDigiColl);

    edm::Handle<T1DigiVfatCollection> myDigiCollStrip;
    iEvent.getByLabel(t1DigiVfatCollectionLabel, myDigiCollStrip);

    edm::Handle<T1RecHit2DCollection>myRecoColl;
    iEvent.getByLabel(t1RecHit2DCollectionLabel, myRecoColl);
  
    edm::Handle<T1T2TrackCollection> trackCollection;
    iEvent.getByLabel(trackLabel,"T1TrackColl",trackCollection);
  
 
 //Fill the T2 ntple

 t1obj.ev_no = iEvent.id().event();
 t1obj.run_no =iEvent.id().run();
 t1obj.timestamp = iEvent.time().unixTime();

 t1obj.DigiStripA.clear();           
 t1obj.DigiStripB.clear(); 
 t1obj.DigiWire.clear();
 t1obj.DigiStripB_Arm.clear();           
 t1obj.DigiStripB_Plane.clear();           
 t1obj.DigiStripB_CSC.clear();           

  t1obj.DigiStripA_Arm.clear();           
 t1obj.DigiStripA_Plane.clear();           
 t1obj.DigiStripA_CSC.clear();           

 t1obj.DigiWire_Arm.clear();           
 t1obj.DigiWire_Plane.clear();           
 t1obj.DigiWire_CSC.clear();           

  
 t1obj.TrkAx.clear();           
 t1obj.TrkAy.clear();          
 t1obj.TrkX0.clear();          
 t1obj.TrkY0.clear();         
 t1obj.TrkPhi.clear();          
 t1obj.TrkEta.clear();  
 t1obj.TrkZatRmin.clear();
 t1obj.TrkRmin.clear();

 t1obj.TrkChi2OverN.clear();
 t1obj.TrkHits.clear();

 
 // track hits?????


 
 t1obj.RecoHitX.clear();
 t1obj.RecoHitY.clear();
 t1obj.RecoHitZ.clear();
 

 
//loop over digis

   T1DigiWireCollection::DigiRangeIterator T1Digi_it;
    for(T1Digi_it=myDigiColl->begin();T1Digi_it!=myDigiColl->end();T1Digi_it++){

      int braccio = (*T1Digi_it).first.Arm();
      int piano = (*T1Digi_it).first.Plane();
      int camera = (*T1Digi_it).first.CSC();

      const T1DigiWireCollection::Range& range = (*T1Digi_it).second;
      for (T1DigiWireCollection::const_iterator digiIt = range.first;
           digiIt!=range.second;++digiIt){

	t1obj.DigiWire.push_back(digiIt->wire());
	t1obj.DigiWire_Arm.push_back(braccio);
	t1obj.DigiWire_Plane.push_back(piano);
	t1obj.DigiWire_CSC.push_back(camera);

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
	
	if(strato==1){
	  t1obj.DigiStripA.push_back(digiIt->strip());
	  t1obj.DigiStripA_Arm.push_back(braccio);
	  t1obj.DigiStripA_Plane.push_back(piano);
	  t1obj.DigiStripA_CSC.push_back(camera);
	}
	else if(strato==2){
	  t1obj.DigiStripB.push_back(digiIt->strip());
	  t1obj.DigiStripB_Arm.push_back(braccio);
	  t1obj.DigiStripB_Plane.push_back(piano);
	  t1obj.DigiStripB_CSC.push_back(camera);
	}else{
	  assert(0>1);
	}


      }
    }




    T1RecHit2DCollection::const_iterator T1Reco_it;

    for(T1Reco_it=myRecoColl->begin();T1Reco_it!=myRecoColl->end();T1Reco_it++){

      T1DetId t1Id((*T1Reco_it).t1DetId());

      float xxxR = layer->xFromLocal2BeamSystem((*T1Reco_it).t1DetId().rawId(),(*T1Reco_it).localPosition().x() );
      float yyyR = layer->yFromLocal2BeamSystem((*T1Reco_it).t1DetId().rawId(),(*T1Reco_it).localPosition().y() );

      float xxxxR=0;
      float yyyyR=0;
      layer->RotationLocal2Global((*T1Reco_it).t1DetId().rawId(),xxxR,yyyR,xxxxR,yyyyR);
      float zzzzR = layer->Zeta((*T1Reco_it).t1DetId().rawId());


      t1obj.RecoHitX.push_back(xxxxR);

      t1obj.RecoHitY.push_back(yyyyR);

      t1obj.RecoHitZ.push_back(zzzzR);

    }


    T1T2TrackCollection::const_iterator TC_it;

    for(TC_it=trackCollection->begin(); TC_it!=trackCollection->end(); TC_it++){
      
      t1obj.TrkAx.push_back( (*TC_it).GetTx() );
      t1obj.TrkAy.push_back( (*TC_it).GetTy() );
      t1obj.TrkX0.push_back( (*TC_it).X0() );
      t1obj.TrkY0.push_back( (*TC_it).Y0() );
      t1obj.TrkPhi.push_back( (*TC_it).Phi() );          
      t1obj.TrkEta.push_back( (*TC_it).Eta() );  
      t1obj.TrkZatRmin.push_back( (*TC_it).Z_at_Rmin() );
      t1obj.TrkRmin.push_back( (*TC_it).Rmin() );
      t1obj.TrkChi2OverN.push_back( (*TC_it).ChiSquaredOverN() );
      t1obj.TrkHits.push_back( (*TC_it).GetHitEntries() );
    }




}








// ------------ method called once each job just before starting event loop  ------------
void T1Ntuplizer::CreateBranches(const edm::EventSetup&, TTree *tree)
{
  TH1::AddDirectory(kFALSE);


 
  tree->Branch("branchT1EV.",&t1obj);
}
