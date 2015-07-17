// Description: Function to retrieve all the classes

// UABaseTree Analysis class decleration
#include "UATree/UABaseTree/interface/UABaseTree.h"

bool AllDebug = false;

void UABaseTree::GetAll(const edm::Event& iEvent , const edm::EventSetup& iSetup){

  // General Info
                                    GetAllBeamSpots(iEvent);
  if(storeEvtId_)                   GetEvtId(iEvent);
  if(storeFwdGap_)                  GetFwdGap(iEvent); 

  if(storeGenKin_)		    GetGenKin (iEvent);
  if(storeGenPart_ || 
     enableGenMetFromGenPart_)      GetGenPart(iEvent,iSetup);
  if(storePUSumInfo_)               GetPUSumInfo(iEvent);

  if(hlt_paths_.size() > 0)         GetHLTrig(iEvent,iSetup);
  //  if(storeL1Trig_)		    GetL1Trig(iEvent,iSetup,iRun);
  if(storeL1TrigOld_)		    GetL1TrigOld(iEvent);
  
  if(storeMITEvtSel_)               GetMITEvtSel(iEvent);

  // Reset vtx id and vector
  vtxid = 0;
  vtxid_xyz.clear();

  GetAllVertices(iEvent); // Needs to be before tracks & PFJets !!
  GetAllTracks(iEvent); 

  GetAllCaloJets(iEvent,iSetup);
  GetAllPFJets(iEvent,iSetup); 
  GetAllGenJets(iEvent); 
  GetAllBasicJets(iEvent); 
  GetAllTrackJets(iEvent, iSetup); 
  
  if(castorrechits_.label().size() > 0) GetCastorRecHit(iEvent);
  if(castorjets_.label().size() > 0 &&
     castorjetid_.label().size() > 0)   GetCastorJet(iEvent); 
  if(castordigis_.label().size() > 0)   GetCastorDigi(iEvent,iSetup);
  
  GetAllElectrons(iEvent); 
  GetAllMuons(iEvent); 
  GetAllPFCands(iEvent); 
  
  GetAllMETs(iEvent); 

  if(storeCaloObjects_) GetCaloTower(iEvent);

  if( (storeZDCInfo_  && zdcrechits_.label().size() > 0) ||
      (storeZDCHits_  && zdcrechits_.label().size() > 0) ||
      (storeZDCDigis_ && zdcdigis_.label().size() > 0) ) GetZDCInfo(iEvent,iSetup);

  if( (storeFSCInfo_  && fscrechits_.label().size() > 0) ||
      (storeFSCHits_  && fscrechits_.label().size() > 0) ||
      (storeFSCDigis_ && fscdigis_.label().size() > 0) ) GetFSCInfo(iEvent,iSetup);
}
