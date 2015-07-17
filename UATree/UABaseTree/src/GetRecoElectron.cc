// Description: Function to retrieve Gsf Electron 


// Gsf Electron
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronIsoCollection.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronIsoNumCollection.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Common/interface/ValueMap.h"

// Candidates
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

// Isolation
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

//Tracks and conversion
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"



#include "UATree/UABaseTree/interface/UABaseTree.h"


Bool_t ElectronDebug = false;

void UABaseTree::GetRecoElectron(const edm::Event& iEvent, const InputTag& gsfelectroncoll_ ,vector<MyElectron>& ElecVector ){

  ElecVector.clear();

  // Get Gsf Electron Collection 
  Handle<reco::GsfElectronCollection> gsfElectronsHandle;
  try {
    iEvent.getByLabel(gsfelectroncoll_,gsfElectronsHandle);
  } catch ( cms::Exception& ex ) {
    printf("Error! can't get GsfElectronCollection_ collection \n");
  }
  const GsfElectronCollection gsfElectrons = *(gsfElectronsHandle.product());


  // Another way of getting electron info
  edm::Handle< edm::View<reco::Candidate> > emObjectHandle;
  try {
    iEvent.getByLabel("pixelMatchGsfElectrons",emObjectHandle);
  }
  catch ( cms::Exception& ex ) {
    printf("Can't get emObjectHandle Collection\n");
  }
 
  ElecVector.assign(gsfElectrons.size() , MyElectron());

  
  // Loop on Gsf Electron
  Int_t i = 0;
  for (GsfElectronCollection::const_iterator iElectron = gsfElectrons.begin() ; iElectron != gsfElectrons.end() ; ++iElectron , ++i){
   
  
    unsigned int iEl = iElectron - gsfElectrons.begin();
    GsfElectronRef eRef(gsfElectronsHandle,iEl);

    ElecVector[i].SetPxPyPzE(iElectron->px() , iElectron->py() , iElectron->pz(), iElectron->superCluster()->energy());
    ElecVector[i].charge = iElectron->charge(); 

    ElecVector[i].eSupClusOverP      = iElectron->eSuperClusterOverP()  	    ;
    ElecVector[i].eSeedClusOverPout  = iElectron->eSeedClusterOverPout()	    ;
    ElecVector[i].PIn		     = iElectron->trackMomentumAtVtx().R();
    ElecVector[i].POut  	     = iElectron->trackMomentumOut().R();
    
    
    ElecVector[i].dEtaSupClusTrVtx   = iElectron->deltaEtaSuperClusterTrackAtVtx()  ;
    ElecVector[i].dEtaSeedClusTrCalo = iElectron->deltaEtaSeedClusterTrackAtCalo()  ; 
    ElecVector[i].dPhiSupClusTrVtx   = iElectron->deltaPhiSuperClusterTrackAtVtx()  ; 
    ElecVector[i].dPhiSeedClusTrCalo = iElectron->deltaPhiSeedClusterTrackAtCalo()  ; 
    ElecVector[i].isEScaleCorr       = iElectron->isEnergyScaleCorrected()	    ; 
    //    ElecVector[i].isMomentumCorr     = iElectron->isMomentumCorrected() 	    ;  deprecated !
    ElecVector[i].isEcalDriven       = iElectron->ecalDrivenSeed()		    ;
    ElecVector[i].isTrackerDriven    = iElectron->trackerDrivenSeed()		    ;
    ElecVector[i].nClus 	     = iElectron->basicClustersSize()		    ; 
    ElecVector[i].classification     = iElectron->classification()		    ; 
    ElecVector[i].fbrem 	     = iElectron->fbrem()			    ; 
    
    //fiducial
    ElecVector[i].isBarrel	     = iElectron->isEB()			    ;
    ElecVector[i].isEndCap	     = iElectron->isEE()			    ;
    
    // Shower Shape variables
    ElecVector[i].E15			  = iElectron->e1x5();
    ElecVector[i].E25Max		  = iElectron->e2x5Max();
    ElecVector[i].E55			  = iElectron->e5x5();
    ElecVector[i].CovEtaEta		  = iElectron->sigmaEtaEta();
    ElecVector[i].CoviEtaiEta		  = iElectron->sigmaIetaIeta();
    ElecVector[i].HadronicOverEm	  = iElectron->hcalOverEcal();
    ElecVector[i].HcalDepth1OverEcal	  = iElectron->hcalDepth1OverEcal();
    ElecVector[i].HcalDepth2OverEcal	  = iElectron->hcalDepth2OverEcal();

    try{
      edm::Handle<edm::ValueMap<int> > vmEl;
      iEvent.getByLabel("expectedHitsEle",vmEl);
    
      // fill corrected expected inner hits
      ElecVector[i].expectedInnerHits      =  (*vmEl)[eRef];
    }
    catch ( ... ) {
      //printf("Can't access expectedHitsEle\n");
    }
   
    
    
    // Gsf Track Info
    if (iElectron->gsfTrack().isNonnull())
      this->FillTrack( *(iElectron->gsfTrack()) , ElecVector[i].GsfTrack , MASS_EL);

    
    // Closest CTF Track Info
    if (iElectron->closestCtfTrackRef().isNonnull())
      this->FillTrack( *(iElectron->closestCtfTrackRef()) , ElecVector[i].GsfTrack , MASS_EL);


    
    //Isolation Variables
    ElecVector[i].EcalRecHitIsoDr04	      = iElectron ->dr04EcalRecHitSumEt();
    ElecVector[i].HcalDepth1TowerSumEtDr04    = iElectron ->dr04HcalDepth1TowerSumEt();
    ElecVector[i].HcalDepth2TowerSumEtDr04    = iElectron ->dr04HcalDepth2TowerSumEt();
    ElecVector[i].TrackIsolationDr04	      = iElectron ->dr04TkSumPt();
    ElecVector[i].EcalRecHitIsoDr03	      = iElectron ->dr03EcalRecHitSumEt();
    ElecVector[i].HcalTowerSumEtDr03	      = iElectron ->dr03HcalTowerSumEt();
    ElecVector[i].HcalDepth1TowerSumEtDr03    = iElectron ->dr03HcalDepth1TowerSumEt();
    ElecVector[i].HcalDepth2TowerSumEtDr03    = iElectron ->dr03HcalDepth2TowerSumEt();
    ElecVector[i].TrackIsolationDr03	      = iElectron ->dr03TkSumPt();
    	 
    
    // Conversion variables
    edm::Handle<TrackCollection> tracks;
    iEvent.getByLabel("generalTracks", tracks);
    ConversionFinder convFinder;
    
    ConversionInfo convInfo = convFinder.getConversionInfo(*iElectron, tracks, 3.8112);
    ElecVector[i].dist_conv   = convInfo.dist();
    ElecVector[i].dcot_conv   = convInfo.dcot();
    
    
    
    // Id boolean (don't use them now)

    try{
      //Read eID results
      std::vector<edm::Handle<edm::ValueMap<float> > > eIDValueMap(4); 
      //Robust-Loose 
      iEvent.getByLabel( "eidRobustLoose" , eIDValueMap[0] ); 
      const edm::ValueMap<float> & eIDmapRL = * eIDValueMap[0] ;
      //Robust-Tight 
      iEvent.getByLabel( "eidRobustTight" , eIDValueMap[1] ); 
      const edm::ValueMap<float> & eIDmapRT = * eIDValueMap[1] ;
      //Loose 
      iEvent.getByLabel( "eidLoose" , eIDValueMap[2] ); 
      const edm::ValueMap<float> & eIDmapL = * eIDValueMap[2] ;
      //Tight 
      iEvent.getByLabel( "eidTight" , eIDValueMap[3] ); 
      const edm::ValueMap<float> & eIDmapT = * eIDValueMap[3] ;
 
      ElecVector[i].eidRobustLoose = eIDmapRL[eRef] ;
      ElecVector[i].eidRobustTight = eIDmapRT[eRef] ;
      ElecVector[i].eidLoose       = eIDmapL[eRef]  ;
      ElecVector[i].eidTight       = eIDmapT[eRef]  ;
    }
    catch ( cms::Exception& ex ){
      //printf("Can't store Tk isolation variables from Majid\n");
    }
    
    
    if(ElectronDebug) ElecVector[i].Print();
    
  } // End Loop on Gsf Electron


}



void UABaseTree::GetAllElectrons( const edm::Event& iEvent ){
  for(vector<InputTag>::iterator icoll = electrons_.begin() ; icoll!= electrons_.end() ; ++icoll)
    this->GetRecoElectron(iEvent , *icoll , allElectrons[icoll->label()] );
}
