// Description: Function to retrieve Gsf Muon 



// Genaral Tracks and Vertex
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// Muon Includes
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"



#include "UATree/UABaseTree/interface/UABaseTree.h"


Bool_t MuonDebug = false;

void UABaseTree::GetRecoMuon(const edm::Event& iEvent, const InputTag& globalmuoncoll_ , vector<MyMuon>& MuonVector ){

   MuonVector.clear();

   Handle<MuonCollection> muonsHandle;
   try {
     iEvent.getByLabel(globalmuoncoll_,muonsHandle);
   } catch ( cms::Exception& ex ) {
    printf("Error! can't get globalMuon collection\n");
   }
   const MuonCollection & muons = *(muonsHandle.product());
   
   MuonVector.assign( muons.size() , MyMuon() );

   Int_t i = 0;
   for (MuonCollection::const_iterator iMuon = muons.begin() ; iMuon != muons.end() ; ++iMuon , ++i) {
     
     MuonVector[i].SetPxPyPzE(iMuon->px() , iMuon->py() , iMuon->pz() , sqrt(iMuon->momentum().mag2()+MASS_MU*MASS_MU));
     MuonVector[i].charge = iMuon->charge();

     // Global Properties
     MuonVector[i].nChambers             = iMuon -> numberOfChambers();
     MuonVector[i].nChambersMatched      = iMuon -> numberOfMatches();
     
     MuonVector[i].AllGlobalMuons                          =  muon::isGoodMuon(*iMuon,muon::AllGlobalMuons);
     MuonVector[i].AllStandAloneMuons                      =  muon::isGoodMuon(*iMuon,muon::AllStandAloneMuons);
     MuonVector[i].AllTrackerMuons                         =  muon::isGoodMuon(*iMuon,muon::AllTrackerMuons);
     MuonVector[i].TrackerMuonArbitrated                   =  muon::isGoodMuon(*iMuon,muon::TrackerMuonArbitrated);
     MuonVector[i].AllArbitrated                           =  muon::isGoodMuon(*iMuon,muon::AllArbitrated);
     MuonVector[i].GlobalMuonPromptTight                   =  muon::isGoodMuon(*iMuon,muon::GlobalMuonPromptTight);
     MuonVector[i].TMLastStationLoose                      =  muon::isGoodMuon(*iMuon,muon::TMLastStationLoose);
     MuonVector[i].TMLastStationTight                      =  muon::isGoodMuon(*iMuon,muon::TMLastStationTight);
     MuonVector[i].TM2DCompatibilityLoose                  =  muon::isGoodMuon(*iMuon,muon::TM2DCompatibilityLoose);
     MuonVector[i].TM2DCompatibilityTight                  =  muon::isGoodMuon(*iMuon,muon::TM2DCompatibilityTight);
     MuonVector[i].TMOneStationLoose                       =  muon::isGoodMuon(*iMuon,muon::TMOneStationLoose);
     MuonVector[i].TMOneStationTight                       =  muon::isGoodMuon(*iMuon,muon::TMOneStationTight);
     MuonVector[i].TMLastStationOptimizedLowPtLoose        =  muon::isGoodMuon(*iMuon,muon::TMLastStationOptimizedLowPtLoose);
     MuonVector[i].TMLastStationOptimizedLowPtTight        =  muon::isGoodMuon(*iMuon,muon::TMLastStationOptimizedLowPtTight);
     MuonVector[i].GMTkChiCompatibility                    =  muon::isGoodMuon(*iMuon,muon::GMTkChiCompatibility);
     MuonVector[i].GMStaChiCompatibility                   =  muon::isGoodMuon(*iMuon,muon::GMStaChiCompatibility);
     MuonVector[i].GMTkKinkTight                           =  muon::isGoodMuon(*iMuon,muon::GMTkKinkTight);
     MuonVector[i].TMLastStationAngLoose                   =  muon::isGoodMuon(*iMuon,muon::TMLastStationAngLoose);
     MuonVector[i].TMLastStationAngTight                   =  muon::isGoodMuon(*iMuon,muon::TMLastStationAngTight);
     MuonVector[i].TMOneStationAngLoose                    =  muon::isGoodMuon(*iMuon,muon::TMOneStationAngLoose);
     MuonVector[i].TMOneStationAngTight                    =  muon::isGoodMuon(*iMuon,muon::TMOneStationAngTight);
     MuonVector[i].TMLastStationOptimizedBarrelLowPtLoose  =  muon::isGoodMuon(*iMuon,muon::TMLastStationOptimizedBarrelLowPtLoose);
     MuonVector[i].TMLastStationOptimizedBarrelLowPtTight  =  muon::isGoodMuon(*iMuon,muon::TMLastStationOptimizedBarrelLowPtTight);
     
 
     MuonVector[i].isoR03sumPt   = iMuon->isolationR03().sumPt   ;
     MuonVector[i].isoR03emEt    = iMuon->isolationR03().emEt    ;
     MuonVector[i].isoR03hadEt   = iMuon->isolationR03().hadEt   ;
     MuonVector[i].isoR03hoEt    = iMuon->isolationR03().hoEt    ;
     MuonVector[i].isoR03nTracks = iMuon->isolationR03().nTracks ;
     MuonVector[i].isoR03nJets   = iMuon->isolationR03().nJets   ;

     MuonVector[i].isoR05sumPt   = iMuon->isolationR05().sumPt   ;
     MuonVector[i].isoR05emEt    = iMuon->isolationR05().emEt    ;
     MuonVector[i].isoR05hadEt   = iMuon->isolationR05().hadEt   ;
     MuonVector[i].isoR05hoEt    = iMuon->isolationR05().hoEt    ;
     MuonVector[i].isoR05nTracks = iMuon->isolationR05().nTracks ;
     MuonVector[i].isoR05nJets   = iMuon->isolationR05().nJets   ;

     MuonVector[i].calEnergyEm   = iMuon->calEnergy().em   ; 
     MuonVector[i].calEnergyHad  = iMuon->calEnergy().had  ;
     MuonVector[i].calEnergyHo   = iMuon->calEnergy().ho   ;
     MuonVector[i].calEnergyEmS9 = iMuon->calEnergy().emS9 ;
     MuonVector[i].calEnergyHadS9= iMuon->calEnergy().hadS9;
     MuonVector[i].calEnergyHoS9 = iMuon->calEnergy().hoS9 ;

     MuonVector[i].IsGlobalMuon      = iMuon->isGlobalMuon()     ;
     MuonVector[i].IsTrackerMuon     = iMuon->isTrackerMuon()    ;
     MuonVector[i].IsStandaloneMuon  = iMuon->isStandAloneMuon() ;
     MuonVector[i].IsCaloMuon        = iMuon->isCaloMuon()       ;

     
     
 

     // Global Muon Track
     if(iMuon->globalTrack().isAvailable())
       this->FillTrack( *(iMuon->globalTrack()) ,  MuonVector[i].globalTrack , MASS_MU);
       
       
     // Inner Muon Track
     if(iMuon->innerTrack().isAvailable())
       this->FillTrack( *(iMuon->innerTrack()) ,  MuonVector[i].innerTrack , MASS_MU);
       
       
     // Outer Muon Track
     if(iMuon->outerTrack().isAvailable())
       this->FillTrack( *(iMuon->outerTrack()) ,  MuonVector[i].outerTrack , MASS_MU);
     


     if(MuonDebug) MuonVector[i].Print();


   } // end for MuonCollection 

}


void UABaseTree::GetAllMuons( const edm::Event& iEvent ){
  for(vector<InputTag>::iterator icoll = muons_.begin() ; icoll!= muons_.end() ; ++icoll)
    this->GetRecoMuon(iEvent , *icoll , allMuons[icoll->label()] );
}
