#ifndef __UABASETREE_H__
#define __UABASETREE_H__


// system include files
#include <string>
#include <vector>
#include <iostream>

// ROOT
#include "TFile.h"
#include "TTree.h"


// CMSSW Include files (Minimal)
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Point 3D
#include "DataFormats/Math/interface/Point3D.h"

// Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// Track
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//L1
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"


// UATree Analysis class declaration

#include "UATree/UADataFormat/src/MassParticles.h"

#include "UATree/UADataFormat/src/MyPart.h"

#include "UATree/UADataFormat/src/MyEvtId.h"
#include "UATree/UADataFormat/src/MyL1Trig.h"
#include "UATree/UADataFormat/src/MyL1TrigOld.h"
#include "UATree/UADataFormat/src/MyHLTrig.h"

#include "UATree/UADataFormat/src/MyGenKin.h"
#include "UATree/UADataFormat/src/MyGenMet.h"
#include "UATree/UADataFormat/src/MyGenPart.h"
#include "UATree/UADataFormat/src/MySimVertex.h"
#include "UATree/UADataFormat/src/MyPUSumInfo.h"

#include "UATree/UADataFormat/src/MyBeamSpot.h"
#include "UATree/UADataFormat/src/MyVertex.h"
#include "UATree/UADataFormat/src/MyTracks.h"

#include "UATree/UADataFormat/src/MyFwdGap.h"

// MIT code
#include "UATree/UADataFormat/src/MyMITEvtSel.h"

// Jets
#include "UATree/UADataFormat/src/MyCaloJet.h"
#include "UATree/UADataFormat/src/MyTrackJet.h"
#include "UATree/UADataFormat/src/MyPFJet.h"
#include "UATree/UADataFormat/src/MyGenJet.h"

//Castor
#include "UATree/UADataFormat/src/MyCastorRecHit.h"
#include "UATree/UADataFormat/src/MyCastorJet.h"
#include "UATree/UADataFormat/src/MyCastorDigi.h"
#include "UATree/UADataFormat/src/MyDiJet.h"

//Identified Particles
#include "UATree/UADataFormat/src/MyElectron.h"
#include "UATree/UADataFormat/src/MyMuon.h"
#include "UATree/UADataFormat/src/MyPFCand.h"

#include "UATree/UADataFormat/src/MyMet.h"

//CaloObjects
#include "UATree/UADataFormat/src/MyCaloTower.h"

// ZDC
#include "UATree/UADataFormat/src/MyZDCHit.h"
#include "UATree/UADataFormat/src/MyZDCDigi.h"
#include "UATree/UADataFormat/src/MyZDCInfo.h"

// FSC
#include "UATree/UADataFormat/src/MyFSCHit.h"
#include "UATree/UADataFormat/src/MyFSCDigi.h"
#include "UATree/UADataFormat/src/MyFSCInfo.h"

using namespace std;
using namespace edm;
using namespace reco;

typedef ParameterSet PSet;

//
// Class declaration
//

class UABaseTree : public EDAnalyzer {
   public:
      explicit UABaseTree(const ParameterSet&);
      ~UABaseTree();

   private:
      virtual void beginJob() ;
      virtual void beginRun(Run const &, EventSetup const&) ;
      virtual void analyze(const Event&, const EventSetup&);
      virtual void endJob() ;

      // --------------------   Getters   --------------------
      virtual void GetAll(          const Event& , const EventSetup& );
      
      virtual void GetBeamSpot(     const Event& , const string& , MyBeamSpot& );
      virtual void GetAllBeamSpots( const Event& );
      virtual void GetEvtId(        const Event& );
      virtual void GetFwdGap(       const Event& );
      
      virtual void GetGenKin(       const Event& );
      virtual void GetGenPart(      const Event& , const EventSetup& ); 
      virtual void FillGenPart(     const GenParticle& , MyGenPart& );
      virtual void GetPUSumInfo(    const Event& );
      
      virtual void GetHLTrig(       const Event& , const EventSetup& );
      virtual void GetL1Trig(       const Event& , const EventSetup& );
      virtual void GetL1TrigOld(    const Event& );
      
      virtual void GetMITEvtSel(    const Event& );
      
      virtual void GetRecoTracks(   const Event& , const string , vector<MyTracks>& );
      virtual void GetAllTracks(    const Event& );
      virtual void FillTrack(       const Track& , MyTracks& , Double_t = MASS_PI );
      virtual void GetRecoVertex(   const Event& , const string , vector<MyVertex>& ); 
      virtual void GetAllVertices(  const Event& ); 
      
      template <class T,class U>
      void FillJetCorrections(      const Event& , const EventSetup& , const vector<T>& , const vector<string>& , vector<U>& );
      virtual void GetRecoCaloJets( const Event& , const EventSetup& , const PSet& , vector<MyCaloJet>& );
      virtual void GetAllCaloJets(  const Event& , const EventSetup& );
      virtual void GetRecoPFJets(   const Event& , const EventSetup& , const PSet& , vector<MyPFJet>& );
      virtual void GetAllPFJets(    const Event& , const EventSetup& );
      
      virtual void GetGenJets(      const Event& , const InputTag& , vector<MyGenJet>& );
      virtual void GetAllGenJets(   const Event& );
     
      template <class T>
      void FillBasicJet(            const vector<T>& , vector<MyBaseJet>& ); 
      virtual void GetBasicJet    ( const Event& , const InputTag& , vector<MyBaseJet>& );
      virtual void GetAllBasicJets( const Event& );

      virtual void GetRecoTrackJets( const Event& , const EventSetup& , const PSet& ,  vector<MyTrackJet>& );
      virtual void GetAllTrackJets ( const Event& , const EventSetup& );
      //      virtual int  GetVertexId     ( const Vertex& );

      virtual void GetCastorRecHit( const Event& ); 
      virtual void GetCastorJet(    const Event& ); 
      virtual void GetCastorDigi(   const Event& , const EventSetup& ); 
      virtual void GetCentralDiJet( const vector<MyJet*>& , const string , MyDiJet& ); 
      
      virtual void GetRecoElectron( const Event& , const InputTag& , vector<MyElectron>& ); 
      virtual void GetAllElectrons( const Event& ); 
      virtual void GetRecoMuon(     const Event& , const InputTag& , vector<MyMuon>& ); 
      virtual void GetAllMuons(     const Event& ); 
      virtual void GetRecoPFCand(   const Event& , const InputTag& , vector<MyPFCand>& ); 
      virtual void GetAllPFCands(   const Event& ); 
      
      virtual void GetMET(          const Event& , const string& , vector<MyMet>& );
      virtual void GetAllMETs(      const Event& );
      template <class T>
      void FillAllMET(              const vector<T>& , vector<MyMet>& );
      
      virtual void GetCaloTower(    const Event& ); 

      virtual void GetZDCInfo( const Event&, const EventSetup& ); 

      virtual void GetFSCInfo( const Event&, const EventSetup& ); 

      // --------------------   Get All Parameters   --------------------
      virtual void GetParameters( const ParameterSet& );
     
      // --------------------   Init All Branches   --------------------
      virtual void Init();
      
      // --------------------      Filter Event      --------------------
      virtual Bool_t FilterEvents();
      
      // --------------------   Other Functions   --------------------
      bool hasFired(const std::string& , const TriggerNames& trigNames, const TriggerResults& ) const;
      Bool_t GetLooseCaloJetId(const MyCaloJet& , const string& );
      Bool_t GetTightCaloJetId(const MyCaloJet& , const string& );
      Bool_t GetLoosePFJetId(  const MyPFJet&   , const string& );
      Bool_t GetTightPFJetId(  const MyPFJet&   , const string& );
      const string GetBranchName(InputTag& , Bool_t = 1);
      const string GetCollName(const string&);
      const string GetColl(const string&);
      const InputTag GetCollInputTag(const string&);

      // -------------------------------------------------------------------
      // --------------------   Vars From Config File   --------------------
      
      vector<InputTag> beamspots_ ;
      Bool_t           storeEvtId_;
      InputTag         calotower_ ;
      Bool_t           storeFwdGap_;
      InputTag         hepmc_ ;
      InputTag         genpart_ ;
      Bool_t           storeGenKin_;
      Bool_t           storeGenPart_;
      InputTag         pusuminfo_;
      Bool_t           storePUSumInfo_;
      vector<string>   hlt_paths_;
      Bool_t           storeL1Trig_;
      Bool_t           storeL1TrigOld_;
      Bool_t           storeMITEvtSel_;
      vector<InputTag> tracks_ ;
      vector<InputTag> vertices_ ;
      vector<PSet>     vcalojets_;
      vector<PSet>     vpfjets_;
      vector<PSet>     vtrackjets_;
      vector<InputTag> genjets_;
      vector<InputTag> basicjets_;
      vector<InputTag> trackjets_;
      string           vtxcoll_for_trackjets_;
      vector<InputTag> electrons_;
      vector<InputTag> muons_;
      vector<InputTag> pfcands_;
      vector<InputTag> mets_;
      InputTag         calotowercoll_;
      Bool_t           storeCaloObjects_;
      
      // Castor Stuff
      InputTag         castorrechits_;
      InputTag         castorjets_;
      InputTag         castorjetid_;
      InputTag         castordigis_;
      
      // ZDC
      Bool_t           storeZDCHits_;
      Bool_t           storeZDCDigis_;
      Bool_t           storeZDCInfo_;
      InputTag         zdcrechits_;
      InputTag         zdcdigis_;

      // FSC
      Bool_t           storeFSCHits_;
      Bool_t           storeFSCDigis_;
      Bool_t           storeFSCInfo_;
      InputTag         fscrechits_;
      InputTag         fscdigis_;

      //for fwdGap
      double energyThresholdHB_ ;
      double energyThresholdHE_ ;
      double energyThresholdHF_ ;
      double energyThresholdEB_ ;
      double energyThresholdEE_ ;
      
      //for genPart
      Bool_t saveMothersAndDaughters_;
      Bool_t onlyStableGenPart_;
      Bool_t onlyChargedGenPart_;
      Bool_t enableGenMetFromGenPart_;
      Bool_t saveGenPartsInDifferentColls_;
      
      //for PFJets
      Bool_t storeTracksInPFJets_;
      
      //for jet ID
      PSet ParaSetLooseCaloJetID_;
      PSet ParaSetTightCaloJetID_;
      PSet ParaSetLoosePFJetID_;
      PSet ParaSetTightPFJetID_;
      
      //for DiJets
      Double_t   jetPtCut_;
      Double_t   jetEtaCut_;
     
      Bool_t     filterEvents_; 
      
      string outputfilename_ ;

      // ----------------------------------------------------------
      // --------------------   Tree Content   --------------------
     
      map<string,MyBeamSpot>        allBeamSpots;
      MyEvtId                       evtId;
      MyFwdGap                      fwdGap;
      
      MyGenKin                      genKin;
      MyPart                        genMetfromGenPartst1;
      MyPart                        genMetfromGenPartst3;
      vector<MyGenPart>             genPart;
      vector<MyGenPart>             genElec;
      vector<MyGenPart>             genMu;
      vector<MyGenPart>             genNu;
      MySimVertex                   simVertex;
      MyPUSumInfo                   pusuminfo;      

      MyL1Trig                      L1Trig;
      MyL1TrigOld                   L1TrigOld; 
      MyHLTrig                      HLTrig;
      
      MyMITEvtSel                   MITEvtSel;
      
      map<string,vector<MyTracks> > allTracks;
      map<string,vector<MyVertex> > allVertices;

      map<string,vector<MyCaloJet> > allCaloJets;
      map<string,vector<MyPFJet> >   allPFJets;
      
      map<string,vector<MyGenJet> >   allGenJets;
      map<string,vector<MyBaseJet> >  allBasicJets;
      map<string,vector<MyTrackJet> > allTrackJets;
 
      vector<MyCastorRecHit>        castorRecHits;
      vector<MyCastorJet>           castorJets;
      vector<MyCastorDigi>          castorDigis;
      map<string , MyDiJet>         allDiJets;

      map<string,vector<MyElectron> > allElectrons;
      map<string,vector<MyMuon> >     allMuons;
      map<string,vector<MyPFCand> >   allPFCands;

      map<string,vector<MyMet> >      allMETs;

      vector<MyCaloTower>             caloTowers;

      // ZDC
      vector<MyZDCHit>  zdcHits;
      vector<MyZDCDigi> zdcDigis;
      MyZDCInfo         zdcInfo;

      // FSC
      vector<MyFSCHit>  fscHits;
      vector<MyFSCDigi> fscDigis;
      MyFSCInfo         fscInfo;

      // -------------------------------------------------------
      // --------------------   Vertex Id   --------------------
      Int_t vtxid;
      vector<math::XYZPoint> vtxid_xyz;

      // map<int,string> HLT_map;

      // --------------------   Needed For HLT   --------------------
      bool isValidHltConfig_;
      HLTConfigProvider hltConfig;
      
      // --------------------   File & Tree   --------------------
      TFile*   fout;
      TTree*   tree;
};

#include "TemplateFunctions_jets.h"

#endif

