// Description: Function to retrieve the Parameters from the python file

// UABaseTree Analysis class decleration
#include "UATree/UABaseTree/interface/UABaseTree.h"

bool ParametersDebug = false;


void UABaseTree::GetParameters(const edm::ParameterSet& iConfig){

   beamspots_      = iConfig.getUntrackedParameter<vector<InputTag> >("beamspots",vector<InputTag>(1,InputTag("offlineBeamSpot")));
   storeEvtId_     = iConfig.getUntrackedParameter<bool>("storeEvtId",false);
   calotower_      = iConfig.getUntrackedParameter<InputTag>("calotower",InputTag("towerMaker"));
   storeFwdGap_    = iConfig.getUntrackedParameter<bool>("storeFwdGap",false);
   hepmc_          = iConfig.getUntrackedParameter<InputTag>("hepmc",InputTag("generator"));
   genpart_        = iConfig.getUntrackedParameter<InputTag>("genpart",InputTag("genParticles"));
   storeGenKin_    = iConfig.getUntrackedParameter<bool>("storeGenKin",false);
   storeGenPart_   = iConfig.getUntrackedParameter<bool>("storeGenPart",false);
   pusuminfo_      = iConfig.getUntrackedParameter<InputTag>("pusuminfo",InputTag("addPileupInfo"));
   storePUSumInfo_ = iConfig.getUntrackedParameter<bool>("storePUSumInfo",false);
   hlt_paths_      = iConfig.getUntrackedParameter<vector<string> >("hlt_paths",vector<string>());
   storeL1Trig_    = iConfig.getUntrackedParameter<bool>("storeL1Trig",false);
   storeL1TrigOld_ = iConfig.getUntrackedParameter<bool>("storeL1TrigOld",false);
   storeMITEvtSel_ = iConfig.getUntrackedParameter<bool>("storeMITEvtSel",false);
   tracks_         = iConfig.getUntrackedParameter<vector<InputTag> >("tracks",vector<InputTag>());
   vertices_       = iConfig.getUntrackedParameter<vector<InputTag> >("vertices",vector<InputTag>());
   vcalojets_      = iConfig.getUntrackedParameter<vector<PSet> >("vcalojets",vector<PSet>());
   vpfjets_        = iConfig.getUntrackedParameter<vector<PSet> >("vpfjets",vector<PSet>());
   genjets_        = iConfig.getUntrackedParameter<vector<InputTag> >("genjets",vector<InputTag>());
   basicjets_      = iConfig.getUntrackedParameter<vector<InputTag> >("basicjets",vector<InputTag>());

   trackjets_      = iConfig.getUntrackedParameter<vector<InputTag> >("trackjets",vector<InputTag>());
   vtxcoll_for_trackjets_ = iConfig.getUntrackedParameter<string>("vtxcoll_for_trackjets","");
   
   castorrechits_  = iConfig.getUntrackedParameter<InputTag>("castorrechits",InputTag());
   castorjets_     = iConfig.getUntrackedParameter<InputTag>("castorjets",InputTag());
   castorjetid_    = iConfig.getUntrackedParameter<InputTag>("castorjetid",InputTag());
   castordigis_    = iConfig.getUntrackedParameter<InputTag>("castordigis",InputTag());
   
   electrons_      = iConfig.getUntrackedParameter<vector<InputTag> >("electrons",vector<InputTag>());
   muons_          = iConfig.getUntrackedParameter<vector<InputTag> >("muons",vector<InputTag>());
   pfcands_        = iConfig.getUntrackedParameter<vector<InputTag> >("pfcands",vector<InputTag>());
   
   mets_           = iConfig.getUntrackedParameter<vector<InputTag> >("mets",vector<InputTag>());

   calotowercoll_     = iConfig.getUntrackedParameter<InputTag>("calotowercoll",InputTag("towerMaker"));
   storeCaloObjects_  = iConfig.getUntrackedParameter<bool>("storeCaloObjects",false);
   
   // ZDC
   storeZDCInfo_  = iConfig.getUntrackedParameter<bool>("storeZDCInfo",false);
   storeZDCHits_  = iConfig.getUntrackedParameter<bool>("storeZDCHits",false);
   storeZDCDigis_ = iConfig.getUntrackedParameter<bool>("storeZDCDigis",false);
   zdcrechits_    = iConfig.getUntrackedParameter<InputTag>("zdcrechits",InputTag());
   zdcdigis_      = iConfig.getUntrackedParameter<InputTag>("zdcdigis",InputTag());

   // FSC
   storeFSCInfo_  = iConfig.getUntrackedParameter<bool>("storeFSCInfo",false);
   storeFSCHits_  = iConfig.getUntrackedParameter<bool>("storeFSCHits",false);
   storeFSCDigis_ = iConfig.getUntrackedParameter<bool>("storeFSCDigis",false);
   fscrechits_    = iConfig.getUntrackedParameter<InputTag>("fscrechits",InputTag());
   fscdigis_      = iConfig.getUntrackedParameter<InputTag>("fscdigis",InputTag());

   //Specific for fwdGap
   energyThresholdHB_ = iConfig.getUntrackedParameter<double>("EnergyThresholdHB",1.5) ;
   energyThresholdHE_ = iConfig.getUntrackedParameter<double>("EnergyThresholdHE",2.0) ;
   energyThresholdHF_ = iConfig.getUntrackedParameter<double>("EnergyThresholdHF",4.0) ;
   energyThresholdEB_ = iConfig.getUntrackedParameter<double>("EnergyThresholdEB",1.5) ;
   energyThresholdEE_ = iConfig.getUntrackedParameter<double>("EnergyThresholdEE",2.5) ;

   //Specific for genPart
   saveMothersAndDaughters_      = iConfig.getUntrackedParameter<bool>("saveMothersAndDaughters",false);
   onlyStableGenPart_            = iConfig.getUntrackedParameter<bool>("onlyStableGenPart",false);
   onlyChargedGenPart_           = iConfig.getUntrackedParameter<bool>("onlyChargedGenPart_",false);
   enableGenMetFromGenPart_      = iConfig.getUntrackedParameter<bool>("enableGenMetFromGenPart",false);
   saveGenPartsInDifferentColls_ = iConfig.getUntrackedParameter<bool>("saveGenPartsInDifferentColls",false);


   //Specific for PFJets
   storeTracksInPFJets_        = iConfig.getUntrackedParameter<bool>("storeTracksInPFJets",false);


   //Specific to DiJets
   jetPtCut_  = iConfig.getUntrackedParameter<double>("jetPtCut"   , 8.0);
   jetEtaCut_ = iConfig.getUntrackedParameter<double>("jetEtaCut" , 2.5);

   //Switch Event Filter On/Off
   filterEvents_ = iConfig.getUntrackedParameter<bool>("filterEvents",false);

   // Get output filename
   outputfilename_ = iConfig.getUntrackedParameter<string>("outputfilename","");

}
