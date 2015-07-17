//-- Description: Function to retrieve PF Jet information

//--DataFormats
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"


#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"

bool TrackJetDebug = false;

void UABaseTree::GetRecoTrackJets(const Event& iEvent , const EventSetup& iSetup , const PSet& trackjets_ , vector<MyTrackJet>& JetVector){

  JetVector.clear();
  
  InputTag jetcoll_           = trackjets_.getUntrackedParameter<InputTag>("jetcoll",InputTag());
  GetBranchName(jetcoll_);
  vector<string> corrections_ = trackjets_.getUntrackedParameter<vector<string> >("corrections",vector<string>());
    
  Handle<TrackJetCollection> raw;
  iEvent.getByLabel(jetcoll_.label(),raw);
  
  //Initialize vector with number of jets
  JetVector.assign(raw->size() , MyTrackJet());
    
  //Filling mapjet with corrections
  FillJetCorrections( iEvent , iSetup , *raw , corrections_ , JetVector);
    
  //filling raw collection
  Int_t i = 0;
  for (TrackJetCollection::const_iterator jet = raw->begin(); jet != raw->end(); ++jet , ++i){
    JetVector[i].mapjet[jetcoll_.label()].SetPxPyPzE(jet->px() , jet->py() , jet->pz() , jet->energy() );
    if(TrackJetDebug) JetVector[i].Print();
  }
  
  //doing the DiJet if needed
  string dijetcoll_ = trackjets_.getUntrackedParameter<string>("dijetcoll","");
  if(JetVector.size() > 1.){
    if(JetVector[0].mapjet.find(GetCollName(dijetcoll_)) != JetVector[0].mapjet.end()){
  
      //making the vector of MyJet*
      vector<MyJet*> tmpJetVector(JetVector.size() , NULL);
      for(unsigned int i = 0 ; i < JetVector.size() ; ++i)
        tmpJetVector[i] = (MyJet*) &(JetVector[i]);
  
      //Running the DiJet algo
      GetCentralDiJet(tmpJetVector , GetCollName(dijetcoll_),  allDiJets[GetColl(dijetcoll_)+"DiJet"] );
    }
  }

}

void UABaseTree::GetAllTrackJets(const edm::Event& iEvent , const edm::EventSetup& iSetup){
  
  if(TrackJetDebug) cout << "Number of Tracket collections " << vtrackjets_.size() << endl;
  InputTag jetcoll_;
  for(vector<PSet>::iterator it = vtrackjets_.begin() ; it != vtrackjets_.end() ; ++it){
    jetcoll_ = it->getUntrackedParameter<InputTag>("jetcoll",InputTag());
    GetBranchName(jetcoll_);
    if(jetcoll_.label().size() > 0)
      this->GetRecoTrackJets(iEvent , iSetup , *it , this->allTrackJets[jetcoll_.label()] );  
  }
}


// -------------------------------------------------------------------------------------------------------------------------------------------




