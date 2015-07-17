#include "DataFormats/JetReco/interface/Jet.h"

//PF Jets
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

//Calo Jets
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

//Track Jets
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"

//Gen Jets
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"

bool BasicJetDebug = false;

void UABaseTree::GetBasicJet(const edm::Event& iEvent , const InputTag& JetCol , vector<MyBaseJet>& JetVector){

  try{
    Handle<TrackJetCollection> jets;
    iEvent.getByLabel(JetCol , jets);
    FillBasicJet(*jets , JetVector );
  } catch(...){;}

  try{
    Handle<PFJetCollection> jets;
    iEvent.getByLabel(JetCol , jets);
    FillBasicJet(*jets , JetVector );
  } catch(...){;}

  try{
    Handle<CaloJetCollection> jets;
    iEvent.getByLabel(JetCol , jets);
    FillBasicJet(*jets , JetVector );
  } catch(...){;}


  try{
    Handle<GenJetCollection> jets;
    iEvent.getByLabel(JetCol , jets);
    FillBasicJet(*jets , JetVector );
  } catch(...){;}

}

void UABaseTree::GetAllBasicJets(const edm::Event& iEvent){

  for(vector<InputTag>::iterator icoll = basicjets_.begin() ; icoll != basicjets_.end() ; ++icoll){
    this->GetBasicJet( iEvent , *icoll , allBasicJets[icoll->label()] );
  }
}


