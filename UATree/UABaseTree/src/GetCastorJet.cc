//-- Description: Function to retrieve castor jet information

//-- Castor Jet
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/CastorJetID.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"

bool CastorJetDebug = false;

void UABaseTree::GetCastorJet(const edm::Event& iEvent) {
  
  castorJets.clear();
  MyCastorJet mycastorjet;
  
  edm::Handle<edm::View< reco::BasicJet > > basicjetcoll;  //-- uncorrected jets
  iEvent.getByLabel(castorjets_,basicjetcoll);

  edm::Handle<reco::CastorJetIDValueMap> jetIdMap;
  iEvent.getByLabel(castorjetid_,jetIdMap);

  if(basicjetcoll.isValid()) {

    for(edm::View<reco::BasicJet>::const_iterator ibegin = basicjetcoll->begin(), iend = basicjetcoll->end(), ijet = ibegin; ijet != iend; ++ijet) {

      unsigned int idx = ijet - ibegin;
      const BasicJet & basicjet = (*basicjetcoll)[idx];

      mycastorjet.energy = basicjet.energy();
      mycastorjet.pt = basicjet.pt();
      mycastorjet.eta = basicjet.eta();
      mycastorjet.phi = basicjet.phi();

      edm::RefToBase<reco::BasicJet> jetRef = basicjetcoll->refAt(idx);
      reco::CastorJetID const & jetId = (*jetIdMap)[jetRef];

      mycastorjet.fem = jetId.fem;
      mycastorjet.eem = jetId.emEnergy;
      mycastorjet.ehad = jetId.hadEnergy;

      mycastorjet.width = jetId.width;
      mycastorjet.depth = jetId.depth;
      mycastorjet.fhot = jetId.fhot;
      mycastorjet.sigmaz = jetId.sigmaz;

      mycastorjet.ntower = jetId.nTowers;

      castorJets.push_back(mycastorjet);
    
      if (CastorJetDebug) mycastorjet.Print();   
    }
  }
}
