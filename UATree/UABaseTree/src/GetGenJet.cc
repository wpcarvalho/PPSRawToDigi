#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

// UAHiggsTree UAHiggs class decleration
#include "UATree/UABaseTree/interface/UABaseTree.h"


bool GenJetDebug = false;

void UABaseTree::GetGenJets(const edm::Event& iEvent , const InputTag& genjetcoll_ , vector<MyGenJet>& JetVector ){

  JetVector.clear();
  
  //FIXME : needed only for name which is for now commented
  // Handle to access PDG data from iSetup
  //ESHandle <ParticleDataTable> pdt;
  //iSetup.getData( pdt );


  // get gen jet collection
  Handle<GenJetCollection> genjets;
  iEvent.getByLabel(genjetcoll_ , genjets);
  
  JetVector.assign(genjets->size() , MyGenJet());
  
  int ijet = 0;  
  for(GenJetCollection::const_iterator genjet = genjets->begin() ; genjet != genjets->end() ; ++genjet , ++ijet){ 
  
    JetVector[ijet].SetPxPyPzE( genjet->px() , genjet->py() , genjet->pz() , genjet->energy() );
    JetVector[ijet].npart = genjet->nConstituents();
	 
    vector <const GenParticle*> mcparts = genjet->getGenConstituents();
    JetVector[ijet].JetPart.assign(mcparts.size() , MyGenPart());
    
    //FIXME : doesn't fill mothers & daughters
    for (unsigned i = 0; i < mcparts.size (); i++)
      this->FillGenPart( *(mcparts[i]) , JetVector[ijet].JetPart[i]);
       
  }

}


void UABaseTree::GetAllGenJets( const edm::Event& iEvent ){
  for(vector<InputTag>::iterator icoll = genjets_.begin() ; icoll!= genjets_.end() ; ++icoll)
    this->GetGenJets(iEvent , *icoll , allGenJets[icoll->label()] );
}

