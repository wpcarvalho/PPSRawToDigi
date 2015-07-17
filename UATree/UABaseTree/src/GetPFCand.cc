// Description: Function to retrieve PFCandidates

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"


Bool_t PFCandDebug = false;

void UABaseTree::GetRecoPFCand(const edm::Event& iEvent , const InputTag& pfcandcoll_ , vector<MyPFCand>& PFCandVector){

   PFCandVector.clear();

   Handle<PFCandidateCollection> PFCands;
   try {
     iEvent.getByLabel(pfcandcoll_,PFCands);
   }
   catch ( ... ) {
     cout << "[UABaseTree::GetRecoPFCand] Can't find the collection " << pfcandcoll_ << endl;
   }
   
   PFCandVector.assign( PFCands->size() , MyPFCand() );

   Int_t i = 0;
   for (PFCandidateCollection::const_iterator iPFCand = PFCands->begin() ; iPFCand != PFCands->end() ; ++iPFCand , ++i) {
     
     //PFCandVector[i].SetPxPyPzE(iPFCand->px() , iPFCand->py() , iPFCand->pz() , sqrt(iPFCand->momentum().mag2()+MASS_MU*MASS_MU));
     PFCandVector[i].SetPxPyPzE(iPFCand->px() , iPFCand->py() , iPFCand->pz() , iPFCand->energy() );
     PFCandVector[i].charge = iPFCand->charge();
     PFCandVector[i].particleId = static_cast<MyPFCand::ParticleType>(iPFCand->particleId());

     if(PFCandDebug) PFCandVector[i].Print();


   } // end for PFCandCollection 

}


void UABaseTree::GetAllPFCands( const edm::Event& iEvent ){
  for(vector<InputTag>::iterator icoll = pfcands_.begin() ; icoll!= pfcands_.end() ; ++icoll)
    this->GetRecoPFCand(iEvent , *icoll , allPFCands[icoll->label()] );
}
