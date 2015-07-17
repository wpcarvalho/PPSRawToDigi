// Description: Function to retrieve Generated particles from HEPEVT

// system include files
#include <memory>
   
// DataFormats
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

// UABaseTree Analysis class declaration
#include "UATree/UABaseTree/interface/UABaseTree.h"

bool GenPartDebug = false ;

void UABaseTree::GetGenPart(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Clear GenPart List
  // (evt->GenPart).clear();
  genPart.clear();

  // Handle to access PDG data from iSetup
  ESHandle <ParticleDataTable> pdt;
  iSetup.getData( pdt );

  // Handle to access GenParticleCollection
  Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel(genpart_, genParticles);
  if (!genParticles.isValid()) {
    cerr << "[UABaseTree::GetGenPart] Error: non valid GenParticleCollection " << endl;
    return;
  }

  // List for Daugther/Mother Search
  vector<const reco::Candidate *> cands;
  vector<const Candidate *>::const_iterator found = cands.begin();
  if(saveMothersAndDaughters_){
    for(GenParticleCollection::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p)
      cands.push_back(&*p);
  }

  //-- position of the simulated primary vertex
  math::XYZPoint PosPVsim = (*genParticles)[2].vertex();
  simVertex.x = PosPVsim.X();
  simVertex.y = PosPVsim.Y();
  simVertex.z = PosPVsim.Z();
  if(GenPartDebug) simVertex.Print();
  
  //-- for the GenMEt
  double met1Px = 0 ;
  double met1Py = 0 ;
  double met1Pz = 0 ;
  double met1E  = 0 ;

  double met3Px = 0 ;
  double met3Py = 0 ;
  double met3Pz = 0 ;
  double met3E  = 0 ;
  
  

  // Loop on generated particle
  if ( GenPartDebug )
    cout << "GenPart # : " << genParticles->size() << endl;
  
  MyGenPart mygenpart;
  for(GenParticleCollection::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) {
    
    mygenpart.Reset();
    
    int st  = p->status();
    int pid = p->pdgId();

    this->FillGenPart(*p , mygenpart);
    
    // Mother Daughter relations
    if(saveMothersAndDaughters_){
      int nMo = p->numberOfMothers();
      int nDa = p->numberOfDaughters();

      found = find(cands.begin(), cands.end(), p->mother(0));
      if(found != cands.end()) mygenpart.mo1 = found - cands.begin() ;
 
      found = find(cands.begin(), cands.end(), p->mother(nMo-1));
      if(found != cands.end()) mygenpart.mo2 = found - cands.begin() ;

      found = find(cands.begin(), cands.end(), p->daughter(0));
      if(found != cands.end()) mygenpart.da1 = found - cands.begin() ;
   
      found = find(cands.begin(), cands.end(), p->daughter(nDa-1));
      if(found != cands.end()) mygenpart.da2 = found - cands.begin() ;
    }


    Bool_t store = true;
    if( onlyStableGenPart_  && st!=1)                      store = false;
    if( onlyChargedGenPart_ && fabs(mygenpart.charge) < 1) store = false;
    if( saveGenPartsInDifferentColls_ && st!=3 )           store = false;


    if(store) genPart.push_back(mygenpart);
    if(GenPartDebug && store) mygenpart.Print();

    //Saves stable electrons , muons and neutrinos in different collections
    if(saveGenPartsInDifferentColls_ && st==1){
       if(fabs(mygenpart.pdgId)==11)
         genElec.push_back(mygenpart);

       if(fabs(mygenpart.pdgId)==13)
         genMu.push_back(mygenpart);

       if( (fabs(mygenpart.pdgId)==12) || (fabs(mygenpart.pdgId)==14) )
         genNu.push_back(mygenpart);
    }
    
    
    
    
    // Gen Met: From GenPart Neutrinos ( no SUSY !!! )
    if(enableGenMetFromGenPart_){
      if( ( abs(pid)==12 || abs(pid)==14 || abs(pid)==16 ) ){
        if ( st == 1 ){
          met1Px += p->px();
          met1Py += p->py();
          met1Pz += p->pz();
          met1E  += p->energy();
        } 
        else if ( st == 3 ){
          met3Px += p->px();
          met3Py += p->py();
          met3Pz += p->pz();
          met3E  += p->energy();
        }
      }
    }

  }//end of loop over genParticles
  
  //storing genMet
  if(enableGenMetFromGenPart_){
    genMetfromGenPartst1.SetPxPyPzE( met1Px , met1Py , met1Pz , met1E ) ;
    genMetfromGenPartst3.SetPxPyPzE( met3Px , met3Py , met3Pz , met3E ) ;
  }

}



void UABaseTree::FillGenPart(const GenParticle& ingp , MyGenPart& outgp){

  //Filling inherited from MyPart
  outgp.SetPxPyPzE( ingp.px() , ingp.py() , ingp.pz() , ingp.energy() );
  outgp.charge  = ingp.charge();

  // Extra properties
  outgp.pdgId   = ingp.pdgId();
  //outgp.name    = (pdt->particle(ingp.pdgId()))->name(); //FIXME : crashes the code ...
  outgp.status  = ingp.status();
}
