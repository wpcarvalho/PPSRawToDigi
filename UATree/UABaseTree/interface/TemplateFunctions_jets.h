//--DataFormats
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"


template <class T,class U>
void UABaseTree::FillJetCorrections(const edm::Event& iEvent , const EventSetup& iSetup , const vector<T>& raw , const vector<string>& corrections_ , vector<U>& JetVector ){

  //start looping over corrections
  const JetCorrector* Jetcorrector = NULL;
  string correction = "";
  for(unsigned int corr=0 ; corr < corrections_.size() ; ++corr){ 
  
    //Retrieving the correction name to store in mapjet
    correction = GetCollName(corrections_[corr]);

    try{  
      Jetcorrector = JetCorrector::getJetCorrector(GetColl(corrections_[corr]),iSetup);
      
      correction = GetCollName(corrections_[corr]);
      
      //Filling corrected collection
      Int_t i = 0;
      for (typename vector<T>::const_iterator jet = raw.begin(); jet != raw.end(); ++jet , ++i){
        //-- correction 
        T corrected_jet = *jet;                              //-- copy orignial jet

	//Old Way
        //Double_t jec = Jetcorrector->correction(jet->p4());  //-- calculate correction 
	//New Way
        int index = jet-raw.begin();
	edm::RefToBase<reco::Jet> jetRef(edm::Ref<vector<T> >(&raw,index));
        Double_t jec = Jetcorrector->correction(corrected_jet,jetRef,iEvent,iSetup);

        corrected_jet.scaleEnergy(jec);                      //-- apply correction
        JetVector[i].mapjet[correction].jec = jec;

        //-- uncertainty (function of the CORRECTED jet)
        //PFJetCorUnc->setJetEta(corrected_jet.eta());
        //PFJetCorUnc->setJetPt(corrected_jet.pt());
        //JetVector[i].mapjet[correction].jec_unc = PFJetCorUnc->getUncertainty(true);
        JetVector[i].mapjet[correction].SetPxPyPzE(corrected_jet.px() , corrected_jet.py() , corrected_jet.pz() , corrected_jet.energy() );
      }
    }
    catch(...){
      cout << "Please provide an ESSource for coll " << correction << endl;
    }
  }
    
}


template <class T>
void UABaseTree::FillBasicJet(const vector<T>& jetCol , vector<MyBaseJet>& JetVector){

  JetVector.assign(jetCol.size() , MyBaseJet());

  Int_t i = 0;
  for(typename vector<T>::const_iterator jet = jetCol.begin() ; jet != jetCol.end() ; ++jet , ++i)
    JetVector[i].SetPxPyPzE(jet->px() , jet->py() , jet->pz() , jet->energy() );

}
