//-- Description: Function to retrieve PF Jet information

//--DataFormats
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"

bool PFJetDebug = false;

void UABaseTree::GetRecoPFJets(const Event& iEvent , const EventSetup& iSetup , const PSet& pfjets_ , vector<MyPFJet>& JetVector){

  JetVector.clear();

  
  InputTag jetcoll_           = pfjets_.getUntrackedParameter<InputTag>("jetcoll",InputTag());
  GetBranchName(jetcoll_);
  vector<string> corrections_ = pfjets_.getUntrackedParameter<vector<string> >("corrections",vector<string>());
    
  Handle<PFJetCollection> raw;
  iEvent.getByLabel(jetcoll_.label(),raw);
  
  //Initialize vector with number of jets
  JetVector.assign(raw->size() , MyPFJet());
    
  //Filling mapjet with corrections
  FillJetCorrections( iEvent , iSetup , *raw , corrections_ , JetVector);
    
  //filling raw collection
  Int_t i = 0;
  for (PFJetCollection::const_iterator jet = raw->begin(); jet != raw->end(); ++jet , ++i){
  
    JetVector[i].mapjet[jetcoll_.label()].SetPxPyPzE(jet->px() , jet->py() , jet->pz() , jet->energy() );
    
    JetVector[i].fhad_ch = jet->chargedHadronEnergyFraction();
    JetVector[i].fhad_ne = jet->neutralHadronEnergyFraction();

    JetVector[i].fem_ch = jet->chargedEmEnergyFraction();
    JetVector[i].fem_ne = jet->neutralEmEnergyFraction();

    JetVector[i].multi_ch = jet->chargedMultiplicity();
    JetVector[i].multi_ne = jet->neutralMultiplicity();

    JetVector[i].multi_ch_had = jet->chargedHadronMultiplicity();
    JetVector[i].multi_ne_had = jet->neutralHadronMultiplicity();

    JetVector[i].multi_gamma = jet->photonMultiplicity();
    JetVector[i].multi_ele = jet->electronMultiplicity();
    JetVector[i].multi_mu = jet->muonMultiplicity();


    //-- number of constituents (PFObject for PFjet, CaloTower for Calojet)
    JetVector[i].nconstituent = jet->getPFConstituents().size();
    
    //-- number of tracks (PFjet only)
    if(jet->getTrackRefs().isAvailable()){
      reco::TrackRefVector vtrackref = jet->getTrackRefs();
      JetVector[i].ntrack = vtrackref.size();
      if(storeTracksInPFJets_){
	JetVector[i].vtracks.assign(JetVector[i].ntrack , MyTracks());
        for(unsigned int i=0; i<vtrackref.size();i++)
          if(vtrackref.at(i).isAvailable())
	    this->FillTrack(*(vtrackref.at(i)) , JetVector[i].vtracks[i]);
      }
    }

    
    //-- jet ID
    //JetVector[i].LooseJetId = GetLoosePFJetId(JetVector[i],jetcoll_.label());
    JetVector[i].TightJetId = GetTightPFJetId(JetVector[i],jetcoll_.label());
    JetVector[i].LooseJetId = JetVector[i].TightJetId; //here loose & tight have the same output ... (used for DiJet selection)

    if(PFJetDebug) JetVector[i].Print();

  }


  

    //JetCorrectionUncertainty *PFJetCorUnc(0);

    //ESHandle<JetCorrectorParametersCollection> PFJetCorParColl;
    //iSetup.get<JetCorrectionsRecord>().get(PFJetJECunc_,PFJetCorParColl);
    
    //JetCorrectorParameters const & PFJetCorPar = (*PFJetCorParColl)["Uncertainty"];
    //PFJetCorUnc = new JetCorrectionUncertainty(PFJetCorPar);

    //nPFJet = PFJetcoll->size();
  
  
 
  //doing the DiJet if needed
  string dijetcoll_ = pfjets_.getUntrackedParameter<string>("dijetcoll","");
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


void UABaseTree::GetAllPFJets(const edm::Event& iEvent , const edm::EventSetup& iSetup){
  
  if(PFJetDebug) cout << "Number of PJet collections " << vpfjets_.size() << endl;
  InputTag jetcoll_;
  for(vector<PSet>::iterator it = vpfjets_.begin() ; it != vpfjets_.end() ; ++it){
    jetcoll_ = it->getUntrackedParameter<InputTag>("jetcoll",InputTag());
    GetBranchName(jetcoll_);
    if(jetcoll_.label().size() > 0)
      this->GetRecoPFJets(iEvent , iSetup , *it , this->allPFJets[jetcoll_.label()] );  
  }
}


// -------------------------------------------------------------------------------------------------------------------------------------------


Bool_t UABaseTree::GetLoosePFJetId(const MyPFJet& myjet , const string& raw_coll) {

  //DOESN'T EXIST
  //ALWAYS RETURNS 1

  return true;
}


Bool_t UABaseTree::GetTightPFJetId(const MyPFJet& myjet , const string& raw_coll) {

  //-- parameters for the selection

  double fhad_ne_CutUp;
  double fem_ne_CutUp;

  double fhad_ch_CutLow;
  double fem_ch_CutUp;

  unsigned int nconstituent_CutLow;
  int multi_ch_CutLow;

  //-- retrieve the parameters

  fhad_ne_CutUp = ParaSetTightPFJetID_.getUntrackedParameter<double>("fhad_ne_CutUp",0);
  fem_ne_CutUp = ParaSetTightPFJetID_.getUntrackedParameter<double>("fem_ne_CutUp",0);

  fhad_ch_CutLow = ParaSetTightPFJetID_.getUntrackedParameter<double>("fhad_ch_CutLow",0);
  fem_ch_CutUp = ParaSetTightPFJetID_.getUntrackedParameter<double>("fem_ch_CutUp",0);

  nconstituent_CutLow = ParaSetTightPFJetID_.getUntrackedParameter<int>("nconstituent_CutLow",0);
  multi_ch_CutLow = ParaSetTightPFJetID_.getUntrackedParameter<int>("multi_ch_CutLow",0);

  //-- thight selection for PF jet

  Bool_t accept = false;

  //-- selection depends on f_nh, f_ne and nconstituent
  if(myjet.fhad_ne < fhad_ne_CutUp && myjet.fem_ne < fem_ne_CutUp && myjet.nconstituent > nconstituent_CutLow) {

    //-- inside tracker acceptance, selection depends on f_ce, f_ch and multi_c
    if(fabs(myjet.mapjet.find(raw_coll)->second.Eta()) < 2.4) {
      if(myjet.fem_ch < fem_ch_CutUp && myjet.fhad_ch > fhad_ch_CutLow && myjet.multi_ch > multi_ch_CutLow)
        accept = true;
    }

    //-- outside tracker acceptance, jet accepted without f_ce, f_ch and multi_c
    else accept = true;
  }

  return(accept);

}


