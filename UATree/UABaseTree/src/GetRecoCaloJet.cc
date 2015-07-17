//-- Description: Function to retrieve Calo Jet information

//--DataFormats
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"

bool CaloJetDebug = false;

void UABaseTree::GetRecoCaloJets(const edm::Event& iEvent, const edm::EventSetup& iSetup, const PSet& calojets_ , vector<MyCaloJet>& JetVector){

  JetVector.clear();
  
  InputTag jetcoll_           = calojets_.getUntrackedParameter<InputTag>("jetcoll",InputTag());
  GetBranchName(jetcoll_);
  vector<string> corrections_ = calojets_.getUntrackedParameter<vector<string> >("corrections",vector<string>());
  InputTag calojetid_  = calojets_.getUntrackedParameter<InputTag>("calojetid",InputTag());
    
  Handle<CaloJetCollection> raw;
  iEvent.getByLabel(jetcoll_.label(),raw);
  
  //Initialize vector with number of jets
  JetVector.assign(raw->size(),MyCaloJet());
  
  //Specific to calojets
  Handle<reco::JetIDValueMap> CaloJetId;
  iEvent.getByLabel(calojetid_,CaloJetId);
  JetIDValueMap jetid = *CaloJetId;

  //Filling mapjet with corrections
  this->FillJetCorrections( iEvent , iSetup , *raw , corrections_ , JetVector);
  
  //filling raw collection
  Int_t i = 0;
  for (CaloJetCollection::const_iterator jet = raw->begin(); jet != raw->end(); ++jet , ++i){
    JetVector[i].mapjet[jetcoll_.label()].SetPxPyPzE(jet->px() , jet->py() , jet->pz() , jet->energy() );
    JetVector[i].fem = jet->emEnergyFraction();
    JetVector[i].eem_EB = jet->emEnergyInEB();
    JetVector[i].eem_EE = jet->emEnergyInEE();
    JetVector[i].eem_HF = jet->emEnergyInHF();

    JetVector[i].fhad = jet->energyFractionHadronic();
    JetVector[i].ehad_HB = jet->hadEnergyInHB();
    JetVector[i].ehad_HE = jet->hadEnergyInHE();
    JetVector[i].ehad_HF = jet->hadEnergyInHF();
    JetVector[i].ehad_HO = jet->hadEnergyInHO();

    JetVector[i].n60 = jet->n60();
    JetVector[i].n90 = jet->n90();

    JetVector[i].emax_ecal = jet->maxEInEmTowers();
    JetVector[i].emax_hcal = jet->maxEInHadTowers();

    CaloJetRef jetref(raw, jet - raw->begin() );

    JetVector[i].n90hits = jetid[jetref].n90Hits;
    JetVector[i].HPD = jetid[jetref].fHPD;
    JetVector[i].RBX = jetid[jetref].fRBX;
    JetVector[i].sigma_eta = sqrt(jet->etaetaMoment());
    JetVector[i].sigma_phi = sqrt(jet->phiphiMoment());

    //-- number of constituents (PFObject for PFjet, CaloTower for Calojet)
    JetVector[i].nconstituent = jet->getCaloConstituents().size();
    
    //-- jet ID
    JetVector[i].LooseJetId = GetLooseCaloJetId(JetVector[i],jetcoll_.label());
    JetVector[i].TightJetId = GetTightCaloJetId(JetVector[i],jetcoll_.label());

    if(CaloJetDebug) JetVector[i].Print();

  }

  
  //doing the DiJet if needed
  string dijetcoll_ = calojets_.getUntrackedParameter<string>("dijetcoll","");
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


void UABaseTree::GetAllCaloJets(const edm::Event& iEvent , const edm::EventSetup& iSetup){
  
  InputTag jetcoll_;
  for(vector<PSet>::iterator it = vcalojets_.begin() ; it != vcalojets_.end() ; ++it){
    jetcoll_ = it->getUntrackedParameter<InputTag>("jetcoll",InputTag());
    GetBranchName(jetcoll_);
    if(jetcoll_.label().size() > 0)
      this->GetRecoCaloJets(iEvent , iSetup , *it , this->allCaloJets[jetcoll_.label()] );  
  }
}


// -------------------------------------------------------------------------------------------------------------------------------------------


Bool_t UABaseTree::GetLooseCaloJetId(const MyCaloJet& myjet , const string& raw_coll) {

  //-- parameters for the selection

  Double_t HPD_CutUp;
  Int_t n90hits_CutLow;
  Double_t fem_CutLow;

  //-- retrieve the parameters

  HPD_CutUp = ParaSetLooseCaloJetID_.getUntrackedParameter<double>("HPD_CutUp",0);
  n90hits_CutLow = ParaSetLooseCaloJetID_.getUntrackedParameter<int>("n90hits_CutLow",0);
  fem_CutLow = ParaSetLooseCaloJetID_.getUntrackedParameter<double>("fem_CutLow",0);

  //-- loose selection for Calo jet

  Bool_t accept = false;

  //-- loose jet id
  if(myjet.HPD < HPD_CutUp && myjet.n90hits > n90hits_CutLow) {
    //cout << &((myjet.mapjet)[string("gf")]);
    if(fabs(myjet.mapjet.find(raw_coll)->second.Eta()) < 2.4) { //FIXME: depends on raw eta here !!
      if(myjet.fem > fem_CutLow) accept = true;
    }

    //-- outside tracker acceptance, jet accepted without fem
    else accept = true;
  }

  return(accept);
}


Bool_t UABaseTree::GetTightCaloJetId(const MyCaloJet& myjet , const string& raw_coll) {

  //-- parameters for the selection

  Double_t HPD_CutUp;
  Int_t    n90hits_CutLow;
  Double_t RBX_CutUp;
  Double_t sigma_eta_CutLow;
  Double_t sigma_phi_CutLow;
  Double_t fem_CutLow;

  //-- retrieve the parameters

  HPD_CutUp = ParaSetTightCaloJetID_.getUntrackedParameter<double>("HPD_CutUp",0);
  n90hits_CutLow = ParaSetTightCaloJetID_.getUntrackedParameter<int>("n90hits_CutLow",0);
  RBX_CutUp = ParaSetTightCaloJetID_.getUntrackedParameter<double>("RBX_CutUp",0);
  sigma_eta_CutLow = ParaSetTightCaloJetID_.getUntrackedParameter<double>("sigma_eta_CutLow",0);
  sigma_phi_CutLow = ParaSetTightCaloJetID_.getUntrackedParameter<double>("sigma_phi_CutLow",0);
  fem_CutLow = ParaSetTightCaloJetID_.getUntrackedParameter<double>("fem_CutLow",0);

  //-- tight selection for Calo jet

  Bool_t accept = false;

  //-- tight jet id
  if(myjet.HPD < HPD_CutUp && myjet.n90hits > n90hits_CutLow && myjet.RBX < RBX_CutUp
     && myjet.sigma_eta > sigma_eta_CutLow && myjet.sigma_phi > sigma_phi_CutLow) {

    //-- inside tracker acceptance depends on fem
    if(fabs(myjet.mapjet.find(raw_coll)->second.Eta()) < 2.4) { //FIXME: depends on raw eta here !!
      if(myjet.fem > fem_CutLow) accept = true;
    }

    //-- outside tracker acceptance, jet accepted without fem
    else accept = true;
  }

  return(accept);
}

