#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"

#include "DataFormats/Common/interface/Handle.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"

Bool_t METDebug = false;

void UABaseTree::GetMET(const edm::Event& iEvent , const string& metcoll_ , vector<MyMet>& METVector ){
  
  
  Handle<GenMETCollection> genmet;
  try{
    iEvent.getByLabel(metcoll_ , genmet);
  } catch(...){}
  
    
  Handle<PFMETCollection> pfmet;
  try{
    iEvent.getByLabel(metcoll_ , pfmet);
  } catch(...){}
    
    
  Handle<CaloMETCollection> calomet;
  try{
    iEvent.getByLabel(metcoll_ , calomet);
  } catch(...){}
    
    
  Handle<METCollection> met;
  try{
    iEvent.getByLabel(metcoll_ , met);
  } catch(...){}
  
  
  if(genmet.isValid()){
    if(METDebug) cout << "Collection " << metcoll_ << "is of type genMet" << endl;
    this->FillAllMET(*genmet , METVector);
  }
  else if(pfmet.isValid()){
    if(METDebug) cout << "Collection " << metcoll_ << "is of type pfMet" << endl;
    this->FillAllMET(*pfmet , METVector);
  }
  else if(calomet.isValid()){
    if(METDebug) cout << "Collection " << metcoll_ << "is of type caloMet" << endl;
    this->FillAllMET(*calomet , METVector);
  }
  else if(met.isValid()){
    if(METDebug) cout << "Collection " << metcoll_ << "is of type Met" << endl;
    this->FillAllMET(*met , METVector);
  }
  else  
    cout << "[UABaseTree::GetMET] Met collection " << metcoll_ << " is either not found, or not of a known type." << endl;
  
}


void UABaseTree::GetAllMETs( const edm::Event& iEvent ){
  for(vector<InputTag>::iterator icoll = mets_.begin() ; icoll!= mets_.end() ; ++icoll)
    this->GetMET(iEvent , icoll->label() , allMETs[icoll->label()] );
}

template <class T>
void UABaseTree::FillAllMET(const vector<T>& in , vector<MyMet>& out){
  out.clear();
  out.assign(in.size() , MyMet());

  for(unsigned int i = 0 ; i < in.size() ; ++i){
    //out[i].SetPxPyPzE( in[i].px() , in[i].py() , in[i].pz() , in[i].energy() );
    if(METDebug) cout << in[i].px() << "  " << in[i].py() << "  " << in[i].pz() << "  " << in[i].energy() << endl;
    out[i].SetPxPyPzE( in[i].px() , in[i].py() , in[i].pz() , sqrt( pow(in[i].pt(),2) + pow(in[i].pz(),2)) );
    out[i].sumet   = in[i].sumEt();
    out[i].elongit = in[i].e_longitudinal();

    if(METDebug) out[i].Print();
  }

}
