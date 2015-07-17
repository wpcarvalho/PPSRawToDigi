// UABaseTree Analysis class decleration
#include "UATree/UABaseTree/interface/UABaseTree.h"

#include "FWCore/Utilities/interface/Parse.h"

void UABaseTree::Init(){


  fout = new TFile(outputfilename_.c_str(), "RECREATE" ) ;
  tree = new TTree("evt","evt");

  //branch name, in case separated by #
  string branch = "";


  // --------------------   Setting all branches   --------------------

  // BeamSpots
  for(vector<InputTag>::iterator icoll = this->beamspots_.begin() ; icoll!= this->beamspots_.end() ; ++icoll){
    branch = GetBranchName(*icoll);
    this->tree->Branch( branch.c_str() , &(this->allBeamSpots[icoll->label()]) );
  }

  // General Event Info 
  if(storeEvtId_)		     tree->Branch("evtId",&evtId);
  if(storeFwdGap_)                   tree->Branch("fwdGap",&fwdGap);
  
  //MC Info
  if(storeGenKin_)                   tree->Branch("genKin",&genKin);
  if(enableGenMetFromGenPart_)       tree->Branch("genMetfromGenPartst1",&genMetfromGenPartst1);
  if(enableGenMetFromGenPart_)       tree->Branch("genMetfromGenPartst3",&genMetfromGenPartst3);
  if(storeGenPart_){                 tree->Branch("genPart",&genPart);
                                     tree->Branch("simVertex",&simVertex);
   if(saveGenPartsInDifferentColls_){
                                     tree->Branch("genElec",&genElec);
                                     tree->Branch("genMu",&genMu);
                                     tree->Branch("genNu",&genNu);
   }
  }

  if(storePUSumInfo_)                tree->Branch("PUSumInfo",&pusuminfo); 
 
  //Triggers
  if(hlt_paths_.size() > 0)	     tree->Branch("HLTrig",&HLTrig);
  if(storeL1Trig_)		     tree->Branch("L1Trig",&L1Trig);
  if(storeL1TrigOld_)		     tree->Branch("L1TrigOld",&L1TrigOld);
  
  //MITEvtSel
  if(storeMITEvtSel_)		     tree->Branch("MITEvtSel",&MITEvtSel);

  //RecoTracks
  for(vector<InputTag>::iterator icoll = this->tracks_.begin() ; icoll!= this->tracks_.end() ; ++icoll){
    branch = GetBranchName(*icoll);
    this->tree->Branch( branch.c_str() , &(this->allTracks[icoll->label()]) );
  }
    
  //RecoVertices
  for(vector<InputTag>::iterator icoll = this->vertices_.begin() ; icoll!= this->vertices_.end() ; ++icoll){
    branch = GetBranchName(*icoll);
    this->tree->Branch( branch.c_str() , &(this->allVertices[icoll->label()]) );
  }
    
     
  //RecoCaloJets
  InputTag       jetcoll_;
  string         dijetcoll_;
  vector<string> corrections_;
  for(vector<PSet>::iterator icoll = vcalojets_.begin() ; icoll != vcalojets_.end() ; ++icoll){
    jetcoll_ = icoll->getUntrackedParameter<InputTag>("jetcoll",InputTag());
    branch = GetBranchName(jetcoll_);
    if(branch.size() > 0)
      this->tree->Branch( branch.c_str() , &(this->allCaloJets[jetcoll_.label()]) );
      
    //DiJets from this Jetcoll
    dijetcoll_   = icoll->getUntrackedParameter<string>("dijetcoll","");
    corrections_ = icoll->getUntrackedParameter<vector<string> >("corrections",vector<string>());
    if(find(corrections_.begin() , corrections_.end() , dijetcoll_) != corrections_.end())
      this->tree->Branch( string(GetCollName(dijetcoll_)+"DiJet").c_str() , &(this->allDiJets[string(GetColl(dijetcoll_)+"DiJet")]) );
    
  }
    
     
  //RecoPFJets
  for(vector<PSet>::iterator icoll = vpfjets_.begin() ; icoll != vpfjets_.end() ; ++icoll){
    jetcoll_ = icoll->getUntrackedParameter<InputTag>("jetcoll",InputTag());
    branch = GetBranchName(jetcoll_);
    if(branch.size() > 0)
      this->tree->Branch( branch.c_str() , &(this->allPFJets[jetcoll_.label()]) );
     
    //DiJets from this Jetcoll
    dijetcoll_   = icoll->getUntrackedParameter<string>("dijetcoll","");
    corrections_ = icoll->getUntrackedParameter<vector<string> >("corrections",vector<string>());
    if(find(corrections_.begin() , corrections_.end() , dijetcoll_) != corrections_.end())
      this->tree->Branch( string(GetCollName(dijetcoll_)+"DiJet").c_str() , &(this->allDiJets[GetColl(dijetcoll_)+"DiJet"]) );

  }
  
    
     
  //RecoGenJets
  for(vector<InputTag>::iterator icoll = genjets_.begin() ; icoll != genjets_.end() ; ++icoll){
    branch = GetBranchName(*icoll);
    this->tree->Branch( branch.c_str() , &(this->allGenJets[icoll->label()]) );
  }
  
  //BasicJets
  for(vector<InputTag>::iterator icoll = basicjets_.begin() ; icoll != basicjets_.end() ; ++icoll){
    branch = GetBranchName(*icoll);
    this->tree->Branch( branch.c_str() , &(this->allBasicJets[icoll->label()]) );
  }
  
  //TrackJets
  for(vector<InputTag>::iterator icoll = trackjets_.begin() ; icoll != trackjets_.end() ; ++icoll){
    branch = GetBranchName(*icoll);
    this->tree->Branch( branch.c_str() , &(this->allTrackJets[icoll->label()]) );
  }


  //Castor
  if(castorrechits_.label().size() > 0)    tree->Branch("castorRecHits",&castorRecHits);
  if(castorjets_.label().size() > 0 &&
     castorjetid_.label().size() > 0)      tree->Branch("castorJets",&castorJets);
  if(castordigis_.label().size() > 0)      tree->Branch("castorDigis",&castorDigis);
  
  
  
  //RecoElectrons
  for(vector<InputTag>::iterator icoll = this->electrons_.begin() ; icoll!= this->electrons_.end() ; ++icoll){
    branch = GetBranchName(*icoll);
    this->tree->Branch( branch.c_str() , &(this->allElectrons[icoll->label()]) );
  }
  
  
  //RecoMuons
  for(vector<InputTag>::iterator icoll = this->muons_.begin() ; icoll!= this->muons_.end() ; ++icoll){
    branch = GetBranchName(*icoll);
    this->tree->Branch( branch.c_str() , &(this->allMuons[icoll->label()]) );
  }


  //PFCands
  for(vector<InputTag>::iterator icoll = this->pfcands_.begin() ; icoll!= this->pfcands_.end() ; ++icoll){
    branch = GetBranchName(*icoll);
    this->tree->Branch( branch.c_str() , &(this->allPFCands[icoll->label()]) );
  }


  //MET
  for(vector<InputTag>::iterator icoll = this->mets_.begin() ; icoll!= this->mets_.end() ; ++icoll){
    branch = GetBranchName(*icoll);
    this->tree->Branch( branch.c_str() , &(this->allMETs[icoll->label()]) );
  }
 

  //CaloObjects
  if(storeCaloObjects_) this->tree->Branch( "caloTowers" , &caloTowers );
   
  // ZDC 
  if(storeZDCInfo_  && zdcrechits_.label().size() > 0) this->tree->Branch( "zdcInfo" , &zdcInfo );
  if(storeZDCHits_  && zdcrechits_.label().size() > 0) this->tree->Branch( "zdcHits" , &zdcHits );
  if(storeZDCDigis_ && zdcdigis_.label().size() > 0)   this->tree->Branch( "zdcDigis" , &zdcDigis );

  // FSC 
  if(storeFSCInfo_  && fscrechits_.label().size() > 0) this->tree->Branch( "fscInfo" , &fscInfo );
  if(storeFSCHits_  && fscrechits_.label().size() > 0) this->tree->Branch( "fscHits" , &fscHits );
  if(storeFSCDigis_ && fscdigis_.label().size() > 0)   this->tree->Branch( "fscDigis" , &fscDigis );
}


const string UABaseTree::GetBranchName(InputTag& itag , Bool_t deleteBranchFromString){
 
  // string is delimited by #
  vector<std::string> tokens = tokenize(itag.label() , "#");
  size_t nwords = tokens.size();
  if(nwords == 0) return "";
  else if(nwords > 2) {
    cout << "[UABaseTree::getBranchName] You have too many # in the InputTag " << itag << endl;
    cout << "                            " << tokens[0] << "returned & stored in InputTag !" << endl;
    itag = InputTag(tokens[0] , itag.instance() , itag.process());
    return tokens[0];
  }
  else if(nwords == 1) return tokens[0];
  else if(nwords == 2){
    if(deleteBranchFromString) itag = InputTag(tokens[0] , itag.instance() , itag.process());
    return tokens[1];
  }
  else return 0;
}

const string UABaseTree::GetCollName(const string& str){
  
  // string is delimited by #
  vector<std::string> tokens  = tokenize(str , "#");
  size_t nwords = tokens.size();
  if(nwords == 0) return "";
  else if(nwords > 2) {
    cout << "[UABaseTree::GetCollName] You have too many # in the string " << str << endl;
    cout << "                          " << tokens[0] << "returned !" << endl;
    return tokens[0];
  }
  else if(nwords == 1) return tokens[0];
  else if(nwords == 2){
    vector<std::string> tokens_itag = tokenize(tokens[1] , ":");
    return tokens_itag[0];
  }
  else return 0;
}


const string UABaseTree::GetColl(const string& str){
  
  // string is delimited by #
  vector<std::string> tokens      = tokenize(str , "#");
  size_t nwords = tokens.size();
  if(nwords == 0) return "";
  if(nwords > 2) {
    cout << "[UABaseTree::GetCollInputTag] You have too many # in the string " << str << endl;
    cout << "                              " << tokens[0] << "returned !" << endl;
  }
  return tokens[0];
}

const InputTag UABaseTree::GetCollInputTag(const string& str){
  
  // string is delimited by #
  vector<std::string> tokens      = tokenize(str , "#");
  size_t nwords = tokens.size();
  if(nwords == 0) return InputTag();
  else if(nwords > 2) {
    cout << "[UABaseTree::GetCollInputTag] You have too many # in the string " << str << endl;
    cout << "                              " << tokens[0] << "returned !" << endl;
    return tokens[0];
  }
  else if(nwords == 1) return InputTag(tokens[0]);
  else if(nwords == 2){
    vector<std::string> tokens_itag = tokenize(tokens[1] , ":");
    
    if(tokens_itag.size() > 3){
      cout << "[UABaseTree::GetCollInputTag] You have too many \":\" in the string, wrong InputTag " << str << endl;
      cout << "                              " << tokens[0] << "returned as InputTag !" << endl;
      return InputTag(tokens[0]);
    }
   
    for(unsigned s = 1 ; s < tokens_itag.size() ; ++s)
      tokens[0] += ":" + tokens_itag[s];
    return tokens[0];
  }
  else return InputTag();
}



//------ NOT NEEDED ==> Parameters should always be input tags !!!

/*const string UABaseTree::getBranchName(string& str , Bool_t deleteBranchFromString){
  // string is delimited by #
  std::vector<std::string> tokens = tokenize(str, "#");
  size_t nwords = tokens.size();
  if(nwords > 2) {
    cout << "[UABaseTree::getBranchName] You have too many # in the InputTag " << str << endl;
  }
  if(nwords == 1) return tokens[0];
  if(nwords == 2){
  
    
    return tokens[1];
  }
}*/
