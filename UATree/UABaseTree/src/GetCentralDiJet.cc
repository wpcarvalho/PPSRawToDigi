//-- Description: Function to retrieve central PF Di Jet information

//--DataFormats
#include "DataFormats/Math/interface/deltaPhi.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"

bool CentralDiJetDebug = false;

void UABaseTree::GetCentralDiJet(const vector<MyJet*>& JetVector, const string correction,  MyDiJet& dijet) {

      
  //-- collection to save information

  dijet.Reset();

  //-- central jet selection

  bool accept = false;
  
  short int posJet1 = -1;
  short int posJet2 = -1;

  double refPt1 = 0;
  double refPt2 = 0;

  //-- find the two highest pt jets (corrected pt)

  for(vector<MyJet*>::const_iterator jet = JetVector.begin(); jet < JetVector.end(); ++jet) {

    Double_t pt = (*jet)->mapjet[correction].Pt();

    if(pt > refPt1) {
      refPt2 = refPt1;
      posJet2 = posJet1;      
      refPt1 = pt;
      posJet1 = jet - JetVector.begin();
    }

    else if(pt > refPt2) {
      refPt2 = pt ;
      posJet2 = jet - JetVector.begin();
    }

  } 

  //-- apply the thight selection to them

  if(posJet1 >= 0 && posJet2 >= 0) {
    
    MyJet* jet1 = JetVector.at(posJet1);
    MyJet* jet2 = JetVector.at(posJet2);

    bool accept_jet1 = false;
    bool accept_jet2 = false;


    //-- jet 1 selection
    if(jet1->mapjet[correction].Pt() > jetPtCut_ && fabs(jet1->mapjet[correction].Eta()) < jetEtaCut_ && jet1->LooseJetId == true) accept_jet1 = true;
  
    //-- jet 2 selection
    if(jet2->mapjet[correction].Pt() > jetPtCut_ && fabs(jet2->mapjet[correction].Eta()) < jetEtaCut_ && jet2->LooseJetId == true) accept_jet2 = true;
  
    //-- final selection (back-to-back)

    if(accept_jet1 == true && accept_jet2 == true) {
      double deltaPhi = fabs( reco::deltaPhi( jet1->mapjet[correction].Phi(), jet2->mapjet[correction].Phi() ) );
      if (fabs(deltaPhi - M_PI) < 1.0) accept = true;
    }

    if(accept == true) {
      dijet.isDiJet = true;
      dijet.posJet1 = posJet1;
      dijet.posJet2 = posJet2;
    }
  
  }

  if(CentralDiJetDebug) dijet.Print();
}

