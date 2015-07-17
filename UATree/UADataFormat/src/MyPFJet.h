#ifndef __MyPFJet_H__
#define __MyPFJet_H__

#include "MyJet.h"
#include "MyTracks.h"
#include <vector>

class MyPFJet : public MyJet {

  public :
    MyPFJet();
    ~MyPFJet();

    virtual void Reset();
    virtual void Print();
  
    //-- estimators (PFjet only)
    Double_t fhad_ch;    //-- chargedHadronEnergyFraction()
    Double_t fhad_ne;    //-- neutralHadronEnergyFraction()

    Double_t fem_ch;     //-- chargedEmEnergyFraction()
    Double_t fem_ne;     //-- neutralEmEnergyFraction()

    Int_t multi_ch;      //-- chargedMultiplicity()
    Int_t multi_ne;      //-- neutralMultiplicity()

    Int_t multi_ch_had;  //-- chargedHadronMultiplicity()
    Int_t multi_ne_had;  //-- neutralHadronMultiplicity()

    Int_t multi_gamma;  //-- photonMultiplicity()
    Int_t multi_ele;    //-- electronMultiplicity()
    Int_t multi_mu;     //-- muonMultiplicity()
    
    //-- number of tracks
    Int_t ntrack;
    vector<MyTracks> vtracks;

  private:

  ClassDef (MyPFJet,1)
};

#endif
