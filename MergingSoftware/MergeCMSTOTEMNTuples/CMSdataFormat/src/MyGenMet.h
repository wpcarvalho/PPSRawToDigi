#ifndef __MyGenMet_H__
#define __MyGenMet_H__

#include "TObject.h"
#include "TLorentzVector.h"

class MyGenMet : public TObject {

  public :
    MyGenMet();
  ~MyGenMet();

  void Reset();
  void Print();

  //GenMet from MET Collection
  Double_t Met ;
  Double_t MetX;
  Double_t MetY;
  Double_t MetPhi;

  //GenMet from GenPart Collection (status 1 and 3)
  TLorentzVector MetGP1 ;
  TLorentzVector MetGP3 ;
 
  private:

  ClassDef (MyGenMet,1)
};

#endif
