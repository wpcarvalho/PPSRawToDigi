#ifndef __MyFwdGap_H__
#define __MyFwdGap_H__

#include "TObject.h"

class MyFwdGap : public TObject {

  public :
    MyFwdGap();
  virtual ~MyFwdGap();
  
  void Reset();
  void Print();

  Int_t nTowersHF_plus ;
  Int_t nTowersHF_minus ;
  Int_t nTowersHE_plus ;
  Int_t nTowersHE_minus ;
  Int_t nTowersHB_plus ;
  Int_t nTowersHB_minus ;
  Int_t nTowersEE_plus ;
  Int_t nTowersEE_minus ;
  Int_t nTowersEB_plus ;
  Int_t nTowersEB_minus ; 
  //Sum(E)
  Double_t sumEHF_plus ;
  Double_t sumEHF_minus ;
  Double_t sumEHE_plus ;
  Double_t sumEHE_minus ;
  Double_t sumEHB_plus ;
  Double_t sumEHB_minus ;
  Double_t sumEEE_plus ;
  Double_t sumEEE_minus ;
  Double_t sumEEB_plus ;
  Double_t sumEEB_minus ;
  // Sum(ET)
  Double_t sumETHF_plus ;
  Double_t sumETHF_minus ;
  Double_t sumETHE_plus ;
  Double_t sumETHE_minus ;
  Double_t sumETHB_plus ;
  Double_t sumETHB_minus ;
  Double_t sumETEE_plus ;
  Double_t sumETEE_minus ;
  Double_t sumETEB_plus ;
  Double_t sumETEB_minus ;

  private:

  ClassDef (MyFwdGap,1)
};

#endif
