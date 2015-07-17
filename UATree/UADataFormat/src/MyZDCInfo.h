#ifndef __MyZDCInfo_H__
#define __MyZDCInfo_H__

#include "TObject.h"

#include <map>

class MyZDCInfo : public TObject {
  
 public :
   MyZDCInfo();
  ~MyZDCInfo();
  
  void Reset();
  void Print();
 
  /*Int_t    nHits;
  Double_t sumEnergyEMPlus;
  Double_t sumEnergyEMMinus;
  Double_t sumEnergyHADPlus;
  Double_t sumEnergyHADMinus;*/
  std::map<int,int>    nHitsPerChannel;
  std::map<int,double> sumEnergyPerChannel;

 private:
  ClassDef (MyZDCInfo,1)
};

#endif
