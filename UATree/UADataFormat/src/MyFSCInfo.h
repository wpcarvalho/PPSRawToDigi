#ifndef __MyFSCInfo_H__
#define __MyFSCInfo_H__

#include "TObject.h"

#include <map>

class MyFSCInfo : public TObject {
  
 public :
   MyFSCInfo();
  ~MyFSCInfo();
  
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
  ClassDef (MyFSCInfo,1)
};

#endif
