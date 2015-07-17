#ifndef __MyZDCHit_H__
#define __MyZDCHit_H__

#include "TObject.h"

class MyZDCHit : public TObject {
  
 public :
   MyZDCHit();
  ~MyZDCHit();
  
  void Reset();
  void Print();
 
  Int_t side; 
  Int_t section;
  Int_t channel;
  Int_t channelId;
  Double_t energy;
  Double_t time;

 private:
  ClassDef (MyZDCHit,1)
};

#endif
