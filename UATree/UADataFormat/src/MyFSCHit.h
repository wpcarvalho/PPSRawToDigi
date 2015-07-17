#ifndef __MyFSCHit_H__
#define __MyFSCHit_H__

#include "TObject.h"

class MyFSCHit : public TObject {
  
 public :
   MyFSCHit();
  ~MyFSCHit();
  
  void Reset();
  void Print();
 
  Int_t side; 
  Int_t section;
  Int_t channel;
  Int_t channelId;
  Double_t energy;
  Double_t time;

 private:
  ClassDef (MyFSCHit,1)
};

#endif
