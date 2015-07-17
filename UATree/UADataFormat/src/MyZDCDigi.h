#ifndef __MyZDCDigi_H__
#define __MyZDCDigi_H__

#include "TObject.h"

#include <vector>

class MyZDCDigi : public TObject {
  
 public :
   MyZDCDigi();
  ~MyZDCDigi();
  
  void Reset();
  void Print();
 
  Int_t side; 
  Int_t section;
  Int_t channel;
  Int_t channelId;
  std::vector<int> digiADC;
  std::vector<float> digifC;

 private:
  ClassDef (MyZDCDigi,1)
};

#endif
