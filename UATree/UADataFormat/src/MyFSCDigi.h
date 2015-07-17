#ifndef __MyFSCDigi_H__
#define __MyFSCDigi_H__

#include "TObject.h"

#include <vector>

class MyFSCDigi : public TObject {
  
 public :
   MyFSCDigi();
  ~MyFSCDigi();
  
  void Reset();
  void Print();
 
  Int_t side; 
  Int_t section;
  Int_t channel;
  Int_t channelId;
  std::vector<int> digiADC;
  std::vector<float> digifC;

 private:
  ClassDef (MyFSCDigi,1)
};

#endif
