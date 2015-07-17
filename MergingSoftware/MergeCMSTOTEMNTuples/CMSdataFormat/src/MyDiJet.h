#ifndef __MyDiJet_H__
#define __MyDiJet_H__

#include <vector>
#include "TObject.h"

class MyDiJet : public TObject {
  
 public :
    MyDiJet();
  ~MyDiJet();
  
  void Reset();
  void Print();
  
  Bool_t isDiJet;
  Short_t posJet1, posJet2;
  
 private:
  
  ClassDef (MyDiJet,1)
};

#endif
