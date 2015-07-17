#ifndef __MyCastorDigi_H__
#define __MyCastorDigi_H__

#include "TObject.h"

using namespace std;

class MyCastorDigi : public TObject {
  
 public :
    MyCastorDigi();
  virtual ~MyCastorDigi();
  
  virtual void Reset();
  virtual void Print();
  
  vector<Double_t> adc;
  vector<Double_t> fC;
  Int_t mod, sec, cha;
  
 private:
  
  ClassDef (MyCastorDigi,1)
};

#endif

