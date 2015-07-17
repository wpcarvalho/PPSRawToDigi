#ifndef __MyCastorRecHit_H__
#define __MyCastorRecHit_H__

#include "TObject.h"
#include <vector>

using namespace std;

class MyCastorRecHit : public TObject {
  
 public :
    MyCastorRecHit();
  virtual ~MyCastorRecHit();
  
  virtual void Reset(){};
  virtual void Print();
  
  Double_t energy, fC, time;
  Int_t mod, sec, cha;
  Double_t smearing, energy_smeared, fC_smeared;

 private:
  
  ClassDef (MyCastorRecHit,1)
};

#endif

