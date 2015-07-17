#ifndef __MyCastorJet_H__
#define __MyCastorJet_H__

#include "TObject.h"

using namespace std;

class MyCastorJet : public TObject {
  
 public :
    MyCastorJet();
  virtual ~MyCastorJet();
  
  virtual void Reset(){};
  virtual void Print();
  
  Double_t energy, pt, eta, phi;
  Double_t fem, eem, ehad;
  Double_t width, depth, fhot, sigmaz;
  Int_t ntower;
  
 private:
  
  ClassDef (MyCastorJet,1)
};

#endif

