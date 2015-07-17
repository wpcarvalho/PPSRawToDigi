#ifndef __MyCaloTower_H__
#define __MyCaloTower_H__

#include <string>
#include "MyPart.h"

using namespace std;

class MyCaloTower : public MyPart {

  public :
    MyCaloTower();
    ~MyCaloTower();

    virtual void Reset();
    virtual void Print();
    
    Double_t emEnergy;
    Double_t hadEnergy;

    Bool_t hasEB;
    Bool_t hasEE;

    Bool_t hasHB;
    Bool_t hasHE;
    Bool_t hasHF;

    Int_t zside;
    
  private:

  ClassDef (MyCaloTower,1)
};

#endif
