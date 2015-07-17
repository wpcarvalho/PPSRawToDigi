#ifndef __MyL1TrigOld_H__
#define __MyL1TrigOld_H__

#include "TObject.h"


class MyL1TrigOld : public TObject {

  public :
    MyL1TrigOld();
  virtual ~MyL1TrigOld();

  bool PhysTrigWord[128];
  bool TechTrigWord[64];

  void Reset();
  void Print(); 

  private:

  ClassDef (MyL1TrigOld,1)
};

#endif
