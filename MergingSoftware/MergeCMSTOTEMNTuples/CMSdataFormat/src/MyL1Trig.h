#ifndef __MyL1Trig_H__
#define __MyL1Trig_H__

#include "TObject.h"


class MyL1Trig : public TObject {

  public :
    MyL1Trig();
  virtual ~MyL1Trig();

  bool PhysTrigWord[128];
  bool TechTrigWord[64];

  void Reset();
  void Print(); 

  private:

  ClassDef (MyL1Trig,1)
};

#endif
