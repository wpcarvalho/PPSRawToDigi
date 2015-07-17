#ifndef __MyEvtId_H__
#define __MyEvtId_H__

#include "TObject.h"

typedef unsigned long long TimeValue_t;


class MyEvtId : public TObject {

  public :
    MyEvtId();
  virtual ~MyEvtId();

  void Reset();
  void Print();
 
  UInt_t       Run ;
  UInt_t       Evt ;
  UInt_t       LumiSect ;
  TimeValue_t  Time;
  Bool_t       IsData;
  UInt_t       ExpType;
  Int_t        Bunch;
  Int_t        Orbit;   

  private:

  ClassDef (MyEvtId,1)
};

#endif
