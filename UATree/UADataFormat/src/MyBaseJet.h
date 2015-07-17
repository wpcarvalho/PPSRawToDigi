#ifndef __MyBaseJet_H__
#define __MyBaseJet_H__

#include "MyPart.h"

class MyBaseJet : public MyPart {

  public :
    MyBaseJet();
    ~MyBaseJet();

    virtual void Print();
    virtual void Reset();

    //-- correction - uncertainty
    Double_t jec, jec_unc;  
  
  private :
  
  ClassDef (MyBaseJet,1)
};

#endif


