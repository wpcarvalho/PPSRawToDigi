#ifndef __MyMet_H__
#define __MyMet_H__

#include "MyPart.h"

class MyMet : public MyPart {

  public :
    MyMet();
    virtual ~MyMet();
    
    virtual void Reset();
    virtual void Print();

    double sumet;
    double elongit;

  private:

    ClassDef (MyMet,1)
};



#endif


