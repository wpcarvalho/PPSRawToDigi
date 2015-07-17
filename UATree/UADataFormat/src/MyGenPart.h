#ifndef __MyGenPart_H__
#define __MyGenPart_H__

#include <string>
#include "MyPart.h"

class MyGenPart : public MyPart {

  public :
    MyGenPart();
    virtual ~MyGenPart();
    
    virtual void Reset();
    virtual void Print();

    Int_t  pdgId,status,mo1,mo2,da1,da2;
    //std::string name;

  private:

  ClassDef (MyGenPart,1)
};

#endif

