#ifndef __MySimVertex_H__
#define __MySimVertex_H__

#include "TObject.h"

class MySimVertex : public TObject {

  public :
    MySimVertex();
    virtual ~MySimVertex();

    void Reset(){};
    void Print();

    Double_t x,y,z;

  private:

  ClassDef (MySimVertex,1)
};

#endif

