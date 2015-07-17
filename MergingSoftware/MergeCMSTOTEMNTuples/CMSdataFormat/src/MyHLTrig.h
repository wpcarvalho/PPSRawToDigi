#ifndef __MyHLTrig_H__
#define __MyHLTrig_H__

#include "TObject.h"
#include <map>
using namespace std;

class MyHLTrig : public TObject {

  public :
  
    MyHLTrig();
    ~MyHLTrig();
  
    void Reset();
    void Print();

    map<string,bool>   HLTmap;
    map<string,double> HLTprescale;

  private:

    ClassDef (MyHLTrig,1)
};

#endif



