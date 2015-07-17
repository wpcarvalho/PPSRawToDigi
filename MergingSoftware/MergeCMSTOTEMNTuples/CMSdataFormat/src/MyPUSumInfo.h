#ifndef __MyPUSumInfo_H__
#define __MyPUSumInfo_H__

#include <vector>
#include "TObject.h"

using namespace std;

class MyPUSumInfo : public TObject {

  public :
    MyPUSumInfo();
    virtual ~MyPUSumInfo();

    void Reset();
    void Print();

    Int_t             nPU          ;
    vector<Float_t>   zposition    ;
    vector<Float_t>   sumpT_lowpT  ;
    vector<Float_t>   sumpT_highpT ;
    vector<Int_t>     ntrks_lowpT  ;
    vector<Int_t>     ntrks_highpT ;

  private:

    ClassDef (MyPUSumInfo,1)
};

#endif

