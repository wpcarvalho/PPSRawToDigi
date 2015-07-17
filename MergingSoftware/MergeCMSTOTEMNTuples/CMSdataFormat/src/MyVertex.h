#ifndef __MyVertex_H__
#define __MyVertex_H__

#include "TObject.h"

class MyVertex : public TObject {

  public :
    MyVertex();
    virtual ~MyVertex();

    void Reset();
    void Print();

    Int_t    id;

    Double_t x,y,z;
    Double_t ex,ey,ez;    

    Bool_t   validity;  
    Bool_t   fake    ;
    Double_t chi2  ;
    Double_t ndof  ;
    
    Int_t    ntracks;
    Double_t SumPtTracks;

    Double_t chi2n();

  private:

    ClassDef (MyVertex,1)
};

#endif

