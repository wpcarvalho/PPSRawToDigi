#ifndef __MyBeamSpot_H__
#define __MyBeamSpot_H__

#include "TObject.h"

class MyBeamSpot : public TObject {

  public :
    MyBeamSpot();
    virtual ~MyBeamSpot();

    void Reset(){};
    void Print();

    Double_t x,y,z;
    Double_t ex,ey,ez;
    Double_t sigmaZ , dxdz , dydz ;
    Double_t esigmaZ , edxdz , edydz ;
    Double_t BeamWidthX , BeamWidthY ;
    Double_t eBeamWidthX , eBeamWidthY ;


  private:

  ClassDef (MyBeamSpot,1)
};

#endif

