#ifndef __MyPFCand_H__
#define __MyPFCand_H__

#include "MyPart.h"

class MyPFCand : public MyPart {

  public :
    MyPFCand();
    virtual ~MyPFCand();
    
    virtual void Reset();
    virtual void Print();
    
    enum ParticleType {
      X=0,         // undefined
      h,           // charged hadron
      e,           // electron 
      mu,          // muon 
      gamma,       // photon
      h0,          // neutral hadron
      h_HF,        // HF tower identified as a hadron
      egamma_HF    // HF tower identified as an EM PFCandicle
    };

   ParticleType  particleId;
    
    
  private:

    ClassDef (MyPFCand,1)
};



#endif


