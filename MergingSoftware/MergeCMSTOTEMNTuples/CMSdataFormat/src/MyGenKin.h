#ifndef __MyGenKin_H__
#define __MyGenKin_H__

#include "TObject.h"
#include "TLorentzVector.h"

class MyGenKin : public TObject {

  public :
    MyGenKin();
    virtual ~MyGenKin();

  void Reset();
  void Print();

  // General Info
  Int_t    MCProcId ;   
  Double_t PtHat  ;
  Double_t genWeight  ;

  //-- Pdf Info
  Double_t x1;        //-- fraction of beam momentum carried by first parton 
  Double_t x2;        //-- fraction of beam momentum carried by second parton
  Double_t Q;         //-- Q-scale used in evaluation of PDFs (in GeV)
  Int_t    Part1Id ;  //-- flavour code of first parton
  Int_t    Part2Id ;  //-- flavour code of second parton

  // Kfactor for Signal
  Double_t kfactor;

  private:

  ClassDef (MyGenKin,1)
};

#endif
