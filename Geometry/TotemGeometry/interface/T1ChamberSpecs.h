#ifndef _T1_CHAMBER_SPEC
#define _T1_CHAMBER_SPEC

//#define _DEBUG_
#include <iostream>
#include <string>
#include "DataFormats/T1DetId/interface/T1DetId.h"

class T1ChamberSpecs {

 public:
  T1ChamberSpecs(){}
  ~T1ChamberSpecs(){}

  // parametri globali delle camere
  double anodeCathodeSpacing() const {return 5.;}
  double wireSpacing() const {return 3.;}
  double wireRadius() const {return 0.015;}
  //controllare il pitch delle strip
  double Pitch() const {return 5.0;}
  const int nNodes() const {return 5;}
  //verificare il guadagno del gas
  //double gasGain() {return 10.e4;}
  double gasGain() const {return 1.e04;} //// !!!!!! CHECK
};


#endif


