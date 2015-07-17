#ifndef __MyGenJet_H__
#define __MyGenJet_H__

#include "MyPart.h"
#include "MyGenPart.h"
#include <vector>


class MyGenJet : public MyPart {

  public :
  MyGenJet();
  ~MyGenJet();
  
  void Reset();
  void Print();

  Int_t                  npart;
  std::vector<MyGenPart> JetPart;

  private:

  ClassDef (MyGenJet,1)
};

#endif
