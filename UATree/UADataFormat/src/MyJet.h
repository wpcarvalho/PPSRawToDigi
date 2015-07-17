#ifndef __MyJet_H__
#define __MyJet_H__

#include <map>
#include <string>
#include "TObject.h"
#include "MyBaseJet.h"

using namespace std;

typedef map<string,MyBaseJet>::iterator map_it;

class MyJet : public TObject {

  public :
    MyJet();
    ~MyJet();

    virtual void Reset();
    virtual void Print();
    
    //map containing raw + wanted corrections
    map<string,MyBaseJet> mapjet;
    
    //number of constituents (PFObject for PFjet, CaloTower for Calojet)
    UInt_t nconstituent;

    //-- jet ID
    Bool_t LooseJetId;
    Bool_t TightJetId;

  private:

  ClassDef (MyJet,1)
};

#endif
