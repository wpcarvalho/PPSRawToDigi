#ifndef __MyTrackJet_H__
#define __MyTrackJet_H__

#include "MyBaseJet.h"
//#include "MyJet.h"
#include "MyTracks.h"
#include <vector>
#include <map>

class MyTrackJet : public MyBaseJet {

  public :
    MyTrackJet();
    ~MyTrackJet();

    virtual void Print();
    virtual void Reset();

    //RefById to associated vertex
    Int_t vtxId;

    //map containing raw + wanted corrections
    map<string,MyBaseJet> mapjet;

    //-- raw jet variables
    Double_t e_raw,pt_raw,eta_raw,phi_raw,px_raw,py_raw,pz_raw;
    //-- corrected jet variables
    Double_t e_cal,pt_cal,eta_cal,phi_cal,px_cal,py_cal,pz_cal;

    //-- number of tracks  
    Int_t ntrack;
    vector<MyTracks> vtracks;

    //-- check jet to be associated to the hard primary vertex
    Bool_t trackjet_pv;

  private :
  
  ClassDef (MyTrackJet,1)
};

#endif


