#ifndef TotemHisto_H
#define TotemHisto_H


#include "SimG4CMS/Forward/interface/TotemHit.h"


// ROOT staff 
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom.h"

#define NUMVARROOT 19
#include <vector>





class TotemHisto : public TObject {

public: 

  TotemHisto();
  TotemHisto(std::string);
  virtual ~TotemHisto();

  void set_EVT(int v){evt=v;}
  void set_X(float v) {x = v;}
  void set_Y(float v) {y = v;}
  void set_Z(float v) {z = v;}
  void set_UID(int v) {UID = v;}
  void set_Ptype(int v) {Ptype = v;}
  void set_TID(int v) {TID = v;}
  void set_PID(int v) {PID = v;}
  void set_ELoss(float v) {ELoss = v;}
  void set_PABS(float v) {PABS = v;}
  void set_VX(float v) {vx = v;}
  void set_VY(float v) {vy = v;}
  void set_VZ(float v) {vz = v;}
  void set_PX(float v){Px = v;}
  void set_PY(float v){Py = v;}
  void set_PZ(float v){Pz = v;}
  void set_VPX(float v){VPx = v;}
  void set_VPY(float v){VPy = v;}
  void set_VPZ(float v){VPz = v;}

  void fillNtuple();


private: 
  TFile *rt_hf;

// TotemHisto();

 TNtuple *ntuple;

 float rootvec[NUMVARROOT];

 float x,y,z,ELoss,PABS,vx,vy,vz,Px,Py,Pz,VPx,VPy,VPz;
 int UID,Ptype,TID,PID,evt;

 std::string nome_file;

};
#endif
