#ifndef _Totem_T1Event_include_h_
#define _Totem_T1Event_include_h_

#include <TObject.h>
#include <vector>

class T1Event : public TObject {
 public:
  T1Event();

  std::vector<int> DigiStripA;
  std::vector<int> DigiStripB;
  std::vector<int> DigiWire;
  std::vector<int> DigiStripA_Arm;
  std::vector<int> DigiStripA_Plane;
  std::vector<int> DigiStripA_CSC;
  std::vector<int> DigiStripB_Arm;
  std::vector<int> DigiStripB_Plane;
  std::vector<int> DigiStripB_CSC;
  std::vector<int> DigiWire_Arm;
  std::vector<int> DigiWire_Plane;
  std::vector<int> DigiWire_CSC;
  std::vector<double> TrkAx;
  std::vector<double> TrkAy;
  std::vector<double> TrkX0;
  std::vector<double> TrkY0;
  std::vector<double> TrkPhi;
  std::vector<double> TrkEta;
  std::vector<double> TrkZatRmin;
  std::vector<double> TrkRmin;
  std::vector<double> TrkChi2OverN;
  std::vector<double> RecoHitX;
  std::vector<double> RecoHitY;
  std::vector<double> RecoHitZ;
  long run_no;
  long ev_no;
  unsigned long long timestamp;
  double its_var1;
  int its_var2;

  ClassDef(T1Event, 10023);
};

#endif
