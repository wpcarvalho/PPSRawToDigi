#ifndef ROOTCLASST1Event
#define ROOTCLASST1Event
#include "TObject.h"
#include <vector>
class T1Event : public TObject
{
 public:
  T1Event();
 ~T1Event();

 
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
 std::vector<unsigned int> TrkHits;

 
 // track hits?????


 
 std::vector<double> RecoHitX;
 std::vector<double> RecoHitY;
 std::vector<double> RecoHitZ;
 





 long int run_no, ev_no;
 unsigned long long timestamp;

 
 //std::vector<double> TrkY0;
 //std::vector<double> TrkY0;


 private:
 double its_var1;
 int its_var2;
 ClassDef(T1Event,1)
};
#endif
