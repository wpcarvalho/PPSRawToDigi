#ifndef __MyL1Trig_H__
#define __MyL1Trig_H__

#include "TObject.h"
#include <map>
using namespace std;

class MyL1Trig : public TObject {
  
  typedef map<string,unsigned> TrigMap;  //-- string name, unsigned bit

 private:

  const static unsigned bit_max = 500;

  //-- L1 Physical Trigger (Algo Trigger)
  
  TrigMap fPhysMap;             
  
  bool fPhysMask[128];
  bool fPhysDecisionBefore[128]; //-- L1 algorithm decision, not considering the mask
  bool fPhysDecisionAfter[128];  //-- L1 algorithm decision, considering the mask
  unsigned fPhysPrescale[128];
  string fPhysAlias[128];
  
  //-- L1 Technical Trigger 

  TrigMap fTechMap; 
  
  bool fTechMask[64];
  bool fTechDecisionBefore[64]; //-- L1 algorithm decision, not considering the mask
  bool fTechDecisionAfter[64];  //-- L1 algorithm decision, considering the mask
  unsigned fTechPrescale[64];
  string fTechAlias[64];
  
 public :
 
    MyL1Trig();
  virtual ~MyL1Trig();
  
  void Reset();
  void Print();  

  //-- L1 Physical Trigger (Algo Trigger)

  void SetPhys(string name,unsigned bit,bool mask,bool decisionBeforeMask,bool decisionAfterMask,unsigned prescale,string alias); 
  
  string GetPhysName(unsigned bit);
  bool GetPhysMask(unsigned bit);
  bool GetPhysDecisionBefore(unsigned bit);
  bool GetPhysDecisionAfter(unsigned bit);
  unsigned GetPhysPrescale(unsigned bit);
  string GetPhysAlias(unsigned bit);

  unsigned GetPhysBit(const std::string & name);
  bool GetPhysMaskByName(const std::string & name);
  bool GetPhysDecisionBeforeByName(const std::string & name);
  bool GetPhysDecisionAfterByName(const std::string & name); 
  unsigned GetPhysPrescaleByName(const std::string & name);
  string GetPhysAliasByName(const std::string & name);

  //-- L1 Technical Trigger

  void SetTech(string name,unsigned bit,bool mask,bool decisionBeforeMask,bool decisionAfterMask,unsigned prescale,string alias);

  string GetTechName(unsigned bit);
  bool GetTechMask(unsigned bit);
  bool GetTechDecisionBefore(unsigned bit);
  bool GetTechDecisionAfter(unsigned bit);
  unsigned GetTechPrescale(unsigned bit);
  string GetTechAlias(unsigned bit);

  unsigned GetTechBit(const std::string & name);
  bool GetTechMaskByName(const std::string & name);
  bool GetTechDecisionBeforeByName(const std::string & name);
  bool GetTechDecisionAfterByName(const std::string & name);
  unsigned GetTechPrescaleByName(const std::string & name);
  string GetTechAliasByName(const std::string & name);

  ClassDef (MyL1Trig,1)
};

#endif
