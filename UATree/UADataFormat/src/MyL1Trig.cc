
#include "MyL1Trig.h"
#include <iostream>

using namespace std;

ClassImp(MyL1Trig)

MyL1Trig::MyL1Trig(){
  this->Reset();
}

MyL1Trig::~MyL1Trig(){}

void MyL1Trig::Reset() {

  fPhysMap.clear();

  for (int i = 0; i < 128; i++) {
    fPhysMask[i] = false;
    fPhysDecisionBefore[i] = false;
    fPhysDecisionAfter[i] = false;
    fPhysPrescale[i] = 0;
    fPhysAlias[i] = "";
  }
  
  fTechMap.clear();

  for (int i = 0; i < 64; i++) {
    fTechMask[i] = false;
    fTechDecisionBefore[i] = false;
    fTechDecisionAfter[i] = false;
    fTechPrescale[i] = 0;
    fTechAlias[i] = "";
  }
}

void MyL1Trig::Print() {
  for (TrigMap::const_iterator it = fPhysMap.begin(); it != fPhysMap.end(); ++it) 
    cout<<"L1 Physical Trigger with name "<<(*it).first<<" has bit number "<<(*it).second<<endl;
  
  for (int i=0 ; i <128 ; i++) 
    cout<<"L1 Physical Trigger bit: "<<i<<" mask: "<<fPhysMask[i]<<" Decision Before Mask: "<<fPhysDecisionBefore[i]
	<<" Decision After Mask: "<<fPhysDecisionAfter[i]<<" Prescale: "<<fPhysPrescale[i]<<" Alias: "<<fPhysAlias[i]<<endl;

  for (TrigMap::const_iterator it = fTechMap.begin(); it != fTechMap.end(); ++it)
    cout<<"L1 Technical Trigger with name "<<(*it).first<<" has bit number "<<(*it).second<<endl;

  for (int i=0 ; i <64 ; i++)
    cout<<"L1 Technical Trigger bit: "<<i<<" mask: "<<fTechMask[i]<<" Decision Before Mask: "<<fTechDecisionBefore[i]
        <<" Decision After Mask: "<<fTechDecisionAfter[i]<<" Prescale: "<<fTechPrescale[i]<<" Alias: "<<fTechAlias[i]<<endl;
}

//-- L1 Physical Trigger (Algo Trigger)

void MyL1Trig::SetPhys(string name,unsigned bit,bool mask,bool decisionBeforeMask,bool decisionAfterMask,unsigned prescale,string alias) {
  fPhysMap[name] = bit;
  fPhysMask[bit] = mask;
  fPhysDecisionBefore[bit] = decisionBeforeMask;
  fPhysDecisionAfter[bit] = decisionAfterMask;
  fPhysPrescale[bit] = prescale;
  fPhysAlias[bit] = alias;
}


string MyL1Trig::GetPhysName(unsigned bit) {

  string name = "";

  for (TrigMap::const_iterator it = fPhysMap.begin(); it != fPhysMap.end(); ++it)
    if (it->second == bit) name = it->first;

  if(name == "") cout<<"Trigger bit not found, put name to empty"<<endl;

  return(name);
}


bool MyL1Trig::GetPhysMask(unsigned bit) { return(fPhysMask[(int) bit]); }

bool MyL1Trig::GetPhysDecisionBefore(unsigned bit) { return(fPhysDecisionBefore[(int) bit]); }

bool MyL1Trig::GetPhysDecisionAfter(unsigned bit) { return(fPhysDecisionAfter[(int) bit]); }

unsigned MyL1Trig::GetPhysPrescale(unsigned bit) { return(fPhysPrescale[(int) bit]); }

string MyL1Trig::GetPhysAlias(unsigned bit) { return(fPhysAlias[(int) bit]); }


unsigned MyL1Trig::GetPhysBit(const std::string & name) {

  TrigMap::const_iterator it = fPhysMap.find(name);
  unsigned bit = bit_max;

  if (it != fPhysMap.end()) bit = it->second;
  else cout<<"Trigger name not found, put bit to "<<bit_max<<endl;
    
  return(bit);
}

bool MyL1Trig::GetPhysMaskByName(const std::string & name) {

  bool mask = false;
  unsigned bit = GetPhysBit(name);
  
  if(bit != bit_max) mask = GetPhysMask(bit);
  else cout<<"Trigger name not found, put mask to false"<<endl;

  return(mask);
}
  
bool MyL1Trig::GetPhysDecisionBeforeByName(const std::string & name) {

  bool decision = false;
  unsigned bit = GetPhysBit(name);
  
  if (bit != bit_max) decision = GetPhysDecisionBefore(bit);
  else cout<<"Trigger name not found, put decision before to false"<<endl;

  return(decision);
}

bool MyL1Trig::GetPhysDecisionAfterByName(const std::string & name) {

  bool decision = false;
  unsigned bit = GetPhysBit(name);

  if (bit != bit_max) decision = GetPhysDecisionAfter(bit);
  else cout<<"Trigger name not found, put decision after to false"<<endl;

  return(decision);
}

unsigned MyL1Trig::GetPhysPrescaleByName(const std::string & name) {

  unsigned prescale = 0;
  unsigned bit = GetPhysBit(name);

  if (bit != bit_max) prescale = GetPhysPrescale(bit);
  else cout<<"Trigger name not found, put prescale to zero"<<endl;

  return(prescale);
}

string MyL1Trig::GetPhysAliasByName(const std::string & name) {

  string alias = "";
  unsigned bit = GetPhysBit(name);

  if (bit != bit_max) alias = GetPhysAlias(bit);
  else cout<<"Trigger name not found, put alias to empty"<<endl;

  return(alias);
}

//-- L1 Technical Trigger 

void MyL1Trig::SetTech(string name,unsigned bit,bool mask,bool decisionBeforeMask,bool decisionAfterMask,unsigned prescale,string alias) {
  fTechMap[name] = bit;
  fTechMask[bit] = mask;
  fTechDecisionBefore[bit] = decisionBeforeMask;
  fTechDecisionAfter[bit] = decisionAfterMask;
  fTechPrescale[bit] = prescale;
  fTechAlias[bit] = alias;
}


string MyL1Trig::GetTechName(unsigned bit) {

  string name = "";

  for (TrigMap::const_iterator it = fTechMap.begin(); it != fTechMap.end(); ++it)
    if (it->second == bit) name = it->first;

  if(name == "") cout<<"Trigger bit not found, put name to empty"<<endl;

  return(name);
}


bool MyL1Trig::GetTechMask(unsigned bit) { return(fTechMask[(int) bit]); }

bool MyL1Trig::GetTechDecisionBefore(unsigned bit) { return(fTechDecisionBefore[(int) bit]); }

bool MyL1Trig::GetTechDecisionAfter(unsigned bit) { return(fTechDecisionAfter[(int) bit]); }

unsigned MyL1Trig::GetTechPrescale(unsigned bit) { return(fTechPrescale[(int) bit]); }

string MyL1Trig::GetTechAlias(unsigned bit) { return(fTechAlias[(int) bit]); }


unsigned MyL1Trig::GetTechBit(const std::string & name) {

  TrigMap::const_iterator it = fTechMap.find(name);
  unsigned bit = bit_max;

  if (it != fTechMap.end()) bit = it->second;
  else cout<<"Trigger name not found, put bit to "<<bit_max<<endl;
    
  return(bit);
}

bool MyL1Trig::GetTechMaskByName(const std::string & name) {

  bool mask = false;
  unsigned bit = GetTechBit(name);
  
  if(bit != bit_max) mask = GetTechMask(bit);
  else cout<<"Trigger name not found, put mask to false"<<endl;

  return(mask);
}
  
bool MyL1Trig::GetTechDecisionBeforeByName(const std::string & name) {

  bool decision = false;
  unsigned bit = GetTechBit(name);
  
  if (bit != bit_max) decision = GetTechDecisionBefore(bit);
  else cout<<"Trigger name not found, put decision before to false"<<endl;

  return(decision);
}

bool MyL1Trig::GetTechDecisionAfterByName(const std::string & name) {

  bool decision = false;
  unsigned bit = GetTechBit(name);

  if (bit != bit_max) decision = GetTechDecisionAfter(bit);
  else cout<<"Trigger name not found, put decision after to false"<<endl;

  return(decision);
}

unsigned MyL1Trig::GetTechPrescaleByName(const std::string & name) {

  unsigned prescale = 0;
  unsigned bit = GetTechBit(name);

  if (bit != bit_max) prescale = GetTechPrescale(bit);
  else cout<<"Trigger name not found, put prescale to zero"<<endl;

  return(prescale);
}

string MyL1Trig::GetTechAliasByName(const std::string & name) {

  string alias = "";
  unsigned bit = GetTechBit(name);

  if (bit != bit_max) alias = GetTechAlias(bit);
  else cout<<"Trigger name not found, put alias to empty"<<endl;

  return(alias);
}

