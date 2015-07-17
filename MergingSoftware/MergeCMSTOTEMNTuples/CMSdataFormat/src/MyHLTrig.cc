#include "MyHLTrig.h"
#include <iostream>
using namespace std;

ClassImp(MyHLTrig)

MyHLTrig::MyHLTrig(){
  this->Reset();
}

MyHLTrig::~MyHLTrig(){}

void MyHLTrig::Reset(){
  this->HLTmap.clear();
  this->HLTprescale.clear();
}

void MyHLTrig::Print(){
  for(map<string,bool>::iterator it = HLTmap.begin() ; it != HLTmap.end() ;  ++it){
    cout << it->first << " :    decision  " << it->second
         << "      prescale  " << HLTprescale[it->first] << endl;
  }
}
