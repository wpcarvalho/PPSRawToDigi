#include "MyEvtId.h"
#include <iostream>

using namespace std;

ClassImp(MyEvtId)

MyEvtId::MyEvtId(){
  this->Reset();
}

MyEvtId::~MyEvtId(){}

void MyEvtId::Reset() {
  Run       = 0;
  Evt       = 0;
  LumiSect  = 0;
  Time      = -1;
  IsData    = false;
  ExpType   = 0;
  Bunch     = -1;
  Orbit     = -1;
}

void MyEvtId::Print() {
  cout<<"Run: "<<this->Run<<endl;
  cout<<"Evt: "<<this->Evt<<endl;
  cout<<"LumiSect: "<<this->LumiSect<<endl;

  cout<<"Time: "<<this->Time<<endl;
  cout<<"Is Data: "<<this->IsData<<endl;
  cout<<"Exp Type: "<<this->ExpType<<endl;
  cout<<"Bx: "<<this->Bunch<<endl;
  cout<<"Orbit: "<<this->Orbit<<endl;
}
