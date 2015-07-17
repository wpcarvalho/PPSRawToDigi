#include "MyGenKin.h"
#include <iostream>
using namespace std;

ClassImp(MyGenKin)

MyGenKin::MyGenKin(){
  this->Reset();
}

MyGenKin::~MyGenKin(){}

void MyGenKin::Reset() {
 MCProcId = 0;
 PtHat    = 0;
 genWeight    = 0.;
 x1       = 0;
 x2       = 0;
 Q        = 0;
 Part1Id  = 0;
 Part2Id  = 0;
 kfactor  = 0;
}

void MyGenKin::Print() {
  cout<<"ProcessId : "<<this->MCProcId<<endl;
  cout<<"PtHat     : "<<this->PtHat<<endl;
  cout<<"weight    : "<<this->genWeight<<endl;
  cout<<"x1        : "<<this->x1<<endl;
  cout<<"x2        : "<<this->x2<<endl;
  cout<<"Q         : "<<this->Q<<endl;
  cout<<"Part1Id   : "<<this->Part1Id<<endl;
  cout<<"Part2Id   : "<<this->Part2Id<<endl;
  cout<<"k Factor  : "<<this->kfactor<<endl;
}




