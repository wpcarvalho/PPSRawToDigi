#include "MyDiJet.h"
#include <iostream>

using namespace std;

ClassImp(MyDiJet)

MyDiJet::MyDiJet(){
  this->Reset();
}

MyDiJet::~MyDiJet(){}

void MyDiJet::Reset() {
  isDiJet = 0;
  posJet1 = 0;
  posJet2 = 0;
}

void MyDiJet::Print() {
  cout<<"isDiJet: "<<this->isDiJet<<endl;
  cout<<"posJet1: "<<this->posJet1<<endl;
  cout<<"posJet2: "<<this->posJet2<<endl;
}
