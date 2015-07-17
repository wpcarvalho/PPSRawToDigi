#include "MyPFCand.h"
#include <iostream>
using namespace std;

ClassImp(MyPFCand)

MyPFCand::MyPFCand(){
  this->Reset();
}

MyPFCand::~MyPFCand(){}


void MyPFCand::Reset(){
  this->MyPart::Reset();
  particleId = X;
}

void MyPFCand::Print(){
  this->MyPart::Print();
  cout << "particleId     : " << this->particleId << endl;
}
