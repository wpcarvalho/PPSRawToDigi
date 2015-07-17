#include "MyMet.h"
#include <iostream>
using namespace std;

ClassImp(MyMet)

MyMet::MyMet(){
  this->Reset();
}

MyMet::~MyMet(){}


void MyMet::Reset(){
  this->MyPart::Reset();
  sumet   = 0;
  elongit = 0;
}

void MyMet::Print(){
  
  this->MyPart::Print();

  cout << "sumet   : " << this->sumet   << endl;
  cout << "elongit : " << this->elongit << endl;
}
