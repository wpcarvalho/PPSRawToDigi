#include "MyGenJet.h"
#include <iostream>

using namespace std;

ClassImp(MyGenJet)

MyGenJet::MyGenJet():MyPart(){
  this->Reset();
}

MyGenJet::~MyGenJet(){}

void MyGenJet::Print(){

  cout << "[MyGenJet::Print()]" << endl;
  this->MyPart::Print();
  
  cout<<"npart : " << npart << endl;
  for(vector<MyGenPart>::iterator it = JetPart.begin() ; it != JetPart.end() ; ++it)
    it->Print();
}

void MyGenJet::Reset(){
  this->MyPart::Reset();
  
  npart   = 0;
  JetPart.clear();
  
}
