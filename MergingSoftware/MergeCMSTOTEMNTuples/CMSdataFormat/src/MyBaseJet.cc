#include "MyBaseJet.h"
#include <iostream>

using namespace std;

ClassImp(MyBaseJet)

MyBaseJet::MyBaseJet():MyPart(){
  this->Reset();
}

MyBaseJet::~MyBaseJet(){}

void MyBaseJet::Print(){
  this->MyPart::Print();
  
  cout<<"correction - uncertainty: "<<endl;
  cout<<"jec: "<<this->jec<<endl;
  cout<<"jec_unc: "<<this->jec_unc<<endl;

}

void MyBaseJet::Reset(){
  this->MyPart::Reset();
 
  jec     = 0;
  jec_unc = 0;
  
}
