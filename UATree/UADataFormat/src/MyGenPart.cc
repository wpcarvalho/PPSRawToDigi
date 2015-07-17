#include "MyGenPart.h"
#include <iostream>
using namespace std;

ClassImp(MyGenPart)

MyGenPart::MyGenPart():MyPart(){
  this->Reset();
}

MyGenPart::~MyGenPart(){}

void MyGenPart::Reset(){
  /*v      = TLorentzVector(0,0,0,0);
  charge = 0;*/
  
  this->MyPart::Reset();
  pdgId  = 0;
  status = 0;
  mo1    = -1;
  mo2    = -1;
  da1    = -1;
  da2    = -1;
}

void MyGenPart::Print(){
  /*cout << "px     : " << this->v.Px() << endl;
  cout << "py     : " << this->v.Py() << endl;
  cout << "pz     : " << this->v.Pz() << endl;
  cout << "E      : " << this->v.E()  << endl;
  cout << "charge : " << this->charge << endl;*/
  
  this->MyPart::Print();

  cout << "pdgId  : " << this->pdgId <<endl;
  cout << "status : " << this->status <<endl;
  //cout << "name: "<<this->name<<endl;
  
  cout << " mother positions   : " << this->mo1 << " and " << this->mo2 << endl;
  cout << " daughter positions : " << this->da1 << " and " << this->da2 << endl;
}
  


