#include "MyPart.h"
#include <iostream>
using namespace std;

ClassImp(MyPart)

MyPart::MyPart(){
  this->Reset();
}

MyPart::~MyPart(){}


void MyPart::Reset(){
  this->SetPtEtaPhiM(0,0,0,0);
  charge = 0;
}

void MyPart::Print(){
  cout << "px     : " << this->Px()   << endl;
  cout << "py     : " << this->Py()   << endl;
  cout << "pz     : " << this->Pz()   << endl;
  cout << "E      : " << this->E()    << endl;
  cout << "Et     : " << this->Et()   << endl;
  cout << "pt     : " << this->Pt()   << endl;
  if(this->Pt()>0)
    cout << "eta    : " << this->Eta()  << endl;
  cout << "phi    : " << this->Phi()  << endl;
  cout << "mass   : " << this->M() << endl;
  cout << "charge : " << this->charge << endl;
}

TLorentzVector MyPart::vmpi()
{
  Double_t mpion = 139.57018;
  TLorentzVector vm;
  vm.SetPtEtaPhiM( this->Pt() , this->Eta() , this->Phi() , mpion );
  return vm;
}

Bool_t MyPart::operator<( const MyPart& part ){
  return this->Pt() < part.Pt();
}
 
//bool operator== ( const MyPart* genpart){
//  return (*this).pt == genpart.pt;
//} 
