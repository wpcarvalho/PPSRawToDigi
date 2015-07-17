#include "MyGenMet.h"
#include <iostream>
using namespace std;

ClassImp(MyGenMet)

MyGenMet::MyGenMet(){
  this->Reset();
}

MyGenMet::~MyGenMet(){}

void MyGenMet::Reset() {
 Met    = 0;
 MetX   = 0;
 MetY   = 0;
 MetPhi = 0;
 MetGP1 = TLorentzVector(0,0,0,0);
 MetGP3 = TLorentzVector(0,0,0,0);
}

void MyGenMet::Print() {
  cout << "From GenMet Collection :"    <<endl;
  cout << "  Met    : " << this->Met    <<endl;
  cout << "  MetX   : " << this->MetX   <<endl;
  cout << "  MetY   : " << this->MetY   <<endl;
  cout << "  MetPhi : " << this->MetPhi <<endl;
  cout << "From GenPart Collection (status 1):"   <<endl;
  cout << "  MetGP1.Et  : " << this->MetGP1.Et()  <<endl;
  cout << "  MetGP1.Eta : " << this->MetGP1.Eta() <<endl;
  cout << "  MetGP1.Phi : " << this->MetGP1.Phi() <<endl;
  cout << "From GenPart Collection (status 3):"   <<endl;
  cout << "  MetGP3.Et  : " << this->MetGP3.Et()  <<endl;
  cout << "  MetGP3.Eta : " << this->MetGP3.Eta() <<endl;
  cout << "  MetGP3.Phi : " << this->MetGP3.Phi() <<endl;
}




