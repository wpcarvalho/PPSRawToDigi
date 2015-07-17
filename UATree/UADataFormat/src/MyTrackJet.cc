#include "MyTrackJet.h"
#include <iostream>

using namespace std;

ClassImp(MyTrackJet)

MyTrackJet::MyTrackJet():MyBaseJet(){
  this->Reset();
}

MyTrackJet::~MyTrackJet(){}

void MyTrackJet::Print(){
  this->MyBaseJet::Print();

  cout << "vtxId: " << vtxId << endl;
}

void MyTrackJet::Reset(){
  this->MyPart::Reset();
 
  vtxId = -1;
  
}
