#include "MyFSCHit.h"
#include <iostream>

using namespace std;

ClassImp(MyFSCHit)

MyFSCHit::MyFSCHit(){
  this->Reset();
}

MyFSCHit::~MyFSCHit(){}

void MyFSCHit::Reset() {
  side = 0; 
  section = -1;
  channel = -1;
  channelId = -1;
  energy = -1.;
  time = -1.;
}

void MyFSCHit::Print() {
  cout << "side: " << this->side << endl;
  cout << "section: " << this->section << endl;
  cout << "channel: " << this->channel << endl;
  cout << "channelId: " << this->channelId << endl;
  cout << "energy: " << this->energy << endl;
  cout << "time: " << this->time << endl;
}
