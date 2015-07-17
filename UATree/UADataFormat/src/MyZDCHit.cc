#include "MyZDCHit.h"
#include <iostream>

using namespace std;

ClassImp(MyZDCHit)

MyZDCHit::MyZDCHit(){
  this->Reset();
}

MyZDCHit::~MyZDCHit(){}

void MyZDCHit::Reset() {
  side = 0; 
  section = -1;
  channel = -1;
  channelId = -1;
  energy = -1.;
  time = -1.;
}

void MyZDCHit::Print() {
  cout << "side: " << this->side << endl;
  cout << "section: " << this->section << endl;
  cout << "channel: " << this->channel << endl;
  cout << "channelId: " << this->channelId << endl;
  cout << "energy: " << this->energy << endl;
  cout << "time: " << this->time << endl;
}
