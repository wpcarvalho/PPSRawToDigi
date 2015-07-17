#include "MyZDCDigi.h"
#include <iostream>
#include <sstream>

using namespace std;

ClassImp(MyZDCDigi)

MyZDCDigi::MyZDCDigi(){
  this->Reset();
}

MyZDCDigi::~MyZDCDigi(){}

void MyZDCDigi::Reset() {
  side = 0; 
  section = -1;
  channel = -1;
  channelId = -1;
  digiADC.clear();
  digifC.clear();
}

void MyZDCDigi::Print() {
  stringstream oss;
  oss << "side: " << this->side << endl
      << "section: " << this->section << endl
      << "channel: " << this->channel << endl
      << "channelId: " << this->channelId << endl;
  for(size_t iTS = 0; iTS < digiADC.size(); ++iTS)
     oss << "   TS: " << iTS << "  ADC: " << digiADC[iTS] << "  charge (fC): " << digifC[iTS] << endl;

  cout << oss.str();
}
