#include "MyZDCInfo.h"
#include <iostream>

using namespace std;

ClassImp(MyZDCInfo)

MyZDCInfo::MyZDCInfo(){
  this->Reset();
}

MyZDCInfo::~MyZDCInfo(){}

void MyZDCInfo::Reset() {
  /*nHits             = 0;
  sumEnergyEMPlus   = 0.;
  sumEnergyEMMinus  = 0.;
  sumEnergyHADPlus  = 0.;
  sumEnergyHADMinus = 0.;*/
  nHitsPerChannel.clear();
  sumEnergyPerChannel.clear();
}

void MyZDCInfo::Print() {
  /*cout << "nHits: " << this->nHits << endl;
  cout << "sumEnergyEMPlus: " << this->sumEnergyEMPlus << endl;
  cout << "sumEnergyEMMinus: " << this->sumEnergyEMMinus << endl;
  cout << "sumEnergyHADPlus: " << this->sumEnergyHADPlus << endl;
  cout << "sumEnergyHADMinus: " << this->sumEnergyHADMinus << endl;*/
  std::map<int,int>::const_iterator it_channel = nHitsPerChannel.begin();
  for(; it_channel != nHitsPerChannel.end(); ++it_channel){
     int channelId    = it_channel->first;
     int nHits        = it_channel->second;
     double sumEnergy = sumEnergyPerChannel[channelId]; 
     cout << "channel Id: " << channelId << endl
          << "     nHits: " << nHits << endl
          << "    energy: " << sumEnergy << endl;
  } 
}
