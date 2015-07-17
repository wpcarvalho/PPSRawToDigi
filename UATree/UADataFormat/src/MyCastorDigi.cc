
#include "MyCastorDigi.h"
#include <iostream>

using namespace std;

ClassImp(MyCastorDigi)

MyCastorDigi::MyCastorDigi(){
  this->Reset();
}

MyCastorDigi::~MyCastorDigi(){}

void MyCastorDigi::Reset() {

  this->adc.clear();
  this->fC.clear();

  this->mod = 0;
  this->sec = 0;
  this->cha = 0;
}


void MyCastorDigi::Print() {

  cout<<"castor digi information: "<<endl;

  cout<<"module: "<<this->mod<<endl;
  cout<<"sector: "<<this->sec<<endl;
  cout<<"channel: "<<this->cha<<endl;

  Int_t ts = 1;
  for(vector<Double_t>::const_iterator it_adc = this->adc.begin(); it_adc != this->adc.end(); ++it_adc) 
    cout<<"time slice : "<<ts++<<"adc: "<<*it_adc<<endl;
  
  ts = 1;
  for(vector<Double_t>::const_iterator it_fC = this->fC.begin(); it_fC != this->fC.end(); ++it_fC) 
    cout<<"time slice : "<<ts++<<"fC: "<<*it_fC<<endl;
}

