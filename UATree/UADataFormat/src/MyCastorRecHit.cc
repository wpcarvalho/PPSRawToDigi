#include "MyCastorRecHit.h"
#include <iostream>

using namespace std;

ClassImp(MyCastorRecHit)

MyCastorRecHit::MyCastorRecHit(){
  this->Reset();
}

MyCastorRecHit::~MyCastorRecHit(){}

void MyCastorRecHit::Print() {

  cout<<"castor rechit information: "<<endl;

  cout<<"module: "<<this->mod<<endl;
  cout<<"sector: "<<this->sec<<endl;
  cout<<"channel: "<<this->cha<<endl;
  cout<<"energy: "<<this->energy<<endl;
  cout<<"fC: "<<this->fC<<endl;
  cout<<"time: "<<this->time<<endl;
  cout<<"smearing: "<<this->smearing<<endl;
  cout<<"energy smeared: "<<this->energy_smeared<<endl;
  cout<<"fC smeared: "<<this->fC_smeared<<endl;
}

