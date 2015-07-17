#include "MyCaloTower.h"
#include <iostream>

ClassImp(MyCaloTower)

MyCaloTower::MyCaloTower(){
  this->Reset();
}

MyCaloTower::~MyCaloTower(){}

void MyCaloTower::Print(){
  MyPart::Print();
  
  cout << "emEnergy  : " << emEnergy  << endl;
  cout << "hadEnergy : " << hadEnergy << endl;

  cout << "hasEB : "<< hasEB << endl;
  cout << "hasEE : "<< hasEE << endl;
       
  cout << "hasHB : "<< hasHB << endl;
  cout << "hasHE : "<< hasHE << endl;
  cout << "hasHF : "<< hasHF << endl;

  cout << "zside : " << zside << endl;
}

void MyCaloTower::Reset(){
  MyPart::Reset();
  
  emEnergy  = 0;
  hadEnergy = 0;

  hasEB = false;
  hasEE = false;
  
  hasHB = false;
  hasHE = false;
  hasHF = false;

  zside= 0;
}
