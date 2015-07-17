#include "MyCastorJet.h"
#include <iostream>

using namespace std;

ClassImp(MyCastorJet)

MyCastorJet::MyCastorJet(){
  this->Reset();
}

MyCastorJet::~MyCastorJet(){}

void MyCastorJet::Print() {

  cout<<"castor jet information: "<<endl;

  cout<<"energy: "<<this->energy<<endl;
  cout<<"pt: "<<this->pt<<endl;
  cout<<"eta: "<<this->eta<<endl;
  cout<<"phi: "<<this->phi<<endl;
  
  cout<<"fem: "<<this->fem<<endl;
  cout<<"eem: "<<this->eem<<endl;
  cout<<"ehad: "<<this->ehad<<endl;

  cout<<"width: "<<this->width<<endl;
  cout<<"depth: "<<this->depth<<endl;
  cout<<"fhot: "<<this->fhot<<endl;
  cout<<"sigmaz: "<<this->sigmaz<<endl;
  
  cout<<"ntower: "<<this->ntower<<endl;
}

