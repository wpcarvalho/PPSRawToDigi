#include "MyJet.h"
#include <iostream>

ClassImp(MyJet)

MyJet::MyJet(){
  this->Reset();
}

MyJet::~MyJet(){}

void MyJet::Print(){
  cout << "Printing mapjet :" << endl;
  
  for(map_it it = mapjet.begin() ; it != mapjet.end() ; ++it){
    cout << it->first << endl;
    it->second.Print(); 
  }
  
  cout << "nconstituent : " << nconstituent << endl;
  cout << "LooseJetId   : " << LooseJetId << endl;
  cout << "TightJetId   : " << TightJetId << endl;

}

void MyJet::Reset(){
  mapjet.clear();
  
  nconstituent = 0;
  LooseJetId   = 0;
  TightJetId   = 0;
}
