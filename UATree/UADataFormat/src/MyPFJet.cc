#include "MyPFJet.h"
#include <iostream>

using namespace std;

ClassImp(MyPFJet)

MyPFJet::MyPFJet():MyJet(){
  this->Reset();
}

MyPFJet::~MyPFJet(){}

void MyPFJet::Print(){

  this->MyJet::Print();
  
  cout << "fhad_ch      : " << fhad_ch	    << endl;  
  cout << "fhad_ne      : " << fhad_ne	    << endl; 
  
  cout << "fem_ch       : " << fem_ch	    << endl; 
  cout << "fem_ne	: " << fem_ne       << endl; 
  cout << "multi_ch     : " << multi_ch	    << endl;   
  cout << "multi_ne     : " << multi_ne	    << endl;
  
  cout << "multi_ch_had : " << multi_ch_had << endl;
  cout << "multi_ne_had : " << multi_ne_had << endl;
  
  cout << "multi_gamma  : " << multi_gamma  << endl; 
  cout << "multi_ele    : " << multi_ele    << endl;
  cout << "multi_mu     : " << multi_mu	    << endl;
  
}

void MyPFJet::Reset(){

  this->MyJet::Reset();
  
  fhad_ch       = 0;
  fhad_ne       = 0;
  
  fem_ch        = 0;
  fem_ne        = 0;
  multi_ch      = 0;
  multi_ne      = 0;
  
  multi_ch_had  = 0;
  multi_ne_had  = 0;
  
  multi_gamma   = 0;
  multi_ele     = 0;
  multi_mu      = 0;

  ntrack = 0;
  vtracks.clear();
}



