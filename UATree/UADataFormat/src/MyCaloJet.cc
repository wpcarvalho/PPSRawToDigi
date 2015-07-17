#include "MyCaloJet.h"
#include <iostream>

using namespace std;

ClassImp(MyCaloJet)

MyCaloJet::MyCaloJet():MyJet(){
  this->Reset();
}

MyCaloJet::~MyCaloJet(){}

void MyCaloJet::Print(){

  this->MyJet::Print();
  
  cout << "fem:       " << fem       << endl;  
  cout << "eem_EB:    " << eem_EB    << endl; 
  cout << "eem_EE:    " << eem_EE    << endl; 
  cout << "eem_HF:    " << eem_HF    << endl; 

  cout << "fhad:      " << fhad      << endl;   
  cout << "ehad_HB:   " << ehad_HB   << endl;
  cout << "ehad_HE:   " << ehad_HE   << endl;
  cout << "ehad_HF:   " << ehad_HF   << endl;
  cout << "ehad_HO:   " << ehad_HO   << endl;
  cout << "n60:       " << n60       << endl; 
  cout << "n90:       " << n90       << endl; 
  cout << "emax_ecal: " << emax_ecal << endl;
  cout << "emax_hcal: " << emax_hcal << endl;
  
  cout << "n90hits:   " << n90hits   << endl;
  cout << "HPD:       " << HPD	     << endl;
  cout << "RBX:       " << RBX	     << endl;
  cout << "sigma_eta: " << sigma_eta << endl;
  cout << "sigma_phi: " << sigma_phi << endl;
}

void MyCaloJet::Reset(){

  this->MyJet::Reset();
  
  fem        = 0;
  eem_EB     = 0;
  eem_EE     = 0;
  eem_HF     = 0;

  fhad       = 0;
  ehad_HB    = 0;
  ehad_HE    = 0;
  ehad_HF    = 0;
  ehad_HO    = 0;
  n60 	     = 0;
  n90 	     = 0;
  emax_ecal  = 0;
  emax_hcal  = 0;
  	     
  n90hits    = 0;
  HPD	     = 0;
  RBX	     = 0;
  sigma_eta  = 0;
  sigma_phi  = 0;
   
}
