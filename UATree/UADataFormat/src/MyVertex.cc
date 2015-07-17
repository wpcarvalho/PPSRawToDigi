#include "MyVertex.h"
#include <iostream>

using namespace std;

ClassImp(MyVertex)

MyVertex::MyVertex(){
  this->Reset();
}

MyVertex::~MyVertex(){}


void MyVertex::Reset(){
  id          = 0;
  
  x	      = 0;
  y	      = 0;
  z	      = 0;

  ex	      = 0;
  ey	      = 0;
  ez	      = 0;

  validity    = 0;
  fake	      = 0;
  chi2	      = 0;
  ndof	      = 0;

  ntracks     = 0;
  SumPtTracks = 0;
}

void MyVertex::Print(){

  cout<<"vertex information: "<<endl;
  cout << "id         : " << this->id          <<endl;

  cout << "x          : " << this->x           <<endl;
  cout << "y          : " << this->y           <<endl;
  cout << "z          : " << this->z           <<endl;

  cout << "error x    : " << this->ex          <<endl;
  cout << "error y    : " << this->ey          <<endl;
  cout << "error z    : " << this->ez          <<endl;

  cout << "validity   : " << this->validity    <<endl;
  cout << "fake       : " << this->fake        <<endl;
  cout << "chi2       : " << this->chi2        <<endl;
  cout << "ndof       : " << this->ndof        <<endl;
  cout << "chi2n      : " << this->chi2n()     <<endl;

  cout << "ntracks    : " << this->ntracks     <<endl;
  cout << "SumPtTracks: " << this->SumPtTracks <<endl;
}


Double_t MyVertex::chi2n(){
  return ndof != 0 ? chi2 / ndof : -99;
}
