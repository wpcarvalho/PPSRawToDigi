#include "MyFwdGap.h"

ClassImp(MyFwdGap)

MyFwdGap::MyFwdGap(){
  this->Reset();
}

MyFwdGap::~MyFwdGap(){}

void MyFwdGap::Print(){}

void MyFwdGap::Reset()
{
  nTowersHF_plus = 0;
  nTowersHF_minus = 0;
  nTowersHE_plus = 0;
  nTowersHE_minus = 0;
  nTowersHB_plus = 0;
  nTowersHB_minus = 0;
  nTowersEE_plus = 0;
  nTowersEE_minus = 0;
  nTowersEB_plus = 0;
  nTowersEB_minus = 0; 
  //Sum(E)
  sumEHF_plus = 0.;
  sumEHF_minus = 0.;
  sumEHE_plus = 0.;
  sumEHE_minus = 0.;
  sumEHB_plus = 0.;
  sumEHB_minus = 0.;
  sumEEE_plus = 0.;
  sumEEE_minus = 0.;
  sumEEB_plus = 0.;
  sumEEB_minus = 0.;
  // Sum(ET)
  sumETHF_plus = 0.;
  sumETHF_minus = 0.;
  sumETHE_plus = 0.;
  sumETHE_minus = 0.;
  sumETHB_plus = 0.;
  sumETHB_minus = 0.;
  sumETEE_plus = 0.;
  sumETEE_minus = 0.;
  sumETEB_plus = 0.;
  sumETEB_minus = 0.;
}

