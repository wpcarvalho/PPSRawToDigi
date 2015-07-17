#include "MyPUSumInfo.h"
#include <iostream>

using namespace std;

ClassImp(MyPUSumInfo)

MyPUSumInfo::MyPUSumInfo(){
  this->Reset();
}

MyPUSumInfo::~MyPUSumInfo(){}


void MyPUSumInfo::Reset(){

  nPU = -1 ;
  zposition.clear();
  sumpT_lowpT.clear();
  sumpT_highpT.clear();
  ntrks_lowpT.clear();  
  ntrks_highpT.clear();

}

void MyPUSumInfo::Print(){

  cout<<"MC PU information: "<<endl;

  cout << "nPU    : " << this->nPU     <<endl;
  if ( nPU >= 0 ) {
    for ( int ivtx = 0; ivtx<nPU ; ++ivtx) {
      cout << "PileUp event #    : " << ivtx     << endl;
      cout << "  --> zposition   : " << zposition.at(ivtx) << endl;
      cout << "  --> sumpT_lowpT : " << sumpT_lowpT.at(ivtx) << endl;
      cout << "  --> sumpT_highpT: " << sumpT_highpT.at(ivtx) << endl;
      cout << "  --> ntrks_lowpT : " << ntrks_lowpT.at(ivtx) << endl;
      cout << "  --> ntrks_highpT: " << ntrks_highpT.at(ivtx) << endl;
    }
  }

}


