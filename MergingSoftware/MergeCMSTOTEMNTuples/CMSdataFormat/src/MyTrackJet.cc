#include "MyTrackJet.h"
#include <iostream>

using namespace std;

ClassImp(MyTrackJet)

MyTrackJet::MyTrackJet():MyBaseJet(){
  this->Reset();
}

MyTrackJet::~MyTrackJet(){}

void MyTrackJet::Print(){
  this->MyBaseJet::Print();

  cout<<endl<<"track jet information: "<<endl;

  cout<<"raw jet variables: "<<endl;
  cout<<"energy: "<<this->e_raw<<endl;
  cout<<"pt: "<<this->pt_raw<<endl;
  cout<<"eta: "<<this->eta_raw<<endl;
  cout<<"phi: "<<this->phi_raw<<endl;
  cout<<"px: "<<this->px_raw<<endl;
  cout<<"py: "<<this->py_raw<<endl;
  cout<<"pz: "<<this->pz_raw<<endl;

  cout<<"correction - uncertainty: "<<endl;
  cout<<"jec: "<<this->jec<<endl;
  cout<<"jec_unc: "<<this->jec_unc<<endl;

  cout<<"corrected jet variables: "<<endl;
  cout<<"energy cal: "<<this->e_cal<<endl;
  cout<<"pt cal: "<<this->pt_cal<<endl;
  cout<<"eta cal: "<<this->eta_cal<<endl;
  cout<<"phi cal: "<<this->phi_cal<<endl;
  cout<<"px cal: "<<this->px_cal<<endl;
  cout<<"py cal: "<<this->py_cal<<endl;
  cout<<"pz cal: "<<this->pz_cal<<endl;

  cout<<"number of tracks: "<<this->ntrack<<endl;
  cout<<"jet associated to the hard primary vertex: "<<this->trackjet_pv<<endl;
  cout << "vtxId: " << vtxId << endl;


}

void MyTrackJet::Reset(){
  this->MyPart::Reset();
 
  vtxId = -1;

  e_raw = 0;
  pt_raw = 0;
  eta_raw = 0;
  phi_raw = 0;
  px_raw = 0;
  py_raw = 0;
  pz_raw = 0;

  e_cal = 0;
  pt_cal = 0;
  eta_cal = 0;
  phi_cal = 0;
  px_cal = 0;
  py_cal = 0;
  pz_cal = 0;

  ntrack = 0;
  trackjet_pv = 0;

  
}
