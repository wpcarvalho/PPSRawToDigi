#include "MyTracks.h"
#

ClassImp(MyTracks)

MyTracks::MyTracks():MyPart(){
  this->Reset();
}

MyTracks::~MyTracks(){}

void MyTracks::Print(){

  cout << "[MyTracks::Print()]" << endl;


  cout << "quality  ( enum , bool ) :" << endl;
  for(int i=0 ; i < 5 ; ++i)
    cout <<"  " << i << "  " << quality[i] << endl;
    
  cout << "trackAlgo         : " << trackAlgo           << endl;
  cout << "nhit              : " << nhit                << endl;
  cout << "nValidTrackerHits : " << nValidTrackerHits() << endl;
  cout << "nValidMuonHits    : " << nValidMuonHits()    << endl;
  cout << "nValidPixelHits   : " << nValidPixelHits     << endl;
  cout << "nValidStripHits   : " << nValidStripHits     << endl;  
  cout << "nValidMuonCSCHits : " << nValidMuonCSCHits   << endl;
  cout << "nValidMuonDTHits  : " << nValidMuonDTHits    << endl;
  cout << "nValidMuonRPCHits : " << nValidMuonRPCHits   << endl;

  cout << "chi2n : " << chi2n << endl;
  cout << "dz    : " << dz    << endl;  
  cout << "d0    : " << d0    << endl;  
  cout << "edz   : " << edz   << endl; 
  cout << "ed0   : " << ed0   << endl; 
  cout << "ept   : " << ept   << endl; 
  cout << "vx    : " << vx    << endl;  
  cout << "vy    : " << vy    << endl;  
  cout << "vz    : " << vz    << endl;  


  cout << " vtxid      vtxdxy      vtxdz :" << endl;
  for(unsigned int i=0 ; i < vtxid.size() ; ++i)
    cout << vtxid[i] << "   "  << vtxdxy[i] << "   " << vtxdz[i] << endl;
  
}

void MyTracks::Reset(){

  for(int i=0 ; i < 5 ; ++i)
    quality[i] = 0;
    
  trackAlgo		   = -1;
  nhit			   = 0;
  nValidPixelHits   = 0;
  nValidStripHits   = 0;
  nValidMuonCSCHits = 0;
  nValidMuonDTHits  = 0;
  nValidMuonRPCHits = 0;
  chi2n = 0;
  dz    = 0;
  d0    = 0;
  edz   = 0;
  ed0   = 0;
  ept   = 0;
  vx    = 0;
  vy    = 0;
  vz    = 0;
  
  vtxid.clear();
  vtxdxy.clear();
  vtxdz.clear();
  
}


