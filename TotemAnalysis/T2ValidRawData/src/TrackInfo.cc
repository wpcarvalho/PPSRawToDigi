#include "TotemAnalysis/T2ValidRawData/interface/TrackInfo.h"

TrackInfo::TrackInfo()
{


 

  eventnumber=0;
  goodhitnumber=0;

  //  std::vector<IdandHit> idactive;
  //T1T2Track thetrack;
  ar=0.;
  br=0.;
  aphi=0.;
  bphi=0.;
  //For R-phi fit
  OLDar=0.;
  OLDbr=0.;
  OLDbphi=0.;
  OLDaphi=0.;

}

 


TrackInfo::~TrackInfo()
{}
