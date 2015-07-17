/**
 * Class AlignParamtype.
 *
 * Author: Mirko Berretti / University of Siena 
 * Email:  mirko.berretti@gmail.com
 * Date:   2007-12-08
 */ //TrackInfo.h

#ifndef _TrackInfo_h
#define _TrackInfo_h
#include <vector>
#include <map>
#include <utility>
//#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/Framework/interface/EDProducer.h"
//#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"
//#include "DataFormats/T2Road/interface/T2Road.h"
#include "DataFormats/T2Cluster/interface/T2PadClusterCollection.h"
#include "DataFormats/T2Cluster/interface/T2StripClusterCollection.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/T2Road/interface/T2Road.h"
#include "DataFormats/T2Road/interface/T2RoadCollection.h"



class TrackInfo {
public:
  /// Constructor
  TrackInfo();
  /// Destructor
  virtual ~TrackInfo();
  

  struct IdandHit
  {
    T2Hit thehit;
    unsigned int iddet;
    double dr;
    double dphi;
//R-phi fit
    double OLDdr;
    double OLDdphi;
  };
 

  unsigned int eventnumber;
  unsigned int goodhitnumber;
  std::vector<IdandHit> idactive;
  //T1T2Track thetrack;
  double ar;
  double br;
  double bphi;
  double aphi;
  double chiR;
  //R-phi fit
  double OLDar;
  double OLDbr;
  double OLDbphi;
  double OLDaphi;
  

};
#endif
