// -*- C++ -*-
//
// Package:    T1Validation
// Class:      T1Validation
//
// Author> F.Ferro - INFN Genova

#ifndef _TotemT1T2ValidationT1ValidationOnlyDataT1ValidationOnlyData_H_
#define _TotemT1T2ValidationT1ValidationOnlyDataT1ValidationOnlyData_H_
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "DataFormats/T1DigiWire/interface/T1DigiWireCollection.h"
#include "DataFormats/T1DigiVfat/interface/T1DigiVfatCollection.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include "DataFormats/T1Road/interface/T1RecHitGlobal.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"
#include "DataFormats/T1Cluster/interface/T1Cluster.h"
#include "DataFormats/T1Cluster/interface/T1ClusterCollection.h"
#include "DataFormats/T1Road/interface/T1Road.h"
#include <boost/shared_ptr.hpp>
#include "TFile.h"
#include "TH1.h"

//
// class decleration
//

class T1ValidationOnlyData : public edm::EDAnalyzer {
public:
  explicit T1ValidationOnlyData(const edm::ParameterSet&);
  ~T1ValidationOnlyData();


private:
  edm::InputTag t1DigiWireCollectionLabel;
  edm::InputTag t1DigiVfatCollectionLabel;
  edm::InputTag t1ClusterCollectionLabel;
  edm::InputTag t1RoadCollectionLabel;
  edm::InputTag rawEventLabel;
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void initHistograms();
  void writeHistograms();

  double _SIM_;
  double _DIGI_;
  double _RECO_;
  int MaxTracks;
  int MinTracks;

  std::vector<PSimHit> *copy_of_hits;


  std::string trackLabel;
  std::string outputFileName;

  //dai digi
  std::auto_ptr<TH1D> ilWire;
//  std::auto_ptr<TH1D> hClusterWidth;
   std::auto_ptr<TH1D>hClDiffBetweenLayers;

  std::auto_ptr<TH1D>hClRMSSex;
  std::auto_ptr<TH1D> ilWireCh0Arm0Pl0;
  std::auto_ptr<TH1D> ilWireCh0Arm0Pl1;
  std::auto_ptr<TH1D> ilWireCh0Arm0Pl2;
  std::auto_ptr<TH1D> ilWireCh0Arm0Pl3;
  std::auto_ptr<TH1D> ilWireCh0Arm0Pl4;
  std::auto_ptr<TH1D> hRoadSize;
  std::auto_ptr<TH1D> hRoadNumber;
  std::auto_ptr<TH1D> hRoadHitSigmaX;
  std::auto_ptr<TH1D> hRoadHitSigmaY;
  std::auto_ptr<TH1D> ilWireCh0Arm1Pl0;
  std::auto_ptr<TH1D> ilWireCh0Arm1Pl1;
  std::auto_ptr<TH1D> ilWireCh0Arm1Pl2;
  std::auto_ptr<TH1D> ilWireCh0Arm1Pl3;
  std::auto_ptr<TH1D> ilWireCh0Arm1Pl4;

  std::auto_ptr<TH1D> ilClWidthCh0Arm0Pl0;
  std::auto_ptr<TH1D> ilClWidthCh0Arm0Pl1;
  std::auto_ptr<TH1D> ilClWidthCh0Arm0Pl2;
  std::auto_ptr<TH1D> ilClWidthCh0Arm0Pl3;
  std::auto_ptr<TH1D> ilClWidthCh0Arm0Pl4;

  std::auto_ptr<TH1D> ilClWidthCh0Arm1Pl0;
  std::auto_ptr<TH1D> ilClWidthCh0Arm1Pl1;
  std::auto_ptr<TH1D> ilClWidthCh0Arm1Pl2;
  std::auto_ptr<TH1D> ilClWidthCh0Arm1Pl3;
  std::auto_ptr<TH1D> ilClWidthCh0Arm1Pl4;

 std::auto_ptr<TH1D> ilClWidthCh0Arm1Pl3_0_20;
 std::auto_ptr<TH1D> ilClWidthCh0Arm1Pl3_20_40;
 std::auto_ptr<TH1D> ilClWidthCh0Arm1Pl3_40_60;
 std::auto_ptr<TH1D> ilClWidthCh0Arm1Pl3_60_80;
 std::auto_ptr<TH1D> ilClWidthCh0Arm1Pl3_80_100;
 std::auto_ptr<TH1D> ilClWidthCh0Arm1Pl3_100_120;
 std::auto_ptr<TH1D> ilClWidthCh0Arm1Pl3_120_140;

 std::auto_ptr<TH1D> ilClWidthCh5Arm1Pl3_0_20;
 std::auto_ptr<TH1D> ilClWidthCh5Arm1Pl3_20_40;
 std::auto_ptr<TH1D> ilClWidthCh5Arm1Pl3_40_60;
 std::auto_ptr<TH1D> ilClWidthCh5Arm1Pl3_60_80;
 std::auto_ptr<TH1D> ilClWidthCh5Arm1Pl3_80_100;
 std::auto_ptr<TH1D> ilClWidthCh5Arm1Pl3_100_120;
 std::auto_ptr<TH1D> ilClWidthCh5Arm1Pl3_120_140;



  std::auto_ptr<TH1D> laStripACh0Arm0Pl0;
  std::auto_ptr<TH1D> laStripACh0Arm0Pl1;
  std::auto_ptr<TH1D> laStripACh0Arm0Pl2;
  std::auto_ptr<TH1D> laStripACh0Arm0Pl3;
  std::auto_ptr<TH1D> laStripACh0Arm0Pl4;

  std::auto_ptr<TH1D> laStripACh0Arm1Pl0;
  std::auto_ptr<TH1D> laStripACh0Arm1Pl1;
  std::auto_ptr<TH1D> laStripACh0Arm1Pl2;
  std::auto_ptr<TH1D> laStripACh0Arm1Pl3;
  std::auto_ptr<TH1D> laStripACh0Arm1Pl4;

  std::auto_ptr<TH1D> laStripBCh0Arm0Pl0;
  std::auto_ptr<TH1D> laStripBCh0Arm0Pl1;
  std::auto_ptr<TH1D> laStripBCh0Arm0Pl2;
  std::auto_ptr<TH1D> laStripBCh0Arm0Pl3;
  std::auto_ptr<TH1D> laStripBCh0Arm0Pl4;

  std::auto_ptr<TH1D> laStripBCh0Arm1Pl0;
  std::auto_ptr<TH1D> laStripBCh0Arm1Pl1;
  std::auto_ptr<TH1D> laStripBCh0Arm1Pl2;
  std::auto_ptr<TH1D> laStripBCh0Arm1Pl3;
  std::auto_ptr<TH1D> laStripBCh0Arm1Pl4;

  std::auto_ptr<TH1D> ilWireCh1Arm0Pl0;
  std::auto_ptr<TH1D> ilWireCh1Arm0Pl1;
  std::auto_ptr<TH1D> ilWireCh1Arm0Pl2;
  std::auto_ptr<TH1D> ilWireCh1Arm0Pl3;
  std::auto_ptr<TH1D> ilWireCh1Arm0Pl4;

  std::auto_ptr<TH1D> ilWireCh1Arm1Pl0;
  std::auto_ptr<TH1D> ilWireCh1Arm1Pl1;
  std::auto_ptr<TH1D> ilWireCh1Arm1Pl2;
  std::auto_ptr<TH1D> ilWireCh1Arm1Pl3;
  std::auto_ptr<TH1D> ilWireCh1Arm1Pl4;

  std::auto_ptr<TH1D> ilClWidthCh1Arm0Pl0;
  std::auto_ptr<TH1D> ilClWidthCh1Arm0Pl1;
  std::auto_ptr<TH1D> ilClWidthCh1Arm0Pl2;
  std::auto_ptr<TH1D> ilClWidthCh1Arm0Pl3;
  std::auto_ptr<TH1D> ilClWidthCh1Arm0Pl4;

  std::auto_ptr<TH1D> ilClWidthCh1Arm1Pl0;
  std::auto_ptr<TH1D> ilClWidthCh1Arm1Pl1;
  std::auto_ptr<TH1D> ilClWidthCh1Arm1Pl2;
  std::auto_ptr<TH1D> ilClWidthCh1Arm1Pl3;
  std::auto_ptr<TH1D> ilClWidthCh1Arm1Pl4;

  std::auto_ptr<TH1D> laStripACh1Arm0Pl0;
  std::auto_ptr<TH1D> laStripACh1Arm0Pl1;
  std::auto_ptr<TH1D> laStripACh1Arm0Pl2;
  std::auto_ptr<TH1D> laStripACh1Arm0Pl3;
  std::auto_ptr<TH1D> laStripACh1Arm0Pl4;

  std::auto_ptr<TH1D> laStripACh1Arm1Pl0;
  std::auto_ptr<TH1D> laStripACh1Arm1Pl1;
  std::auto_ptr<TH1D> laStripACh1Arm1Pl2;
  std::auto_ptr<TH1D> laStripACh1Arm1Pl3;
  std::auto_ptr<TH1D> laStripACh1Arm1Pl4;

  std::auto_ptr<TH1D> laStripBCh1Arm0Pl0;
  std::auto_ptr<TH1D> laStripBCh1Arm0Pl1;
  std::auto_ptr<TH1D> laStripBCh1Arm0Pl2;
  std::auto_ptr<TH1D> laStripBCh1Arm0Pl3;
  std::auto_ptr<TH1D> laStripBCh1Arm0Pl4;

  std::auto_ptr<TH1D> laStripBCh1Arm1Pl0;
  std::auto_ptr<TH1D> laStripBCh1Arm1Pl1;
  std::auto_ptr<TH1D> laStripBCh1Arm1Pl2;
  std::auto_ptr<TH1D> laStripBCh1Arm1Pl3;
  std::auto_ptr<TH1D> laStripBCh1Arm1Pl4;


  std::auto_ptr<TH1D> ilWireCh2Arm0Pl0;
  std::auto_ptr<TH1D> ilWireCh2Arm0Pl1;
  std::auto_ptr<TH1D> ilWireCh2Arm0Pl2;
  std::auto_ptr<TH1D> ilWireCh2Arm0Pl3;
  std::auto_ptr<TH1D> ilWireCh2Arm0Pl4;

  std::auto_ptr<TH1D> ilWireCh2Arm1Pl0;
  std::auto_ptr<TH1D> ilWireCh2Arm1Pl1;
  std::auto_ptr<TH1D> ilWireCh2Arm1Pl2;
  std::auto_ptr<TH1D> ilWireCh2Arm1Pl3;
  std::auto_ptr<TH1D> ilWireCh2Arm1Pl4;

  std::auto_ptr<TH1D> ilClWidthCh2Arm0Pl0;
  std::auto_ptr<TH1D> ilClWidthCh2Arm0Pl1;
  std::auto_ptr<TH1D> ilClWidthCh2Arm0Pl2;
  std::auto_ptr<TH1D> ilClWidthCh2Arm0Pl3;
  std::auto_ptr<TH1D> ilClWidthCh2Arm0Pl4;

  std::auto_ptr<TH1D> ilClWidthCh2Arm1Pl0;
  std::auto_ptr<TH1D> ilClWidthCh2Arm1Pl1;
  std::auto_ptr<TH1D> ilClWidthCh2Arm1Pl2;
  std::auto_ptr<TH1D> ilClWidthCh2Arm1Pl3;
  std::auto_ptr<TH1D> ilClWidthCh2Arm1Pl4;

  std::auto_ptr<TH1D> laStripACh2Arm0Pl0;
  std::auto_ptr<TH1D> laStripACh2Arm0Pl1;
  std::auto_ptr<TH1D> laStripACh2Arm0Pl2;
  std::auto_ptr<TH1D> laStripACh2Arm0Pl3;
  std::auto_ptr<TH1D> laStripACh2Arm0Pl4;

  std::auto_ptr<TH1D> laStripACh2Arm1Pl0;
  std::auto_ptr<TH1D> laStripACh2Arm1Pl1;
  std::auto_ptr<TH1D> laStripACh2Arm1Pl2;
  std::auto_ptr<TH1D> laStripACh2Arm1Pl3;
  std::auto_ptr<TH1D> laStripACh2Arm1Pl4;

  std::auto_ptr<TH1D> laStripBCh2Arm0Pl0;
  std::auto_ptr<TH1D> laStripBCh2Arm0Pl1;
  std::auto_ptr<TH1D> laStripBCh2Arm0Pl2;
  std::auto_ptr<TH1D> laStripBCh2Arm0Pl3;
  std::auto_ptr<TH1D> laStripBCh2Arm0Pl4;

  std::auto_ptr<TH1D> laStripBCh2Arm1Pl0;
  std::auto_ptr<TH1D> laStripBCh2Arm1Pl1;
  std::auto_ptr<TH1D> laStripBCh2Arm1Pl2;
  std::auto_ptr<TH1D> laStripBCh2Arm1Pl3;
  std::auto_ptr<TH1D> laStripBCh2Arm1Pl4;


  std::auto_ptr<TH1D> ilWireCh3Arm0Pl0;
  std::auto_ptr<TH1D> ilWireCh3Arm0Pl1;
  std::auto_ptr<TH1D> ilWireCh3Arm0Pl2;
  std::auto_ptr<TH1D> ilWireCh3Arm0Pl3;
  std::auto_ptr<TH1D> ilWireCh3Arm0Pl4;

  std::auto_ptr<TH1D> ilWireCh3Arm1Pl0;
  std::auto_ptr<TH1D> ilWireCh3Arm1Pl1;
  std::auto_ptr<TH1D> ilWireCh3Arm1Pl2;
  std::auto_ptr<TH1D> ilWireCh3Arm1Pl3;
  std::auto_ptr<TH1D> ilWireCh3Arm1Pl4;

  std::auto_ptr<TH1D> ilClWidthCh3Arm0Pl0;
  std::auto_ptr<TH1D> ilClWidthCh3Arm0Pl1;
  std::auto_ptr<TH1D> ilClWidthCh3Arm0Pl2;
  std::auto_ptr<TH1D> ilClWidthCh3Arm0Pl3;
  std::auto_ptr<TH1D> ilClWidthCh3Arm0Pl4;

  std::auto_ptr<TH1D> ilClWidthCh3Arm1Pl0;
  std::auto_ptr<TH1D> ilClWidthCh3Arm1Pl1;
  std::auto_ptr<TH1D> ilClWidthCh3Arm1Pl2;
  std::auto_ptr<TH1D> ilClWidthCh3Arm1Pl3;
  std::auto_ptr<TH1D> ilClWidthCh3Arm1Pl4;

  std::auto_ptr<TH1D> laStripACh3Arm0Pl0;
  std::auto_ptr<TH1D> laStripACh3Arm0Pl1;
  std::auto_ptr<TH1D> laStripACh3Arm0Pl2;
  std::auto_ptr<TH1D> laStripACh3Arm0Pl3;
  std::auto_ptr<TH1D> laStripACh3Arm0Pl4;

  std::auto_ptr<TH1D> laStripACh3Arm1Pl0;
  std::auto_ptr<TH1D> laStripACh3Arm1Pl1;
  std::auto_ptr<TH1D> laStripACh3Arm1Pl2;
  std::auto_ptr<TH1D> laStripACh3Arm1Pl3;
  std::auto_ptr<TH1D> laStripACh3Arm1Pl4;

  std::auto_ptr<TH1D> laStripBCh3Arm0Pl0;
  std::auto_ptr<TH1D> laStripBCh3Arm0Pl1;
  std::auto_ptr<TH1D> laStripBCh3Arm0Pl2;
  std::auto_ptr<TH1D> laStripBCh3Arm0Pl3;
  std::auto_ptr<TH1D> laStripBCh3Arm0Pl4;

  std::auto_ptr<TH1D> laStripBCh3Arm1Pl0;
  std::auto_ptr<TH1D> laStripBCh3Arm1Pl1;
  std::auto_ptr<TH1D> laStripBCh3Arm1Pl2;
  std::auto_ptr<TH1D> laStripBCh3Arm1Pl3;
  std::auto_ptr<TH1D> laStripBCh3Arm1Pl4;


  std::auto_ptr<TH1D> ilWireCh4Arm0Pl0;
  std::auto_ptr<TH1D> ilWireCh4Arm0Pl1;
  std::auto_ptr<TH1D> ilWireCh4Arm0Pl2;
  std::auto_ptr<TH1D> ilWireCh4Arm0Pl3;
  std::auto_ptr<TH1D> ilWireCh4Arm0Pl4;

  std::auto_ptr<TH1D> ilWireCh4Arm1Pl0;
  std::auto_ptr<TH1D> ilWireCh4Arm1Pl1;
  std::auto_ptr<TH1D> ilWireCh4Arm1Pl2;
  std::auto_ptr<TH1D> ilWireCh4Arm1Pl3;
  std::auto_ptr<TH1D> ilWireCh4Arm1Pl4;


  std::auto_ptr<TH1D> ilClWidthCh4Arm0Pl0;
  std::auto_ptr<TH1D> ilClWidthCh4Arm0Pl1;
  std::auto_ptr<TH1D> ilClWidthCh4Arm0Pl2;
  std::auto_ptr<TH1D> ilClWidthCh4Arm0Pl3;
  std::auto_ptr<TH1D> ilClWidthCh4Arm0Pl4;

  std::auto_ptr<TH1D> ilClWidthCh4Arm1Pl0;
  std::auto_ptr<TH1D> ilClWidthCh4Arm1Pl1;
  std::auto_ptr<TH1D> ilClWidthCh4Arm1Pl2;
  std::auto_ptr<TH1D> ilClWidthCh4Arm1Pl3;
  std::auto_ptr<TH1D> ilClWidthCh4Arm1Pl4;

  std::auto_ptr<TH1D> laStripACh4Arm0Pl0;
  std::auto_ptr<TH1D> laStripACh4Arm0Pl1;
  std::auto_ptr<TH1D> laStripACh4Arm0Pl2;
  std::auto_ptr<TH1D> laStripACh4Arm0Pl3;
  std::auto_ptr<TH1D> laStripACh4Arm0Pl4;

  std::auto_ptr<TH1D> laStripACh4Arm1Pl0;
  std::auto_ptr<TH1D> laStripACh4Arm1Pl1;
  std::auto_ptr<TH1D> laStripACh4Arm1Pl2;
  std::auto_ptr<TH1D> laStripACh4Arm1Pl3;
  std::auto_ptr<TH1D> laStripACh4Arm1Pl4;

  std::auto_ptr<TH1D> laStripBCh4Arm0Pl0;
  std::auto_ptr<TH1D> laStripBCh4Arm0Pl1;
  std::auto_ptr<TH1D> laStripBCh4Arm0Pl2;
  std::auto_ptr<TH1D> laStripBCh4Arm0Pl3;
  std::auto_ptr<TH1D> laStripBCh4Arm0Pl4;

  std::auto_ptr<TH1D> laStripBCh4Arm1Pl0;
  std::auto_ptr<TH1D> laStripBCh4Arm1Pl1;
  std::auto_ptr<TH1D> laStripBCh4Arm1Pl2;
  std::auto_ptr<TH1D> laStripBCh4Arm1Pl3;
  std::auto_ptr<TH1D> laStripBCh4Arm1Pl4;


  std::auto_ptr<TH1D> ilWireCh5Arm0Pl0;
  std::auto_ptr<TH1D> ilWireCh5Arm0Pl1;
  std::auto_ptr<TH1D> ilWireCh5Arm0Pl2;
  std::auto_ptr<TH1D> ilWireCh5Arm0Pl3;
  std::auto_ptr<TH1D> ilWireCh5Arm0Pl4;

  std::auto_ptr<TH1D> ilWireCh5Arm1Pl0;
  std::auto_ptr<TH1D> ilWireCh5Arm1Pl1;
  std::auto_ptr<TH1D> ilWireCh5Arm1Pl2;
  std::auto_ptr<TH1D> ilWireCh5Arm1Pl3;
  std::auto_ptr<TH1D> ilWireCh5Arm1Pl4;


  std::auto_ptr<TH1D> ilClWidthCh5Arm0Pl0;
  std::auto_ptr<TH1D> ilClWidthCh5Arm0Pl1;
  std::auto_ptr<TH1D> ilClWidthCh5Arm0Pl2;
  std::auto_ptr<TH1D> ilClWidthCh5Arm0Pl3;
  std::auto_ptr<TH1D> ilClWidthCh5Arm0Pl4;

  std::auto_ptr<TH1D> ilClWidthCh5Arm1Pl0;
  std::auto_ptr<TH1D> ilClWidthCh5Arm1Pl1;
  std::auto_ptr<TH1D> ilClWidthCh5Arm1Pl2;
  std::auto_ptr<TH1D> ilClWidthCh5Arm1Pl3;
  std::auto_ptr<TH1D> ilClWidthCh5Arm1Pl4;

  std::auto_ptr<TH1D> laStripACh5Arm0Pl0;
  std::auto_ptr<TH1D> laStripACh5Arm0Pl1;
  std::auto_ptr<TH1D> laStripACh5Arm0Pl2;
  std::auto_ptr<TH1D> laStripACh5Arm0Pl3;
  std::auto_ptr<TH1D> laStripACh5Arm0Pl4;

  std::auto_ptr<TH1D> laStripACh5Arm1Pl0;
  std::auto_ptr<TH1D> laStripACh5Arm1Pl1;
  std::auto_ptr<TH1D> laStripACh5Arm1Pl2;
  std::auto_ptr<TH1D> laStripACh5Arm1Pl3;
  std::auto_ptr<TH1D> laStripACh5Arm1Pl4;

  std::auto_ptr<TH1D> laStripBCh5Arm0Pl0;
  std::auto_ptr<TH1D> laStripBCh5Arm0Pl1;
  std::auto_ptr<TH1D> laStripBCh5Arm0Pl2;
  std::auto_ptr<TH1D> laStripBCh5Arm0Pl3;
  std::auto_ptr<TH1D> laStripBCh5Arm0Pl4;

  std::auto_ptr<TH1D> laStripBCh5Arm1Pl0;
  std::auto_ptr<TH1D> laStripBCh5Arm1Pl1;
  std::auto_ptr<TH1D> laStripBCh5Arm1Pl2;
  std::auto_ptr<TH1D> laStripBCh5Arm1Pl3;
  std::auto_ptr<TH1D> laStripBCh5Arm1Pl4;

///////////////////////////////////

  std::auto_ptr<TH1D> NWireCh0Arm0Pl0;
  std::auto_ptr<TH1D> NWireCh0Arm0Pl1;
  std::auto_ptr<TH1D> NWireCh0Arm0Pl2;
  std::auto_ptr<TH1D> NWireCh0Arm0Pl3;
  std::auto_ptr<TH1D> NWireCh0Arm0Pl4;

  std::auto_ptr<TH1D> NWireCh0Arm1Pl0;
  std::auto_ptr<TH1D> NWireCh0Arm1Pl1;
  std::auto_ptr<TH1D> NWireCh0Arm1Pl2;
  std::auto_ptr<TH1D> NWireCh0Arm1Pl3;
  std::auto_ptr<TH1D> NWireCh0Arm1Pl4;

  std::auto_ptr<TH1D> NStripACh0Arm0Pl0;
  std::auto_ptr<TH1D> NStripACh0Arm0Pl1;
  std::auto_ptr<TH1D> NStripACh0Arm0Pl2;
  std::auto_ptr<TH1D> NStripACh0Arm0Pl3;
  std::auto_ptr<TH1D> NStripACh0Arm0Pl4;

  std::auto_ptr<TH1D> NStripACh0Arm1Pl0;
  std::auto_ptr<TH1D> NStripACh0Arm1Pl1;
  std::auto_ptr<TH1D> NStripACh0Arm1Pl2;
  std::auto_ptr<TH1D> NStripACh0Arm1Pl3;
  std::auto_ptr<TH1D> NStripACh0Arm1Pl4;

  std::auto_ptr<TH1D> NStripBCh0Arm0Pl0;
  std::auto_ptr<TH1D> NStripBCh0Arm0Pl1;
  std::auto_ptr<TH1D> NStripBCh0Arm0Pl2;
  std::auto_ptr<TH1D> NStripBCh0Arm0Pl3;
  std::auto_ptr<TH1D> NStripBCh0Arm0Pl4;

  std::auto_ptr<TH1D> NStripBCh0Arm1Pl0;
  std::auto_ptr<TH1D> NStripBCh0Arm1Pl1;
  std::auto_ptr<TH1D> NStripBCh0Arm1Pl2;
  std::auto_ptr<TH1D> NStripBCh0Arm1Pl3;
  std::auto_ptr<TH1D> NStripBCh0Arm1Pl4;


  std::auto_ptr<TH1D> NWireCh1Arm0Pl0;
  std::auto_ptr<TH1D> NWireCh1Arm0Pl1;
  std::auto_ptr<TH1D> NWireCh1Arm0Pl2;
  std::auto_ptr<TH1D> NWireCh1Arm0Pl3;
  std::auto_ptr<TH1D> NWireCh1Arm0Pl4;

  std::auto_ptr<TH1D> NWireCh1Arm1Pl0;
  std::auto_ptr<TH1D> NWireCh1Arm1Pl1;
  std::auto_ptr<TH1D> NWireCh1Arm1Pl2;
  std::auto_ptr<TH1D> NWireCh1Arm1Pl3;
  std::auto_ptr<TH1D> NWireCh1Arm1Pl4;

  std::auto_ptr<TH1D> NStripACh1Arm0Pl0;
  std::auto_ptr<TH1D> NStripACh1Arm0Pl1;
  std::auto_ptr<TH1D> NStripACh1Arm0Pl2;
  std::auto_ptr<TH1D> NStripACh1Arm0Pl3;
  std::auto_ptr<TH1D> NStripACh1Arm0Pl4;

  std::auto_ptr<TH1D> NStripACh1Arm1Pl0;
  std::auto_ptr<TH1D> NStripACh1Arm1Pl1;
  std::auto_ptr<TH1D> NStripACh1Arm1Pl2;
  std::auto_ptr<TH1D> NStripACh1Arm1Pl3;
  std::auto_ptr<TH1D> NStripACh1Arm1Pl4;

  std::auto_ptr<TH1D> NStripBCh1Arm0Pl0;
  std::auto_ptr<TH1D> NStripBCh1Arm0Pl1;
  std::auto_ptr<TH1D> NStripBCh1Arm0Pl2;
  std::auto_ptr<TH1D> NStripBCh1Arm0Pl3;
  std::auto_ptr<TH1D> NStripBCh1Arm0Pl4;

  std::auto_ptr<TH1D> NStripBCh1Arm1Pl0;
  std::auto_ptr<TH1D> NStripBCh1Arm1Pl1;
  std::auto_ptr<TH1D> NStripBCh1Arm1Pl2;
  std::auto_ptr<TH1D> NStripBCh1Arm1Pl3;
  std::auto_ptr<TH1D> NStripBCh1Arm1Pl4;


  std::auto_ptr<TH1D> NWireCh2Arm0Pl0;
  std::auto_ptr<TH1D> NWireCh2Arm0Pl1;
  std::auto_ptr<TH1D> NWireCh2Arm0Pl2;
  std::auto_ptr<TH1D> NWireCh2Arm0Pl3;
  std::auto_ptr<TH1D> NWireCh2Arm0Pl4;

  std::auto_ptr<TH1D> NWireCh2Arm1Pl0;
  std::auto_ptr<TH1D> NWireCh2Arm1Pl1;
  std::auto_ptr<TH1D> NWireCh2Arm1Pl2;
  std::auto_ptr<TH1D> NWireCh2Arm1Pl3;
  std::auto_ptr<TH1D> NWireCh2Arm1Pl4;

  std::auto_ptr<TH1D> NStripACh2Arm0Pl0;
  std::auto_ptr<TH1D> NStripACh2Arm0Pl1;
  std::auto_ptr<TH1D> NStripACh2Arm0Pl2;
  std::auto_ptr<TH1D> NStripACh2Arm0Pl3;
  std::auto_ptr<TH1D> NStripACh2Arm0Pl4;

  std::auto_ptr<TH1D> NStripACh2Arm1Pl0;
  std::auto_ptr<TH1D> NStripACh2Arm1Pl1;
  std::auto_ptr<TH1D> NStripACh2Arm1Pl2;
  std::auto_ptr<TH1D> NStripACh2Arm1Pl3;
  std::auto_ptr<TH1D> NStripACh2Arm1Pl4;

  std::auto_ptr<TH1D> NStripBCh2Arm0Pl0;
  std::auto_ptr<TH1D> NStripBCh2Arm0Pl1;
  std::auto_ptr<TH1D> NStripBCh2Arm0Pl2;
  std::auto_ptr<TH1D> NStripBCh2Arm0Pl3;
  std::auto_ptr<TH1D> NStripBCh2Arm0Pl4;

  std::auto_ptr<TH1D> NStripBCh2Arm1Pl0;
  std::auto_ptr<TH1D> NStripBCh2Arm1Pl1;
  std::auto_ptr<TH1D> NStripBCh2Arm1Pl2;
  std::auto_ptr<TH1D> NStripBCh2Arm1Pl3;
  std::auto_ptr<TH1D> NStripBCh2Arm1Pl4;


  std::auto_ptr<TH1D> NWireCh3Arm0Pl0;
  std::auto_ptr<TH1D> NWireCh3Arm0Pl1;
  std::auto_ptr<TH1D> NWireCh3Arm0Pl2;
  std::auto_ptr<TH1D> NWireCh3Arm0Pl3;
  std::auto_ptr<TH1D> NWireCh3Arm0Pl4;

  std::auto_ptr<TH1D> NWireCh3Arm1Pl0;
  std::auto_ptr<TH1D> NWireCh3Arm1Pl1;
  std::auto_ptr<TH1D> NWireCh3Arm1Pl2;
  std::auto_ptr<TH1D> NWireCh3Arm1Pl3;
  std::auto_ptr<TH1D> NWireCh3Arm1Pl4;

  std::auto_ptr<TH1D> NStripACh3Arm0Pl0;
  std::auto_ptr<TH1D> NStripACh3Arm0Pl1;
  std::auto_ptr<TH1D> NStripACh3Arm0Pl2;
  std::auto_ptr<TH1D> NStripACh3Arm0Pl3;
  std::auto_ptr<TH1D> NStripACh3Arm0Pl4;

  std::auto_ptr<TH1D> NStripACh3Arm1Pl0;
  std::auto_ptr<TH1D> NStripACh3Arm1Pl1;
  std::auto_ptr<TH1D> NStripACh3Arm1Pl2;
  std::auto_ptr<TH1D> NStripACh3Arm1Pl3;
  std::auto_ptr<TH1D> NStripACh3Arm1Pl4;

  std::auto_ptr<TH1D> NStripBCh3Arm0Pl0;
  std::auto_ptr<TH1D> NStripBCh3Arm0Pl1;
  std::auto_ptr<TH1D> NStripBCh3Arm0Pl2;
  std::auto_ptr<TH1D> NStripBCh3Arm0Pl3;
  std::auto_ptr<TH1D> NStripBCh3Arm0Pl4;

  std::auto_ptr<TH1D> NStripBCh3Arm1Pl0;
  std::auto_ptr<TH1D> NStripBCh3Arm1Pl1;
  std::auto_ptr<TH1D> NStripBCh3Arm1Pl2;
  std::auto_ptr<TH1D> NStripBCh3Arm1Pl3;
  std::auto_ptr<TH1D> NStripBCh3Arm1Pl4;


  std::auto_ptr<TH1D> NWireCh4Arm0Pl0;
  std::auto_ptr<TH1D> NWireCh4Arm0Pl1;
  std::auto_ptr<TH1D> NWireCh4Arm0Pl2;
  std::auto_ptr<TH1D> NWireCh4Arm0Pl3;
  std::auto_ptr<TH1D> NWireCh4Arm0Pl4;

  std::auto_ptr<TH1D> NWireCh4Arm1Pl0;
  std::auto_ptr<TH1D> NWireCh4Arm1Pl1;
  std::auto_ptr<TH1D> NWireCh4Arm1Pl2;
  std::auto_ptr<TH1D> NWireCh4Arm1Pl3;
  std::auto_ptr<TH1D> NWireCh4Arm1Pl4;

  std::auto_ptr<TH1D> NStripACh4Arm0Pl0;
  std::auto_ptr<TH1D> NStripACh4Arm0Pl1;
  std::auto_ptr<TH1D> NStripACh4Arm0Pl2;
  std::auto_ptr<TH1D> NStripACh4Arm0Pl3;
  std::auto_ptr<TH1D> NStripACh4Arm0Pl4;

  std::auto_ptr<TH1D> NStripACh4Arm1Pl0;
  std::auto_ptr<TH1D> NStripACh4Arm1Pl1;
  std::auto_ptr<TH1D> NStripACh4Arm1Pl2;
  std::auto_ptr<TH1D> NStripACh4Arm1Pl3;
  std::auto_ptr<TH1D> NStripACh4Arm1Pl4;

  std::auto_ptr<TH1D> NStripBCh4Arm0Pl0;
  std::auto_ptr<TH1D> NStripBCh4Arm0Pl1;
  std::auto_ptr<TH1D> NStripBCh4Arm0Pl2;
  std::auto_ptr<TH1D> NStripBCh4Arm0Pl3;
  std::auto_ptr<TH1D> NStripBCh4Arm0Pl4;

  std::auto_ptr<TH1D> NStripBCh4Arm1Pl0;
  std::auto_ptr<TH1D> NStripBCh4Arm1Pl1;
  std::auto_ptr<TH1D> NStripBCh4Arm1Pl2;
  std::auto_ptr<TH1D> NStripBCh4Arm1Pl3;
  std::auto_ptr<TH1D> NStripBCh4Arm1Pl4;


  std::auto_ptr<TH1D> NWireCh5Arm0Pl0;
  std::auto_ptr<TH1D> NWireCh5Arm0Pl1;
  std::auto_ptr<TH1D> NWireCh5Arm0Pl2;
  std::auto_ptr<TH1D> NWireCh5Arm0Pl3;
  std::auto_ptr<TH1D> NWireCh5Arm0Pl4;

  std::auto_ptr<TH1D> NWireCh5Arm1Pl0;
  std::auto_ptr<TH1D> NWireCh5Arm1Pl1;
  std::auto_ptr<TH1D> NWireCh5Arm1Pl2;
  std::auto_ptr<TH1D> NWireCh5Arm1Pl3;
  std::auto_ptr<TH1D> NWireCh5Arm1Pl4;

  std::auto_ptr<TH1D> NStripACh5Arm0Pl0;
  std::auto_ptr<TH1D> NStripACh5Arm0Pl1;
  std::auto_ptr<TH1D> NStripACh5Arm0Pl2;
  std::auto_ptr<TH1D> NStripACh5Arm0Pl3;
  std::auto_ptr<TH1D> NStripACh5Arm0Pl4;

  std::auto_ptr<TH1D> NStripACh5Arm1Pl0;
  std::auto_ptr<TH1D> NStripACh5Arm1Pl1;
  std::auto_ptr<TH1D> NStripACh5Arm1Pl2;
  std::auto_ptr<TH1D> NStripACh5Arm1Pl3;
  std::auto_ptr<TH1D> NStripACh5Arm1Pl4;

  std::auto_ptr<TH1D> NStripBCh5Arm0Pl0;
  std::auto_ptr<TH1D> NStripBCh5Arm0Pl1;
  std::auto_ptr<TH1D> NStripBCh5Arm0Pl2;
  std::auto_ptr<TH1D> NStripBCh5Arm0Pl3;
  std::auto_ptr<TH1D> NStripBCh5Arm0Pl4;

  std::auto_ptr<TH1D> NStripBCh5Arm1Pl0;
  std::auto_ptr<TH1D> NStripBCh5Arm1Pl1;
  std::auto_ptr<TH1D> NStripBCh5Arm1Pl2;
  std::auto_ptr<TH1D> NStripBCh5Arm1Pl3;
  std::auto_ptr<TH1D> NStripBCh5Arm1Pl4;




  std::auto_ptr<TH1D> RecoX_Arm1Pl4 ; // = new  TH1D("Reco X Arm 1 Pl 4", "Reco X Arm 1 Pl 4", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoX_Arm1Pl3 ; // = new  TH1D("Reco X Arm 1 Pl 3", "Reco X Arm 1 Pl 3", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoX_Arm1Pl2 ; // = new  TH1D("Reco X Arm 1 Pl 2", "Reco X Arm 1 Pl 2", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoX_Arm1Pl1 ; // = new  TH1D("Reco X Arm 1 Pl 1", "Reco X Arm 1 Pl 1", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoX_Arm1Pl0 ; // = new  TH1D("Reco X Arm 1 Pl 0", "Reco X Arm 1 Pl 0", 300, -2000.,2000.);

  std::auto_ptr<TH1D> RecoX_Arm0Pl4 ; // = new  TH1D("Reco X arm 0 Pl 4", "Reco X arm 0 Pl 4", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoX_Arm0Pl3 ; // = new  TH1D("Reco X arm 0 Pl 3", "Reco X arm 0 Pl 3", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoX_Arm0Pl2 ; // = new  TH1D("Reco X arm 0 Pl 2", "Reco X arm 0 Pl 2", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoX_Arm0Pl1 ; // = new  TH1D("Reco X arm 0 Pl 1", "Reco X arm 0 Pl 1", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoX_Arm0Pl0 ; // = new  TH1D("Reco X arm 0 Pl 0", "Reco X arm 0 Pl 0", 300, -2000.,2000.);


  std::auto_ptr<TH1D> RecoY_Arm1Pl4 ; // = new  TH1D("Reco X Arm 1 Pl 4", "Reco X Arm 1 Pl 4", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoY_Arm1Pl3 ; // = new  TH1D("Reco X Arm 1 Pl 3", "Reco X Arm 1 Pl 3", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoY_Arm1Pl2 ; // = new  TH1D("Reco X Arm 1 Pl 2", "Reco X Arm 1 Pl 2", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoY_Arm1Pl1 ; // = new  TH1D("Reco X Arm 1 Pl 1", "Reco X Arm 1 Pl 1", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoY_Arm1Pl0 ; // = new  TH1D("Reco X Arm 1 Pl 0", "Reco X Arm 1 Pl 0", 300, -2000.,2000.);

  std::auto_ptr<TH1D> RecoY_Arm0Pl4 ; // = new  TH1D("Reco X arm 0 Pl 4", "Reco X arm 0 Pl 4", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoY_Arm0Pl3 ; // = new  TH1D("Reco X arm 0 Pl 3", "Reco X arm 0 Pl 3", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoY_Arm0Pl2 ; // = new  TH1D("Reco X arm 0 Pl 2", "Reco X arm 0 Pl 2", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoY_Arm0Pl1 ; // = new  TH1D("Reco X arm 0 Pl 1", "Reco X arm 0 Pl 1", 300, -2000.,2000.);
  std::auto_ptr<TH1D> RecoY_Arm0Pl0 ; // = new  TH1D("Reco X arm 0 Pl 0", "Reco X arm 0 Pl 0", 300, -2000.,2000.);







  std::auto_ptr<TH1D> DeltaX_Arm1Pl4 ; // = new  TH1D("Delta X Arm 1 Pl 4", "Delta X Arm 1 Pl 4", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaX_Arm1Pl3 ; // = new  TH1D("Delta X Arm 1 Pl 3", "Delta X Arm 1 Pl 3", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaX_Arm1Pl2 ; // = new  TH1D("Delta X Arm 1 Pl 2", "Delta X Arm 1 Pl 2", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaX_Arm1Pl1 ; // = new  TH1D("Delta X Arm 1 Pl 1", "Delta X Arm 1 Pl 1", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaX_Arm1Pl0 ; // = new  TH1D("Delta X Arm 1 Pl 0", "Delta X Arm 1 Pl 0", 300, -2000.,2000.);

  std::auto_ptr<TH1D> DeltaX_Arm0Pl4 ; // = new  TH1D("Delta X arm 0 Pl 4", "Delta X arm 0 Pl 4", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaX_Arm0Pl3 ; // = new  TH1D("Delta X arm 0 Pl 3", "Delta X arm 0 Pl 3", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaX_Arm0Pl2 ; // = new  TH1D("Delta X arm 0 Pl 2", "Delta X arm 0 Pl 2", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaX_Arm0Pl1 ; // = new  TH1D("Delta X arm 0 Pl 1", "Delta X arm 0 Pl 1", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaX_Arm0Pl0 ; // = new  TH1D("Delta X arm 0 Pl 0", "Delta X arm 0 Pl 0", 300, -2000.,2000.);


  std::auto_ptr<TH1D> DeltaY_Arm1Pl4 ; // = new  TH1D("Delta X Arm 1 Pl 4", "Delta X Arm 1 Pl 4", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaY_Arm1Pl3 ; // = new  TH1D("Delta X Arm 1 Pl 3", "Delta X Arm 1 Pl 3", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaY_Arm1Pl2 ; // = new  TH1D("Delta X Arm 1 Pl 2", "Delta X Arm 1 Pl 2", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaY_Arm1Pl1 ; // = new  TH1D("Delta X Arm 1 Pl 1", "Delta X Arm 1 Pl 1", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaY_Arm1Pl0 ; // = new  TH1D("Delta X Arm 1 Pl 0", "Delta X Arm 1 Pl 0", 300, -2000.,2000.);

  std::auto_ptr<TH1D> DeltaY_Arm0Pl4 ; // = new  TH1D("Delta X arm 0 Pl 4", "Delta X arm 0 Pl 4", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaY_Arm0Pl3 ; // = new  TH1D("Delta X arm 0 Pl 3", "Delta X arm 0 Pl 3", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaY_Arm0Pl2 ; // = new  TH1D("Delta X arm 0 Pl 2", "Delta X arm 0 Pl 2", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaY_Arm0Pl1 ; // = new  TH1D("Delta X arm 0 Pl 1", "Delta X arm 0 Pl 1", 300, -2000.,2000.);
  std::auto_ptr<TH1D> DeltaY_Arm0Pl0 ; // = new  TH1D("Delta X arm 0 Pl 0", "Delta X arm 0 Pl 0", 300, -2000.,2000.);



std::auto_ptr<TH1D> RecoLX_Arm0Plane0CSC0;
std::auto_ptr<TH1D> RecoLY_Arm0Plane0CSC0;
std::auto_ptr<TH1D> NReco_Arm0Plane0CSC0;
std::auto_ptr<TH1D> RecoLX_Arm0Plane0CSC1;
std::auto_ptr<TH1D> RecoLY_Arm0Plane0CSC1;
std::auto_ptr<TH1D> NReco_Arm0Plane0CSC1;
std::auto_ptr<TH1D> RecoLX_Arm0Plane0CSC2;
std::auto_ptr<TH1D> RecoLY_Arm0Plane0CSC2;
std::auto_ptr<TH1D> NReco_Arm0Plane0CSC2;
std::auto_ptr<TH1D> RecoLX_Arm0Plane0CSC3;
std::auto_ptr<TH1D> RecoLY_Arm0Plane0CSC3;
std::auto_ptr<TH1D> NReco_Arm0Plane0CSC3;
std::auto_ptr<TH1D> RecoLX_Arm0Plane0CSC4;
std::auto_ptr<TH1D> RecoLY_Arm0Plane0CSC4;
std::auto_ptr<TH1D> NReco_Arm0Plane0CSC4;
std::auto_ptr<TH1D> RecoLX_Arm0Plane0CSC5;
std::auto_ptr<TH1D> RecoLY_Arm0Plane0CSC5;
std::auto_ptr<TH1D> NReco_Arm0Plane0CSC5;
std::auto_ptr<TH1D> RecoLX_Arm0Plane1CSC0;
std::auto_ptr<TH1D> RecoLY_Arm0Plane1CSC0;
std::auto_ptr<TH1D> NReco_Arm0Plane1CSC0;
std::auto_ptr<TH1D> RecoLX_Arm0Plane1CSC1;
std::auto_ptr<TH1D> RecoLY_Arm0Plane1CSC1;
std::auto_ptr<TH1D> NReco_Arm0Plane1CSC1;
std::auto_ptr<TH1D> RecoLX_Arm0Plane1CSC2;
std::auto_ptr<TH1D> RecoLY_Arm0Plane1CSC2;
std::auto_ptr<TH1D> NReco_Arm0Plane1CSC2;
std::auto_ptr<TH1D> RecoLX_Arm0Plane1CSC3;
std::auto_ptr<TH1D> RecoLY_Arm0Plane1CSC3;
std::auto_ptr<TH1D> NReco_Arm0Plane1CSC3;
std::auto_ptr<TH1D> RecoLX_Arm0Plane1CSC4;
std::auto_ptr<TH1D> RecoLY_Arm0Plane1CSC4;
std::auto_ptr<TH1D> NReco_Arm0Plane1CSC4;
std::auto_ptr<TH1D> RecoLX_Arm0Plane1CSC5;
std::auto_ptr<TH1D> RecoLY_Arm0Plane1CSC5;
std::auto_ptr<TH1D> NReco_Arm0Plane1CSC5;
std::auto_ptr<TH1D> RecoLX_Arm0Plane2CSC0;
std::auto_ptr<TH1D> RecoLY_Arm0Plane2CSC0;
std::auto_ptr<TH1D> NReco_Arm0Plane2CSC0;
std::auto_ptr<TH1D> RecoLX_Arm0Plane2CSC1;
std::auto_ptr<TH1D> RecoLY_Arm0Plane2CSC1;
std::auto_ptr<TH1D> NReco_Arm0Plane2CSC1;
std::auto_ptr<TH1D> RecoLX_Arm0Plane2CSC2;
std::auto_ptr<TH1D> RecoLY_Arm0Plane2CSC2;
std::auto_ptr<TH1D> NReco_Arm0Plane2CSC2;
std::auto_ptr<TH1D> RecoLX_Arm0Plane2CSC3;
std::auto_ptr<TH1D> RecoLY_Arm0Plane2CSC3;
std::auto_ptr<TH1D> NReco_Arm0Plane2CSC3;
std::auto_ptr<TH1D> RecoLX_Arm0Plane2CSC4;
std::auto_ptr<TH1D> RecoLY_Arm0Plane2CSC4;
std::auto_ptr<TH1D> NReco_Arm0Plane2CSC4;
std::auto_ptr<TH1D> RecoLX_Arm0Plane2CSC5;
std::auto_ptr<TH1D> RecoLY_Arm0Plane2CSC5;
std::auto_ptr<TH1D> NReco_Arm0Plane2CSC5;
std::auto_ptr<TH1D> RecoLX_Arm0Plane3CSC0;
std::auto_ptr<TH1D> RecoLY_Arm0Plane3CSC0;
std::auto_ptr<TH1D> NReco_Arm0Plane3CSC0;
std::auto_ptr<TH1D> RecoLX_Arm0Plane3CSC1;
std::auto_ptr<TH1D> RecoLY_Arm0Plane3CSC1;
std::auto_ptr<TH1D> NReco_Arm0Plane3CSC1;
std::auto_ptr<TH1D> RecoLX_Arm0Plane3CSC2;
std::auto_ptr<TH1D> RecoLY_Arm0Plane3CSC2;
std::auto_ptr<TH1D> NReco_Arm0Plane3CSC2;
std::auto_ptr<TH1D> RecoLX_Arm0Plane3CSC3;
std::auto_ptr<TH1D> RecoLY_Arm0Plane3CSC3;
std::auto_ptr<TH1D> NReco_Arm0Plane3CSC3;
std::auto_ptr<TH1D> RecoLX_Arm0Plane3CSC4;
std::auto_ptr<TH1D> RecoLY_Arm0Plane3CSC4;
std::auto_ptr<TH1D> NReco_Arm0Plane3CSC4;
std::auto_ptr<TH1D> RecoLX_Arm0Plane3CSC5;
std::auto_ptr<TH1D> RecoLY_Arm0Plane3CSC5;
std::auto_ptr<TH1D> NReco_Arm0Plane3CSC5;
std::auto_ptr<TH1D> RecoLX_Arm0Plane4CSC0;
std::auto_ptr<TH1D> RecoLY_Arm0Plane4CSC0;
std::auto_ptr<TH1D> NReco_Arm0Plane4CSC0;
std::auto_ptr<TH1D> RecoLX_Arm0Plane4CSC1;
std::auto_ptr<TH1D> RecoLY_Arm0Plane4CSC1;
std::auto_ptr<TH1D> NReco_Arm0Plane4CSC1;
std::auto_ptr<TH1D> RecoLX_Arm0Plane4CSC2;
std::auto_ptr<TH1D> RecoLY_Arm0Plane4CSC2;
std::auto_ptr<TH1D> NReco_Arm0Plane4CSC2;
std::auto_ptr<TH1D> RecoLX_Arm0Plane4CSC3;
std::auto_ptr<TH1D> RecoLY_Arm0Plane4CSC3;
std::auto_ptr<TH1D> NReco_Arm0Plane4CSC3;
std::auto_ptr<TH1D> RecoLX_Arm0Plane4CSC4;
std::auto_ptr<TH1D> RecoLY_Arm0Plane4CSC4;
std::auto_ptr<TH1D> NReco_Arm0Plane4CSC4;
std::auto_ptr<TH1D> RecoLX_Arm0Plane4CSC5;
std::auto_ptr<TH1D> RecoLY_Arm0Plane4CSC5;
std::auto_ptr<TH1D> NReco_Arm0Plane4CSC5;
std::auto_ptr<TH1D> RecoLX_Arm1Plane0CSC0;
std::auto_ptr<TH1D> RecoLY_Arm1Plane0CSC0;
std::auto_ptr<TH1D> NReco_Arm1Plane0CSC0;
std::auto_ptr<TH1D> RecoLX_Arm1Plane0CSC1;
std::auto_ptr<TH1D> RecoLY_Arm1Plane0CSC1;
std::auto_ptr<TH1D> NReco_Arm1Plane0CSC1;
std::auto_ptr<TH1D> RecoLX_Arm1Plane0CSC2;
std::auto_ptr<TH1D> RecoLY_Arm1Plane0CSC2;
std::auto_ptr<TH1D> NReco_Arm1Plane0CSC2;
std::auto_ptr<TH1D> RecoLX_Arm1Plane0CSC3;
std::auto_ptr<TH1D> RecoLY_Arm1Plane0CSC3;
std::auto_ptr<TH1D> NReco_Arm1Plane0CSC3;
std::auto_ptr<TH1D> RecoLX_Arm1Plane0CSC4;
std::auto_ptr<TH1D> RecoLY_Arm1Plane0CSC4;
std::auto_ptr<TH1D> NReco_Arm1Plane0CSC4;
std::auto_ptr<TH1D> RecoLX_Arm1Plane0CSC5;
std::auto_ptr<TH1D> RecoLY_Arm1Plane0CSC5;
std::auto_ptr<TH1D> NReco_Arm1Plane0CSC5;
std::auto_ptr<TH1D> RecoLX_Arm1Plane1CSC0;
std::auto_ptr<TH1D> RecoLY_Arm1Plane1CSC0;
std::auto_ptr<TH1D> NReco_Arm1Plane1CSC0;
std::auto_ptr<TH1D> RecoLX_Arm1Plane1CSC1;
std::auto_ptr<TH1D> RecoLY_Arm1Plane1CSC1;
std::auto_ptr<TH1D> NReco_Arm1Plane1CSC1;
std::auto_ptr<TH1D> RecoLX_Arm1Plane1CSC2;
std::auto_ptr<TH1D> RecoLY_Arm1Plane1CSC2;
std::auto_ptr<TH1D> NReco_Arm1Plane1CSC2;
std::auto_ptr<TH1D> RecoLX_Arm1Plane1CSC3;
std::auto_ptr<TH1D> RecoLY_Arm1Plane1CSC3;
std::auto_ptr<TH1D> NReco_Arm1Plane1CSC3;
std::auto_ptr<TH1D> RecoLX_Arm1Plane1CSC4;
std::auto_ptr<TH1D> RecoLY_Arm1Plane1CSC4;
std::auto_ptr<TH1D> NReco_Arm1Plane1CSC4;
std::auto_ptr<TH1D> RecoLX_Arm1Plane1CSC5;
std::auto_ptr<TH1D> RecoLY_Arm1Plane1CSC5;
std::auto_ptr<TH1D> NReco_Arm1Plane1CSC5;
std::auto_ptr<TH1D> RecoLX_Arm1Plane2CSC0;
std::auto_ptr<TH1D> RecoLY_Arm1Plane2CSC0;
std::auto_ptr<TH1D> NReco_Arm1Plane2CSC0;
std::auto_ptr<TH1D> RecoLX_Arm1Plane2CSC1;
std::auto_ptr<TH1D> RecoLY_Arm1Plane2CSC1;
std::auto_ptr<TH1D> NReco_Arm1Plane2CSC1;
std::auto_ptr<TH1D> RecoLX_Arm1Plane2CSC2;
std::auto_ptr<TH1D> RecoLY_Arm1Plane2CSC2;
std::auto_ptr<TH1D> NReco_Arm1Plane2CSC2;
std::auto_ptr<TH1D> RecoLX_Arm1Plane2CSC3;
std::auto_ptr<TH1D> RecoLY_Arm1Plane2CSC3;
std::auto_ptr<TH1D> NReco_Arm1Plane2CSC3;
std::auto_ptr<TH1D> RecoLX_Arm1Plane2CSC4;
std::auto_ptr<TH1D> RecoLY_Arm1Plane2CSC4;
std::auto_ptr<TH1D> NReco_Arm1Plane2CSC4;
std::auto_ptr<TH1D> RecoLX_Arm1Plane2CSC5;
std::auto_ptr<TH1D> RecoLY_Arm1Plane2CSC5;
std::auto_ptr<TH1D> NReco_Arm1Plane2CSC5;
std::auto_ptr<TH1D> RecoLX_Arm1Plane3CSC0;
std::auto_ptr<TH1D> RecoLY_Arm1Plane3CSC0;
std::auto_ptr<TH1D> NReco_Arm1Plane3CSC0;
std::auto_ptr<TH1D> RecoLX_Arm1Plane3CSC1;
std::auto_ptr<TH1D> RecoLY_Arm1Plane3CSC1;
std::auto_ptr<TH1D> NReco_Arm1Plane3CSC1;
std::auto_ptr<TH1D> RecoLX_Arm1Plane3CSC2;
std::auto_ptr<TH1D> RecoLY_Arm1Plane3CSC2;
std::auto_ptr<TH1D> NReco_Arm1Plane3CSC2;
std::auto_ptr<TH1D> RecoLX_Arm1Plane3CSC3;
std::auto_ptr<TH1D> RecoLY_Arm1Plane3CSC3;
std::auto_ptr<TH1D> NReco_Arm1Plane3CSC3;
std::auto_ptr<TH1D> RecoLX_Arm1Plane3CSC4;
std::auto_ptr<TH1D> RecoLY_Arm1Plane3CSC4;
std::auto_ptr<TH1D> NReco_Arm1Plane3CSC4;
std::auto_ptr<TH1D> RecoLX_Arm1Plane3CSC5;
std::auto_ptr<TH1D> RecoLY_Arm1Plane3CSC5;
std::auto_ptr<TH1D> NReco_Arm1Plane3CSC5;
std::auto_ptr<TH1D> RecoLX_Arm1Plane4CSC0;
std::auto_ptr<TH1D> RecoLY_Arm1Plane4CSC0;
std::auto_ptr<TH1D> NReco_Arm1Plane4CSC0;
std::auto_ptr<TH1D> RecoLX_Arm1Plane4CSC1;
std::auto_ptr<TH1D> RecoLY_Arm1Plane4CSC1;
std::auto_ptr<TH1D> NReco_Arm1Plane4CSC1;
std::auto_ptr<TH1D> RecoLX_Arm1Plane4CSC2;
std::auto_ptr<TH1D> RecoLY_Arm1Plane4CSC2;
std::auto_ptr<TH1D> NReco_Arm1Plane4CSC2;
std::auto_ptr<TH1D> RecoLX_Arm1Plane4CSC3;
std::auto_ptr<TH1D> RecoLY_Arm1Plane4CSC3;
std::auto_ptr<TH1D> NReco_Arm1Plane4CSC3;
std::auto_ptr<TH1D> RecoLX_Arm1Plane4CSC4;
std::auto_ptr<TH1D> RecoLY_Arm1Plane4CSC4;
std::auto_ptr<TH1D> NReco_Arm1Plane4CSC4;
std::auto_ptr<TH1D> RecoLX_Arm1Plane4CSC5;
std::auto_ptr<TH1D> RecoLY_Arm1Plane4CSC5;
std::auto_ptr<TH1D> NReco_Arm1Plane4CSC5;








  std::auto_ptr<TH1D>  hAllEtaRec;// = new TH1D("AllEtaRec","AllEtaRec",480,-15,15);
  std::auto_ptr<TH1D> hAllPhiRec;// = new TH1D("AllPhiRec","AllPhiRec",50,-4,4);

  std::auto_ptr<TH1D>  hChiSquaredOverN;// = new TH1D("ChiSquaredOverN","ChiSquaredOverN",100,0,30);
  std::auto_ptr<TH1D>  hChiSquaredProb;
  std::auto_ptr<TH1D>  hChiSquaredXOverN;// = new TH1D("ChiSquaredOverN","ChiSquaredOverN",100,0,30);
  std::auto_ptr<TH1D>  hChiSquaredXProb;
  std::auto_ptr<TH1D>  hChiSquaredYOverN;// = new TH1D("ChiSquaredOverN","ChiSquaredOverN",100,0,30);
  std::auto_ptr<TH1D>  hChiSquaredYProb;
  std::auto_ptr<T1Geometry> layer;

  std::auto_ptr<TH1D> lEta;
  //  std::auto_ptr<TH1D> lEtaH;

  std::auto_ptr<TH1D> ilPT;
  //  std::auto_ptr<TH1D> ilPTH;

  // dai simhits
  std::auto_ptr<TH1D> ilMomento;
  std::auto_ptr<TH1D> ilMomentoLog;
  std::auto_ptr<TH1D> ilTipo;
  //  std::auto_ptr<TH1D> ilProcesso;
  std::auto_ptr<TH1D> lEnergia;
  std::auto_ptr<TH1D> ilDetector;

  TFile* theFile;

};
#endif // TotemT1T2ValidationT1ValidationOnlyDataT1ValidationOnlyData
