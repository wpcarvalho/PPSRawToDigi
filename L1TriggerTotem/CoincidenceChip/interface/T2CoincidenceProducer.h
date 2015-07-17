// -*- C++ -*-
//
// Package:    CoincidenceChip
// Class:      T2CoincidenceProducer
//
/**\class T2CoincidenceProducer T2CoincidenceProducer.cc L1TriggerTotem/CoincidenceChip/src/CoincidenceChip.cc

 Description: <one line class summary>

 Implementation:
 <Notes on implementation>
 */
//
// Original Author:  pts/1
//         Created:  Tue Aug 19 19:02:48 CEST 2008
// $Id: T2CoincidenceProducer.h,v 1.4 2009/04/01 15:27:30 oljemark Exp $
//
//
// Copied and pasted by Fredrik Oljemark, Mon Oct 27 2008
//

// system include files
#include <memory>
#include <iomanip>
#include <ios>
#include <iostream>

#include <cstdlib>
#include <bitset>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "L1TriggerTotem/CoincidenceChip/interface/CoincidenceChip.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"

//needed for the geometry:
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

//include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/T2Digi/interface/T2PadDigi.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"
//include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
//include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/TotemData/interface/TotemDigiCollection.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/TotemL1Trigger/interface/T2TriggerBits.h"
#define numOfMetaPadRows 8
#define numOfMetaPadColumns 13
#define minNumOfPlanes 7

//
// class declaration
//
//typedef std::bitset<16> FOTestOut;

class T2CoincidenceProducer: public edm::EDProducer {
  public:
    explicit T2CoincidenceProducer(const edm::ParameterSet&);
    ~T2CoincidenceProducer();

  private:
    virtual void beginJob();

    virtual void produce(edm::Event&, const edm::EventSetup&);

  virtual void endJob();

  void createMetaPadMap();
  void resetMetaPadMap();
  void deleteMetaPadMap();

  // just for debugging purpose
  //  void printoutStuff();

  // arrays to store number of found paths through detector
  int simplePaths[4];
  //  int stepPaths[4];
  int hh; //debugging

  // pattern algorithms
  void simplePath(int* count);
  //  void stepPath(int* count);


    // ----------member data ---------------------------
       
    CoincidenceChip cchip;

    edm::InputTag src_;
    unsigned int verbose_;
    unsigned int quarter_;
  //bool verbose_;
    std::map<int, int*> metaPadMap;


};
