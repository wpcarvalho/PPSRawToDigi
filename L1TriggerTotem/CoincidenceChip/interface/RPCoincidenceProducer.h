// -*- C++ -*-
//
// Package:    CoincidenceChip
// Class:      RPCoincidenceProducer
//
// Original Author:  Leszek Grzanka
//         Created:  Tue Aug 19 19:02:48 CEST 2008
// $Id: RPCoincidenceProducer.h,v 1.2.2.2 2009/08/21 13:24:16 lgrzanka Exp $
//
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

#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrameCollection.h"

//needed for the geometry:
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/TotemL1Trigger/interface/RPCCBits.h"

//
// class declaration
//

class RPCoincidenceProducer: public edm::EDProducer {
  public:
    explicit RPCoincidenceProducer(const edm::ParameterSet&);
    ~RPCoincidenceProducer();

  private:
    virtual void beginJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    
    void printBits( std::bitset<80> );

    // ----------member data ---------------------------

    CoincidenceChip cchip;

    unsigned int verbose_;
	std::string productLabelSimu;
	edm::InputTag detTriggerLabel;
};
