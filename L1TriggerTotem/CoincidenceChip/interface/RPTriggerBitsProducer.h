// -*- C++ -*-
//
// Package:    CoincidenceChip
// Class:      RPTriggerBitsProducer
//
// Original Author:  Leszek Grzanka
// $Id: RPTriggerBitsProducer.h,v 1.1.2.1 2009/07/13 12:23:35 jkaspar Exp $
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

#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Wrapper.h"

#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrameCollection.h"

//
// class declaration
//

class RPTriggerBitsProducer: public edm::EDProducer {
  public:
    explicit RPTriggerBitsProducer(const edm::ParameterSet&);
    ~RPTriggerBitsProducer();

  private:
    virtual void beginJob();

    virtual void produce(edm::Event&, const edm::EventSetup&);

    virtual void endJob();

    // ----------member data ---------------------------
    bool verbose_;
    std::vector<edm::DetSet<RPDetTrigger> > theTriggerVector;    
    edm::InputTag stripDigiLabel;
};
