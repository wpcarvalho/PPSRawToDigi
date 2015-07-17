/*
 * modified by Marcin Borratynski (mborratynski@gmail.com)
 */
#ifndef _T1DIGIPRODUCER_H
#define _T1DIGIPRODUCER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/T1DigiWire/interface/T1DigiWire.h"
#include "DataFormats/T1DigiWire/interface/T1DigiWireCollection.h"
#include "DataFormats/T1DigiSantiard/interface/T1DigiSantiard.h"
#include "DataFormats/T1DigiSantiard/interface/T1DigiSantiardCollection.h"
#include "DataFormats/T1DigiVfat/interface/T1DigiVfat.h"
#include "DataFormats/T1DigiVfat/interface/T1DigiVfatCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/Provenance/interface/BranchDescription.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "SimTotem/T1Digitizer/interface/T1Digitizer.h"
#include "SimTotem/T1Digitizer/interface/T1DeadChannelManager.h"
#include "TotemCondFormats/DAQInformation/interface/AnalysisMask.h"
#include "TotemCondFormats/DataRecord/interface/TotemDAQMappingRecord.h"



#include <fstream>
#include <string>
#include <cctype>
#include <map>
//#define _DEBUG_
class T1DigiProducer : public edm::EDProducer
{
 public:

  // The following is not yet used, but will be the primary
  // constructor when the parameter set system is available.
  //
  explicit T1DigiProducer(const edm::ParameterSet& params);
  virtual ~T1DigiProducer();

  /**Produces the EDM products,*/
  virtual void produce(edm::Event& event, const edm::EventSetup& eventSetup);
  virtual void endJob() ;
  virtual void beginRun(edm::Run&, edm::EventSetup const&);
 private:
 
  int evento;
  int theVerbosity;
  T1Digitizer  *_digitizer;

  std::string _electronics;

  T1DeadChannelManager _deadChannelManager;
  /**
   * this variable indicates whether we take into account dead channels or simulate as if all
   * channels work ok (by default we do not simulate dead channels)
   */
  bool simulateDeadChannels;

};

#endif 
