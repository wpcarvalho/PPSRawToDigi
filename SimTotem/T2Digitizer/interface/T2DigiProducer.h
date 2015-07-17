#ifndef T2DigiProducer_h
#define T2DigiProducer_h

/** 
 * Class T2DigiProducer for digitization of T2 GEM detector
 * 
 * Author: Erik Br??cken / University of Helsinki
 *         Mirko Berretti / University of Siena & Pisa INFN
 * Email:  brucken@cc.helsinki.fi
 *         mirko.berretti@gmail.com  
 * Date:   2007-11-26
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimTotem/T2Digitizer/interface/T2Digitizer.h"


class T2DigiProducer : public edm::EDProducer {
  
 public:
  
  
  explicit T2DigiProducer(const edm::ParameterSet& paraSet);
  
  virtual ~T2DigiProducer() {}
  
  virtual void produce(edm::Event& ev, const edm::EventSetup& evSet);
  virtual void beginJob() ;
  


 private:
 

  //Parameter set for the digitizer
  edm::ParameterSet DigiProdParSet_;
  bool saveDigiVFAT;
  std::string previousModule;
  std::string instanceLabel;
  // T2Digitizer theT2Digitizer;
  CLHEP::HepRandomEngine* rndEnginePR;
  edm::InputTag simTrackContainerLabel;
  edm::InputTag SimVertexContainerLabel;
};

#endif

