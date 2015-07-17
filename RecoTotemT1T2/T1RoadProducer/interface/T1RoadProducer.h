/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
*/
#ifndef _T1RoadProducer_h
#define _T1RoadProducer_h


#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"

#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"
#include "DataFormats/T1Road/interface/T1Road.h"
#include <fstream>
#include <sstream>


namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}



class T1RoadProducer : public edm::EDProducer {
 public:
  /// Constructor
  T1RoadProducer(const edm::ParameterSet&);

  float Eta(float,float,float);
  float Phi(float,float);

  /// Destructor
  virtual ~T1RoadProducer();

  /// The method which produces the rechits
  virtual void produce(edm::Event& event, const edm::EventSetup& setup);

 private:

  edm::InputTag t1RecHit2DCollectionLabel;
  float theDeltaEta, theDeltaPhi;
  int theMinHits, theMaxHits;
  int index;
  T1Geometry * theT1Geometry;
  int theVerbosity;
  bool theAlignment;

};
#endif
