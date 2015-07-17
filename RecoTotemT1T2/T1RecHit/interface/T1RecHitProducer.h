#ifndef RecoLocal_T1RecHitProducer_h
#define RecoLocal_T1RecHitProducer_h


#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"

//#include "RecoTotemT1T2/T1RecHit/interface/T1Clusterizer.h"
#include "DataFormats/T1Cluster/interface/T1Cluster.h"
#include "DataFormats/T1Cluster/interface/T1ClusterCollection.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

//class T1RecHitBaseAlgo;

class T1RecHitProducer : public edm::EDProducer {
 public:
  /// Constructor
  T1RecHitProducer(const edm::ParameterSet&);

  /// Destructor
  virtual ~T1RecHitProducer();

  /// The method which produces the rechits
  virtual void produce(edm::Event& event, const edm::EventSetup& setup);

 private:

  // The label to be used to retrieve T1 digis from the event
  //  std::string theT1DigiLabel;
  // The reconstruction algorithm
  //  T1RecHitBaseAlgo *theAlgo;
  //   static string theAlgoName;

  edm::InputTag t1ClusterCollectionLabel;
  edm::InputTag t1DigiWireCollectionLabel;
  float theChiSquare;
  int index;
  T1Geometry * theT1Geometry;
  int theVerbosity;
  //std::vector< edm::OwnVector<T1RecHit2D> > * RecHitMatrix;
};
#endif
