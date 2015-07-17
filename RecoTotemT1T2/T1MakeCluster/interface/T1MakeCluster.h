#ifndef RecoLocal_T1MakeCluster_h
#define RecoLocal_T1MakeCluster_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"

#include "DataFormats/T1Cluster/interface/T1Cluster.h"
#include "DataFormats/T1Cluster/interface/T1ClusterCollection.h"
#include "TotemCondFormats/DAQInformation/interface/AnalysisMask.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "SimTotem/T1Digitizer/interface/T1DeadChannelManager.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "TotemCondFormats/DAQInformation/interface/DAQMapping.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "TotemCondFormats/DAQInformation/interface/DAQInformationT1.h"
#include "TotemCondFormats/DataRecord/interface/TotemDAQMappingRecord.h"


namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

//class T1RecHitBaseAlgo;

class T1MakeCluster : public edm::EDProducer {
 public:
  /// Constructor
  T1MakeCluster(const edm::ParameterSet&);

  /// Destructor
  virtual ~T1MakeCluster();

  /// The method which produces the rechits
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  void t1DataBeginJob(const edm::EventSetup &es);

  virtual void produce(edm::Event& event, const edm::EventSetup& setup);
  float Sigma(int fstrip, int lstrip, float center, int * st);
  float Center(int fstrip, int lstrip, int * st);
 private:
  edm::InputTag t1DigiSantiardCollectionLabel;
  edm::InputTag t1DigiVfatCollectionLabel;
  string _electronics;
  // The label to be used to retrieve T1 digis from the event
  //  std::string theT1DigiLabel;
  // The reconstruction algorithm
  //  T1RecHitBaseAlgo *theAlgo;
  //   static string theAlgoName;

  float theChiSquare;
  int index;
  T1Geometry * theT1Geometry;
  //std::vector< edm::OwnVector<T1RecHit2D> > * RecHitMatrix;
  float _weight1;
  float _weight2;
  int theVerbosity;

  T1DeadChannelManager _T1deadChannelManager;
  bool _ActivateDeadChannels;

};
#endif
