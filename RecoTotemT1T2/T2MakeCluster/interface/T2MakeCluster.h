#ifndef T2MakeCluster_h
#define T2MakeCluster_h

/** 
 * Class T2MakeCluster for clusterization 
 * 
 * Author: Mirko Berretti / University of Siena
 * Email:  mirko.berretti@gmail.com
 * Date:   2007-12-09
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoTotemT1T2/T2MakeCluster/interface/T2DetClustReconst.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CLHEP/Random/RandomEngine.h"
#include <TH2.h>
#include <TH1.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TFile.h>
#include <boost/utility.hpp>

class T2MakeCluster : public edm::EDProducer {
 public:

  explicit
    T2MakeCluster (const edm::ParameterSet & paraSet);

  virtual ~T2MakeCluster (){}

  virtual void produce (edm::Event & ev, const edm::EventSetup & evSet);
  void GetEffiSector(double Rcl, double Phicl, uint32_t thedet, int &rCell, int &phiCell, int &symbId, int numRsectEffi, int numPhisectEffi);

  virtual void beginJob() ;

  CLHEP::HepRandomEngine* rndEnginePR;

  
  //std::auto_ptr<TH2D> EffiStrip[40];
  // std::auto_ptr<TH2D> EffiPad[40];
  //boost::shared_ptr<TH2D> EffiPad[40];
  //boost::shared_ptr<TH2D> EffiStrip[40];
  std::vector<TH2D> EffiPad;//[40];
  std::vector<TH2D> EffiStrip;//[40];

  boost::shared_ptr<TFile> GeometryEffiFile;
 private:
  std::vector<unsigned int> maskvect;
  T2DetClustReconst theT2Clusterizer;
  bool TakeCleanEventOnly;

  int numrbin;
  int numphibin;

  double Projectionthreshold;
  int BlobMinSize;
  bool SimuClusterEfficiency;
  string EffiGeoRootFileName;
  edm::InputTag T2PadDigiCollectionLabel;
  edm::InputTag T2StripDigiCollectionLabel;
};

#endif
