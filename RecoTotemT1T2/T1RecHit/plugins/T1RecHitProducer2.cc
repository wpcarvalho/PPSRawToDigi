// -*- C++ -*-
//
// Package:    T1RecHit
// Class:      T1RecHitProducer2
// 
/**\class T1RecHitProducer2 T1RecHitProducer2.cc RecoTotemT1T2/T1RecHit/src/T1RecHitProducer2.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Enrico Robutti
//         Created:  Wed Mar  7 15:37:12 CET 2012
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <list>
#include <vector>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include "Geometry/TotemGeometry/interface/T1GeomPars.h"
#include "DataFormats/T1DigiWire/interface/T1DigiWireCollection.h"
#include "DataFormats/T1DigiVfat/interface/T1DigiVfatCollection.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"

using namespace std;
using namespace edm;

//
// class declaration
//

class T1RecHitProducer2 : public edm::EDProducer {
public:
  explicit T1RecHitProducer2(const edm::ParameterSet&);
  ~T1RecHitProducer2();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  virtual void endRun(edm::Run&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  
  // ----------member data ---------------------------
  string _hitProductLabel;
  string _t1DigiModuleLabel;
  string _t1DigiProductLabel;
  T1Geometry _geometry;
  T1GeomPars _geomPars;
  unsigned int _verbosity;
  unsigned int _printProgressFrequency;
  float _tripleTolerance;
};
  //
  // constants, enums and typedefs
  //
  
  
  //
  // static data member definitions
  //
  
  //
  // constructors and destructor
  //

T1RecHitProducer2::T1RecHitProducer2(const edm::ParameterSet& iConfig) :
  _hitProductLabel(iConfig.getUntrackedParameter<string>("hitProductLabel", "3c-clusters")),
  _t1DigiModuleLabel(iConfig.getParameter<string>("t1DigiModuleLabel")),
  _t1DigiProductLabel(iConfig.getParameter<string>("t1DigiProductLabel")),
  _geometry(iConfig.getParameter<string>("geometryFile")),
  _geomPars(iConfig.getParameter<string>("geomParsFile")),
  _verbosity(iConfig.getUntrackedParameter<unsigned int>("verbosity", 0)),
  _printProgressFrequency(iConfig.getUntrackedParameter<unsigned int>("printProgressFrequency", 1000)),
  _tripleTolerance(iConfig.getParameter<double>("tripleTolerance"))
{
   // Register products
  produces<T1RecHit2DCollection>(_hitProductLabel);

}


T1RecHitProducer2::~T1RecHitProducer2()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
T1RecHitProducer2::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Create pointer to the rechits collection
  auto_ptr<T1RecHit2DCollection> recHitCollection(new T1RecHit2DCollection());

  // Get digis from event  
  Handle<T1DigiWireCollection> wireCollection;
  iEvent.getByLabel(_t1DigiModuleLabel, _t1DigiProductLabel, wireCollection);
  Handle<T1DigiVfatCollection> stripCollection;
  iEvent.getByLabel(_t1DigiModuleLabel, _t1DigiProductLabel, stripCollection);

  // Loop over digis
  unsigned int evtNum = iEvent.id().event();
  if (_verbosity > 0 && evtNum%_printProgressFrequency == 0) {
    cout << "Event n. " << evtNum << endl;
    LogInfo("T1RecHitProducer2") << "Event n. " << evtNum;
  }
  const float deltaV = 0.5 + _tripleTolerance/5.;
  for (T1DigiWireCollection::DigiRangeIterator itWC = wireCollection->begin(); itWC != wireCollection->end(); itWC++) {
    list<int> wHits, uHits, vHits;
    bool wDone = false;
    bool uDone = false;
    bool vDone = false;
    const T1DetId& wCsc = (*itWC).first;
    int arm = wCsc.Arm();
    int layer = wCsc.Plane();
    int sextant = wCsc.CSC();
    for (T1DigiVfatCollection::DigiRangeIterator itSC = stripCollection->begin(); itSC != stripCollection->end() && !(uDone && vDone); itSC++) {
      const T1DetId& sCsc = (*itSC).first;
      if (sCsc.Arm() == arm && sCsc.Plane() == layer && sCsc.CSC() == sextant) {
	// Build sorted digi lists
	if (!wDone) {
	  const T1DigiWireCollection::Range& wRange = (*itWC).second;
	  for (T1DigiWireCollection::const_iterator itW = wRange.first; itW != wRange.second; itW++)
	    wHits.push_back(itW->wire());
	  wHits.sort();
	  wDone = true;
	}
	const T1DigiVfatCollection::Range& sRange = (*itSC).second;
	if (sCsc.Layer() == 1) {
	  for (T1DigiVfatCollection::const_iterator itS = sRange.first; itS != sRange.second; itS++)
	    uHits.push_back(itS->strip());
	  uHits.sort();
	  uDone = true;
	} else if (sCsc.Layer() == 2) {
	  for (T1DigiVfatCollection::const_iterator itS = sRange.first; itS != sRange.second; itS++)
	    vHits.push_back(itS->strip());
	  vHits.sort();
	  vDone = true;
	}
      }
      if (_verbosity > 2)
	cout << "Found " << wHits.size() << " wires, " << uHits.size() << " 'u' strips and " << vHits.size() << " 'v' strips." << endl;	
      // Find strip crossings and make clusters
      if (uDone && vDone && wDone) {
	list<pair<unsigned int, pair<float, float> > > clusterList;
	// Loop on hit wires
	for (list<int>::iterator itW = wHits.begin(); itW != wHits.end(); itW++) {
	  int previousU = -10;
	  // Scan hit v-strips from last to first
	  list<int>::iterator itV = vHits.end();
	  itV--;
	  // Loop on hit u-strip
	  float prevVTriple = 200.;
	  for (list<int>::iterator itU = uHits.begin(); itU != uHits.end(); itU++) {
	    float vTriple = 3./5.*(*itW - _geomPars.getWireOffset(layer, sextant)) - (*itU);
	    while (*itV > vTriple + deltaV && itV != vHits.begin()) {
	      itV--;
	      // Check for a "leftover" v-strip on previous cluster
	      if ((int)(*itU) > previousU + 1 && fabs(prevVTriple - *itV) < deltaV) {
		clusterList.back().first++;  // increment last cluster size
		clusterList.back().second.first -= 2.5/sqrt(3.);  // shift x-coord of last cluster
	      }
	    }
	    if (*itV > vTriple + deltaV)
	      break;
	    if (fabs(vTriple - *itV) < deltaV) {  // found new triple coincidence
	      if ((int)(*itU) > previousU + 1)  // new cluster
		clusterList.push_back(pair<unsigned int, pair<float, float> >(1, pair<float, float>(_geometry.xOfStripCrossing(sCsc.rawId(), *itU, *itV), _geometry.yOfWire(wCsc.rawId(), *itW))));
	      else
		clusterList.back().first++;  // increment last cluster size
	      clusterList.back().second.first -= 5./sqrt(3.);  // shift x-coord of last cluster
	      previousU = *itU;
	      prevVTriple = vTriple;
	    }
	  }
	}
	// Loop on clusters and build hit2D vector
	OwnVector<T1RecHit2D> vHit2D;
	for (list<pair<unsigned int, pair<float, float> > >::iterator itClust = clusterList.begin(); itClust != clusterList.end(); itClust++) {
	  LocalPoint lp(itClust->second.first, itClust->second.second, 0.);
	  LocalError le(0.85, 0.85, 0.);  // approximate isotropical error
	  T1RecHit2D* hit2D = new T1RecHit2D(wCsc, lp, le, deltaV);
	  vHit2D.push_back(hit2D);
	}
	if (_verbosity > 2)
	  cout << vHit2D.size() << " hits reconstructed." << endl;
	// Insert 2D hits into event
	recHitCollection->put(wCsc, vHit2D.begin(), vHit2D.end());
	break;
      }
    }
  }
  if (_verbosity > 2)
    cout << recHitCollection->size() << " hits written in collection" << endl;
  iEvent.put(recHitCollection, _hitProductLabel);
}

// ------------ method called once each job just before starting event loop  ------------
void 
T1RecHitProducer2::beginJob()
{
  if (_verbosity > 0)
    LogInfo("T1RecHitProducer2") << "Finding 3-cordinate points";
}

// ------------ method called once each job just after ending the event loop  ------------
void 
T1RecHitProducer2::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
T1RecHitProducer2::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
T1RecHitProducer2::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
T1RecHitProducer2::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
T1RecHitProducer2::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
T1RecHitProducer2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(T1RecHitProducer2);
