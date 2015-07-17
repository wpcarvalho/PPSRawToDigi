#ifndef T2_DIGITIZER
#define T2_DIGITIZER

/** 
 * Class T2Digitizer for digitization of T2 GEM detector
 * 
 * Author: Erik Br??cken / University of Helsinki
 *         Mirko Berretti / University of Siena & Pisa INFN
 * Email:  brucken@cc.helsinki.fi
 *         mirko.berretti@gmail.com  
 * Date:   2007-11-26
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/T2Digi/interface/T2StripDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"
#include "Geometry/TotemGeometry/interface/T2Geometry.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimTotem/T2Digitizer/interface/PSimHitMisaligner.h"
#include <DataFormats/T2DigiVfat/interface/T2DigiVfat.h>
#include "DataFormats/T2DigiVfat/interface/T2DigiVfatCollection.h"
#include <boost/utility.hpp>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CLHEP/Random/RandomEngine.h"

class T2StripHitSim;
class T2PadHitSim;
class T2StripElectronicsSim;
class T2PadElectronicsSim;

class T2Digitizer : public boost::noncopyable {
  
 public:
  
  explicit T2Digitizer(const edm::ParameterSet & paraSet,const edm::EventSetup& iSetup,CLHEP::HepRandomEngine* rndEngine);
  
  ~T2Digitizer();
  
  void process(MixCollection<PSimHit> & simHits,
	       T2StripDigiCollection & stripDigis,
	       T2PadDigiCollection & padDigis,
	       const HepMC::GenEvent* hepmcevt,
	       std::auto_ptr<edm::SimVertexContainer> theSimVertices,
	       std::auto_ptr<edm::SimTrackContainer> theSimTracks
	       /*std::map<int, std::map<unsigned int, T2DigiVfat> > & VfatsCollection*/);
  
  std::map<int, std::map<unsigned int, T2DigiVfat> >* StripVFatsdig ;
  std::map<int, std::map<unsigned int, T2DigiVfat> >* PadVFatsdig ;


  bool TakeOnlyPrim;
  bool TakeOnlySec;
  std::vector<int> knowPart; std::vector<int> AllowedDecayToCount;

 private:
  
  T2Geometry* theT2Geometry;
  T2StripHitSim* theStripHitSim;
  T2PadHitSim* thePadHitSim;
  T2StripElectronicsSim* theStripElectronicsSim;
  T2PadElectronicsSim* thePadElectronicsSim;
  unsigned int RawtoSymb(uint32_t thedet);
  boost::shared_ptr<PSimHitMisaligner> thePSimHitMisaligner; 
  std::map<int, int*> stripChargeMap;

  void createStripChargeMap();
  void resetStripChargeMap();
  void deleteStripChargeMap();

  std::map<int, int*> padChargeMap;
  void createPadChargeMap();
  void resetPadChargeMap();
  void deletePadChargeMap();

  /////////////////////////////////////////////////////////////////////////

  bool PrimaryTrackCondition(SimTrack atrk,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt,int &barcodeMother);
  
  bool getTrack(const unsigned int trackId,  const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, SimTrack &track);

  bool getVertex(const  int vertexId,  const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, SimVertex &vertex);

  bool isTrackPrimary(const int trackId, const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const std::auto_ptr<edm::SimVertexContainer>& theSimVertices);
  
  int GetTrkOlderMotherPid(SimTrack aTrack,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt, int &motherbarcode,int &thePidGeneratorOfThisTrk);


  // just for debugging purpose
  void printoutStuff();

  void HitNoisePosition (double &x, double &y,double RNoiseSlope);
};

#endif

