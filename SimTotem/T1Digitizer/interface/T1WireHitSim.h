#ifndef T1_WIRE_HIT_SIM_H
#define T1_WIRE_HIT_SIM_H

#include <vector>
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimTotem/T1Digitizer/interface/T1DetectorHit.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "SimTotem/T1Digitizer/interface/T1DriftSim.h"
#include "SimTotem/T1Digitizer/interface/T1WireHitSim.h"
#include "SimTotem/T1Digitizer/interface/T1StripHitSim.h"
#include "SimTotem/T1Digitizer/interface/T1GasCollisions.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"

//class T1DriftSim;
//class T1Geometry;
class T1G3Hit;
//class T1GasCollisions;
//class T1LayerGeometry;


class T1WireHitSim
{
 public:
  explicit T1WireHitSim(T1DriftSim* driftSim);
  ~T1WireHitSim();

  // makes wire hits from the given g3hits
  std::vector<T1DetectorHit> & 
    simulate(uint32_t ,  const edm::PSimHitContainer & simHits, T1Geometry *);
  //  void   simulate();

 private:
  // Helper functions
  std::vector<Local3DPoint> getIonizationClusters(const PSimHit & hit,  T1Geometry *);
  T1DetectorHit driftElectronsToWire();

  // member data
  T1DriftSim*  pDriftSim;
  T1GasCollisions* theGasIonizer;
  std::vector<T1DetectorHit> newWireHits;
};

#endif
