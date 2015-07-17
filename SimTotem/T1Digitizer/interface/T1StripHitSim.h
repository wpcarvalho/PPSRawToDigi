#ifndef T1_STRIP_HIT_SIM_H
#define T1_STRIP_HIT_SIM_H

/** \class CSCStripHitSim
 *
 * Class which builds simulated strip hits from wire
 * hits during digitization of Endcap Muon CSCs.
 *
 * \author Rick Wilkinson
 *
 */

//#include "SimTotem/T1Digitizer/interface/T1GattiFunction.h"
#include <vector>
#include "SimTotem/T1Digitizer/interface/T1DetectorHit.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"

// declarations
//class T1Geometry;
//class T1DetectorHit;

class T1StripHitSim
{
 public:
  // make strip hits from the given wire hits
  std::vector<T1DetectorHit> & simulate( uint32_t myDetId, 
					 const std::vector<T1DetectorHit> & wireHits, T1Geometry *);

  void initChamberSpecs();

  ///  returns the fraction of charge on a strip centered
  ///  a distance of x away from the center of the shower,
  ///  at zero.  Note that the user is responsible for making
  ///  sure the constants have been initialized using the chamber specs.
  double binValue( double x, double stripWidth) const;

 private:
  //  T1GattiFunction theGattiFunction;
  std::vector<T1DetectorHit> newStripHits;

  double k1, k2, k3, h;

  double norm, sqrtk3;

};

#endif
