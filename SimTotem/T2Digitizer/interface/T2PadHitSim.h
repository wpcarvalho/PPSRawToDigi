#ifndef T2_Pad_Hit_Sim_h
#define T2_Pad_Hit_Sim_h

/** 
 * Class simulates hit on a pad in the T2 detector
 *
 * Author: Erik Brücken / University of Helsinki
 * email:  brucken@cc.helsinki.fi
 * Date    2007-11-26
 */

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTotem/T2Digitizer/interface/T2DetectorHit.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimTotem/T2Digitizer/interface/PSimHitMisaligner.h"
#include <vector>
#include "CLHEP/Random/Randomize.h"

#define NUMBER_OF_ROWS   24
#define NUMBER_OF_COLS   65

class T2Geometry;
class T2ROGeometry;

class T2PadHitSim {
  
 public:
  
  explicit T2PadHitSim(const edm::ParameterSet & paraSet);

  ~T2PadHitSim();
  
  std::vector<T2DetectorHit> & simulate(T2Geometry* geom, 
                                        const edm::PSimHitContainer & simHits,
					std::vector<double> & chargeDiv,boost::shared_ptr<PSimHitMisaligner> thePSimHitMisaligner);
 
 

 private:

  T2ROGeometry* theT2ROGeometry;
  
  std::vector<T2DetectorHit> theNewPadHits;

  std::vector<double> diffCoeff_;
  int NUMBER_OF_SIM_STEPS_;
  double eIonE_;
  std::vector<double> gain_;
  double z_max_;
  double z_min_;

  unsigned int RawtoSymb(uint32_t thedet);

  // function to calculate the sigma of the electroncloud
  double calculateSigma(double diffCoeff, double zPosition);


  // function to calculate the the M3 factor for the erf function
  double calculateM3(double sigma);
/* 
  // function to calculate the charge collected on a pad using the erf function
  // It has been shown that the fitted erf function just depends on M3 
  // with a small error. n(x) = 50 ( 1 - erf( M3 * x ) )
  double calculatePadChargeCollPart(double m3, double distance);
*/

 
  // function to calculate the charge collected on a pad using the erf function.
  // In this new version the charge collected is obtained considering the borders of pads.
  // It has been shown that the fitted erf function just depends on M3 
  // with a small error. 
  // n(x) = fabs( 50*(1 - erf( m3 * (fabs(DMin)) )) - 50*( 1 - erf( m3 * (fabs(DMax)) )) )
  // or
  // n(x) = fabs( 100 - 50*(1 - erf( m3 * (fabs(DMin)) )) - 50*( 1 - erf( m3 * (fabs(DMax)) )) )
  // depending on the position of the cloud respect with the readout electrode.
  double calculatePadChargeCollPart(double m3, double DMin , double DMax);

  // function to calculate the charge percentage on all involved pads
  void calculate_pad_charge(T2ROGeometry* geometry, double (*charges)[NUMBER_OF_COLS], 
			    int row, int col, double r, double phi,  double z);
};

#endif
