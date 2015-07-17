#ifndef T2_Strip_Hit_Sim_h
#define T2_Strip_Hit_Sim_h

/** 
 * Class simulates hit on a strip in the T2 detector. Polynom that calculates
 * the charge distribution has been developed by Eraldo Olivieri
 *
 * Author: Erik Brücken / University of Helsinki
 * email:  brucken@cc.helsinki.fi
 * Date:   2007-11-26
 */

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTotem/T2Digitizer/interface/T2DetectorHit.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimTotem/T2Digitizer/interface/PSimHitMisaligner.h"
#include <vector>
#include "CLHEP/Random/Randomize.h"

class T2Geometry;
class T2ROGeometry;

class T2StripHitSim {
  
 public:
  
  explicit T2StripHitSim(const edm::ParameterSet & paraSet);

  ~T2StripHitSim();
  
  std::vector<T2DetectorHit> & simulate(T2Geometry* geom, 
  					const edm::PSimHitContainer & simHits,
					std::vector<double> & chargeDiv,boost::shared_ptr<PSimHitMisaligner> thePSimHitMisaligner);

 private:
  
  T2ROGeometry* theT2ROGeometry;

  std::vector<T2DetectorHit> theNewStripHits;

 
  std::vector<double> diffCoeff_;
  int NUMBER_OF_SIM_STEPS_;
  double eIonE_;
  std::vector<double> gain_;
  double z_max_;
  double z_min_;
  std::vector<double> StripWidth_;
  // function to calculate the sigma of the electroncloud  
  double calculateSigma(double diffCoeff, double zPosition);
   unsigned int RawtoSymb(uint32_t thedet);
   // function to calculate the the M3 factor for the erf function
  double calculateM3(double sigma);

/**  
  // function to calculate the scaling factor of 
  // charge-distribution polynomial
  double calculateScalingFactor(double sigma); 

  // function to calculate the stretching factor of 
  // charge-distribution polynomial
  double calculateStretchingFactor(double sigma);

  // function to calculate the charge collected on strips
  // using the above calculated terms for one certain position in SD
  double calculateStripChargesPart(double scalF, double stretchF, double dist);

*/

  // function to calculate the charge collected on a strip using the erf function
  // It has been shown that the fitted erf function just depends on M3 
  // with a small error. 
  // n(x) = fabs( 50*(1 - erf( m3 * (fabs(DMin)) )) - 50*( 1 - erf( m3 * (fabs(DMax)) )) )
  // or
  // n(x) = fabs( 100 - 50*(1 - erf( m3 * (fabs(DMin)) )) - 50*( 1 - erf( m3 * (fabs(DMax)) )) )
  // depending on the position of the cloud respect with the readout electrode.
  
  double calculateStripChargesCollPart(double m3, double DMin , double DMax);

  // function to calculate stripcharges for all strips
  void calculateStripCharges(double* charges, double r, double z, 
			     double r_min, double r_max,T2Geometry* geom);

};

#endif
