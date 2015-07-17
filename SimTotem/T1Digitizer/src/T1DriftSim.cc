#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimTotem/T1Digitizer/interface/T1DriftSim.h"
#include "SimTotem/T1Digitizer/interface/T1DetectorHit.h"

#include "Geometry/TotemGeometry/interface/T1ChamberSpecs.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include <cmath>
#include <iostream>

using CLHEP::RandGaussQ;
//#define _LOGINFO_

T1DriftSim::T1DriftSim() 
  : STEP_SIZE(0.01), 
    ELECTRON_DIFFUSION_COEFF(0.0161), ///??????
    theMagneticField(0)
{
  // just initialize avalanche sim.  There has to be a better
  // way to take the integral of a function!
  static const int   N_INTEGRAL_STEPS = 700;
  double sum = 0.;
  // tmpMap is the unnormalized dNdEIntegral
  dNdEIntegral.resize(N_INTEGRAL_STEPS);
  int i;
  for(i = 0; i < N_INTEGRAL_STEPS; ++i) {
    double xx = STEP_SIZE * (double(i) - 0.5 );
    double dNdE = pow( xx, 0.38) * exp(-1.38*xx);
    if(i > 1) {
      sum += dNdE;
    }
    // store this value in the map
    dNdEIntegral[i] = sum;
  }
  // now normalize the whole map
  for(i =  0; i < N_INTEGRAL_STEPS; ++i) {
    dNdEIntegral[i] /= sum;
  }
}


T1DetectorHit 
T1DriftSim::getWireHit(const Local3DPoint & pos, T1Geometry * layer,
		       int nearestWire, const PSimHit & simHit) {
  //carica le specifiche della camera 

  //  std::cout << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ " << simHit.detUnitId() << " " << simHit.tof() << std::endl;

  T1ChamberSpecs specs;  

  // pos.x e' come il simHit

  HepGeom::Point3D<double> clusterPos( pos.x()-layer->xOfWire(simHit.detUnitId()), pos.y()-layer->yOfWire(simHit.detUnitId(),nearestWire), pos.z());
#ifdef _LOGINFO_
  edm::LogInfo("T1DriftSim") << "Ionization cluster at: " <<  pos.x()<<" "<< pos.y()<<" "<<pos.z();
  edm::LogInfo("T1DriftSim") << "ClusterPos: " <<pos.x()<<"-"<<(layer->xOfWire(simHit.detUnitId()) ) <<" "<<pos.y()<<"-"<<(layer->yOfWire(simHit.detUnitId(),nearestWire)) <<" "<<pos.z();
#endif


  // set the coordinate system with the x-axis along the nearest wire,
  // with the origin in the center of the chamber, on that wire.

  //HepTranslateY3D yShift(-1.*geom->yOfWire(nearestWire));
  //HepRotateZ3D rotation(-1.*geom->wireAngle());
  //clusterPos = yShift * clusterPos;
  // clusterPos = rotation * clusterPos;
  //usato per il campo magnetico....
  // GlobalPoint globalPosition = layer->surface().toGlobal(pos);
  //  assert(theMagneticField != 0);
 
  //  bz = theMagneticField->inTesla(globalPosition).z() * 10.;

  // We need magnetic field in _local_ coordinates
  // Interface now allows access in kGauss directly.
  bz=0.0000001;
  //  bz = layer->toLocal( theMagneticField->inKGauss( globalPosition ) ).z();

  // these subroutines label the coordinates in GEANT coords...
  //parametri che sono usati dalla funzioni avgPathLengthLowB ecc...
  //  ycell = clusterPos.z() / specs->anodeCathodeSpacing(); ///
  // zcell = 2.*clusterPos.y() / specs->wireSpacing();  ///

  ycell = clusterPos.z() / specs.anodeCathodeSpacing(); ///
  zcell = 2.*clusterPos.y() / specs.wireSpacing();  ///

#ifdef _LOGINFO_
  edm::LogInfo("T1DriftSim") << "bz " << bz <<" avgDrift " << avgDrift()
			     << " ycell " << ycell << " zcell " << zcell;
#endif

  double avgPathLength, pathSigma, avgDriftTime, driftTimeSigma;
  static const float B_FIELD_CUT = 15.f;
  if(fabs(bz) < B_FIELD_CUT) {
    //modificare i parametri in prima apprx senza il campo magnetico
    avgPathLength  = avgPathLengthLowB();
    pathSigma      = pathSigmaLowB();
    avgDriftTime   = avgDriftTimeLowB();
    driftTimeSigma = driftTimeSigmaLowB();
  }
  else {
    avgPathLength  = avgPathLengthHighB();
    pathSigma      = pathSigmaHighB();
    avgDriftTime   = avgDriftTimeHighB();
    driftTimeSigma = driftTimeSigmaHighB();
  }

  // electron drift path length 
  double pathLength = RandGaussQ::shoot(avgPathLength, pathSigma);

  // electron drift distance along the anode wire, including diffusion
  double diffusionSigma = ELECTRON_DIFFUSION_COEFF * sqrt(pathLength);
  double x = clusterPos.x() + RandGaussQ::shoot(avgDrift(), driftSigma())
    + RandGaussQ::shoot(0., diffusionSigma);
#ifdef _LOGINFO_
  edm::LogInfo("T1DriftSim")<<" clusterPos.x() + RandGaussQ::shoot(avgDrift(), driftSigma()) + RandGaussQ::shoot(0., diffusionSigma) "
		  << clusterPos.x()<< " " <<  RandGaussQ::shoot(avgDrift(), driftSigma()) << " " << RandGaussQ::shoot(0., diffusionSigma)<<" "<< avgDrift()<<" "<< driftSigma()<<" "<< diffusionSigma<< " " <<ELECTRON_DIFFUSION_COEFF << " "<<pathLength <<" "<<avgPathLength<<" "<< pathSigma;
#endif

  // electron drift time
  double driftTime  = RandGaussQ::shoot(avgDriftTime, driftTimeSigma);

  //@@ Parameters which should be defined outside the code
  // f_att is the fraction of drift electrons lost due to attachment

  //noti...
  //  static const double f_att = 0.5;
  // static const double f_collected = 0.82;

  static const double f_att = 0.9;

  static const double f_sha = 1.;
  // Avalanche charge, with fluctuation ('avalancheCharge()' is the fluctuation generator!)

  //double charge = avalancheCharge() * f_att * f_collected * specs->gasGain() * e_SI * 1.e15;////
  //        f_ind ??????????????????????????????????
  double charge = avalancheCharge() * f_att * f_sha * specs.gasGain() * CLHEP::e_SI * 1.e15;////
  ///1.e15 e' probabilmente dovuto al fatto che la carica viene espressa, per comodita', in fC e non in C
 
  //%// std::cout << " CHARGE " << charge << std::endl;

  float t = simHit.tof() + driftTime;
  


#ifdef _LOGINFO_
  T1DetectorHit DH(nearestWire, charge, x, t, &simHit, simHit.detUnitId());
  edm::LogInfo("T1DriftSim") << "--------------------------------------------------------------------- T1DetectorHit("<<nearestWire<<" " << charge<<" " << x<<" " << t<<" " << simHit.detUnitId();
  edm::LogInfo("T1DriftSim") << "T1DriftSim: tof = " << simHit.tof() << 
    " driftTime = " << driftTime <<
    " MEDH = "<< "element: " << DH.getElement()
			     << "  charge: " << DH.getCharge()
			     << "   pos:  " << DH.getPosition()
			     << "   time: " << DH.getTime()   
			     <<   "  DetId "  << DH.getDetId();

#endif
  return T1DetectorHit(nearestWire, charge, x, t, &simHit,simHit.detUnitId() );
}



// Generate avalanche fluctuation
#include <algorithm>
double T1DriftSim::avalancheCharge() {
  double returnVal = 0.;
  // pick a random value along the dNdE integral
  double x = CLHEP::RandFlat::shoot();
  size_t i;
  size_t isiz = dNdEIntegral.size();
  /*
    for(i = 0; i < isiz-1; ++i) {
    if(dNdEIntegral[i] > x) break;
    }
  */
  // return position of first element with a value >= x
  std::vector<double>::const_iterator p=lower_bound(dNdEIntegral.begin(),dNdEIntegral.end(),x);
  if (p==dNdEIntegral.end()) i=isiz-1;
  else i = p-dNdEIntegral.begin();
			      

  // now extrapolate between values
  if( i ==  isiz-1 ) {
    edm::LogInfo("T1DriftSim") << "Funky integral in T1DriftSim ";
    returnVal = STEP_SIZE * double(i) * dNdEIntegral[i];
  }
  else {
    double x1 = dNdEIntegral[i];
    double x2 = dNdEIntegral[i+1];
    returnVal = STEP_SIZE * (double(i) + (x-x1)/(x2-x1)); 
  }
#ifdef _LOGINFO_
  edm::LogInfo("T1DriftSim") << "avalanche fluc " << returnVal << "  " << x ;
#endif

  return returnVal;  
}

