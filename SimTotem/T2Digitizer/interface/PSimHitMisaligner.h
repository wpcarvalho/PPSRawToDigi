#ifndef PSimHitMisaligner_h
#define PSimHitMisaligner_h

/** 
 * Class PSimHitMisaligner
 * 
 * Author: Mirko Berretti / University of Siena
 * Email:  mirko.berretti@gmail.com
 * Date:   2008-08-10
 */


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include "DataFormats/T2Hit/interface/T2HitCollection.h"
#include "Geometry/TotemT2AlignmentDataFormats/interface/T2AlignmentCorrections.h"
//#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "DetectorDescription/Base/interface/DDRotationMatrix.h"
#include "DetectorDescription/Base/interface/DDTranslation.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"

#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
#include <math.h>
#include "TMatrix.h"
#include "TVector.h"
#include "TH1F.h"
#include "CLHEP/Random/RandGaussQ.h"
#include <cstdlib>
#include <boost/lexical_cast.hpp> 
class PSimHit;

class PSimHitMisaligner{

  
 public:
  
  //explicit PSimHitMisaligner(const edm::ParameterSet& paraSet);
  explicit PSimHitMisaligner(const edm::ParameterSet& paraSet, const edm::EventSetup &iSetup);
 
  virtual ~PSimHitMisaligner();


  PSimHit DisplacePSimHit(const PSimHit &,uint32_t rawdetid);


  Local3DPoint DisplacePoint(const Local3DPoint &,uint32_t rawid);


  bool SimulationActivated(){return simulatemisalign;} 

  typedef std::map<uint32_t,TMatrix> T1T2DetMatrixMap; //Associate Detector RawId to a Matrix Rotation/Translation
  
  //void DisalignHit(T2Hit* t2h);
  

 
  T1T2DetMatrixMap T2DetMatrixRotation;
  T1T2DetMatrixMap T2DetMatrixTranslation;
  T1T2DetMatrixMap T1DetMatrixRotation;
  T1T2DetMatrixMap T1DetMatrixTranslation;


  TH1F TestDx;

private:

  
  std::string inputFileNameMisal;
  bool generaterandom,simulatemisalign,verbosity;
  double sigmaDispl;
  double sigmaPhiDispl;
  double sigmaGlobalThetaXY;
  double sigmaGlobalShiftXY;
  
  // Store T2  detectors disalignment information
  std::vector<std::vector<double> >dispvect;
 
  void ChangeDXDYSignAccordingPsimHitErikConv(TMatrix* matrixtrsl,uint32_t rawiddet);
  void ChangePhiSignAccordingPsimHitErikConv(TMatrix* matrixrot,uint32_t rawiddet);
 
};

#endif

