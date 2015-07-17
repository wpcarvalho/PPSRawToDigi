#ifndef VFatEfficiency_h
#define VFatEfficiency_h

/** 
 * Class VFatEfficiency
 * 
 * Author: Mirko Berretti / University of Siena
 * Email:  mirko.berretti@gmail.com
 * Date:   2008-08-10
 */


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"
#include <math.h>
#include "TMatrix.h"
#include "TVector.h"
#include "TH1F.h"
#include <cstdlib>
#include <boost/lexical_cast.hpp> 
#include "boost/shared_ptr.hpp"
#include <TH2.h>
#include <TFile.h>
#include <THnSparse.h>

class VFatEfficiency{

  
 public:
  
  //explicit PSimHitMisaligner(const edm::ParameterSet& paraSet);
  explicit VFatEfficiency(bool SetVfatEfficiency,std::string inputFileNameEffi,std::string inputFileNameCorrupt,std::string inputFileNameDeadChannels,std::string  inputFileNameNoisyChannels);
 
  virtual ~VFatEfficiency();


  // PSimHit DisplacePSimHit(const PSimHit &,uint32_t rawdetid);
  //Local3DPoint DisplacePoint(const Local3DPoint &,uint32_t rawid);
  //bool SimulationActivated(){return simulatemisalign;} 
  //typedef std::map<uint32_t,TMatrix> T1T2DetMatrixMap; //Associate Detector RawId to a Matrix Rotation/Translation  
  //void DisalignHit(T2Hit* t2h); 
  //T1T2DetMatrixMap T2DetMatrixRotation;
  //T1T2DetMatrixMap T2DetMatrixTranslation;
  // T1T2DetMatrixMap T1DetMatrixRotation;
  //T1T2DetMatrixMap T1DetMatrixTranslation;


  // TH1F TestDx;
  std::map<unsigned int,double> EffiMap;
  std::map<unsigned int,double> CorruptMap;
  std::map<unsigned int,std::vector<unsigned int> > VfatID_ToDeadChannelList;
  
  
  bool SetVfatEfficiency;

private:  
  
  std::string inputFileNameNoisyChannels;
  std::string inputFileNameEffi;
  std::string inputFileNameCorrupt;  
  std::string inputFileNameDeadChannels;
  
  bool verbosity;
  // Store T2  detectors  information
  std::vector<std::vector<double> > effvect;
  std::vector<std::vector<double> > corruptvect;
  //std::vector<std::vector<int> > deadChannelvect;
  //void ChangeDXDYSignAccordingPsimHitErikConv(TMatrix* matrixtrsl,uint32_t rawiddet);


 
};

#endif

