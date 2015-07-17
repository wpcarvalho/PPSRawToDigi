#ifndef T2RecHit_h
#define T2RecHit_h

/** 
 * Class T2RecHit
 * 
 * Author: Mirko Berretti / University of Siena
 * Email:  mirko.berretti@gmail.com
 * Date:   2008-08-10
 */

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include "T2DetHitReconst.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoTotemT1T2/T2RecHit/interface/T2DetHitReconst.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"
#include "Geometry/TotemT2AlignmentDataFormats/interface/T2AlignmentCorrections.h"
#include <cstdlib>
#include <boost/lexical_cast.hpp> 


class T2RecHit : public edm::EDProducer {

  
 public:
  
  explicit T2RecHit(const edm::ParameterSet& paraSet);

  virtual ~T2RecHit() {}

  virtual void produce(edm::Event& ev, const edm::EventSetup& evSet);

  //typedef std::map<uint32_t,TMatrix>
  //std::auto_ptr<T2AlignmentCorrections> T2CorrectionMap;
  T2AlignmentCorrections T2CorrectionMap;
  
  void CorrectHitPosition(T2Hit* ahit);

  virtual void beginJob() ;
  double myround(double d);


private:
  // boost::shared_ptr<T2DetHitReconst> theT2ClusterMatching; 
  //T2DetHitReconst* theT2ClusterMatching;
 bool includeClass0;
 bool InsertAlignmentbyCFG;
 bool CorrectWithResolution;

 bool checkdispsize;
 bool verbosity;
 unsigned int Cl1MaxPad;
 unsigned int Cl1MaxStrip;
 
 edm::ParameterSet RecHitProdParSet_;

 std::string ModuleLabelInput;
 std::string LabelproductInstanceName;

 DDTranslation S_m;        // zero translation by default  //ROOT::Math::DisplacementVector3D
 DDRotationMatrix R_m;     // identity rotation by default  //ROOT::Math::Rotation3D
 T2GeometryUtil conv;
 T2GeometryUtil::T2DetInfo planeinformation;
 std::string inputFileNameMisal;
 bool useTXTfile;
 
 std::vector<double>  DXdisp;
 std::vector<double>  DYdisp;
 std::vector<double > ParRminV;
 std::vector<double > ParRmaxV;


};

#endif
