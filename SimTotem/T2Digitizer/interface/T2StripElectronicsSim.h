#ifndef T2_Strip_Electronics_Sim_h
#define T2_Strip_Electronics_Sim_h

/**
 * Class to simulate the electronics of the T2 strip
 * detector.
 *
 * Author:  Erik Br??cken / University of Helsinki
 * Email:   brucken@cc.helsinki.fi
 * Updated: 2008-05-27
 */



#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/T2Digi/interface/T2StripDigiCollection.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"
#include <DataFormats/T2DigiVfat/interface/T2DigiVfat.h>
//#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "SimTotem/T2Digitizer/interface/VFatEfficiency.h"

#include <vector>
#include <map>
#include "CLHEP/Random/RandFlat.h"

//Random Number
#include <cstdlib> // I need it for random numbers
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CLHEP/Random/RandomEngine.h"
#include <TRandom3.h>








#define NUMBER_OF_CHANNELS 256

class T2DetectorHit;
class T2StripDigi;
class T2Geometry;

class T2StripElectronicsSim {

 public:

  T2StripElectronicsSim(const edm::ParameterSet & parameterSet,const edm::EventSetup& iSet,CLHEP::HepRandomEngine* rndEngine);

  void simulate(std::map<int, int*> & chargeMap);  

  void fillDigis(T2StripDigiCollection & stripDigis, 
		 std::map<int, int*> & chargeMap);
 
  std::map<int, std::map<unsigned int, T2DigiVfat> > StripVFats;
  boost::shared_ptr<VFatEfficiency> theVFatEfficiency; 
  
  std::vector<unsigned int> VectDeadSect_Sector;
  std::vector<unsigned int> VectDeadSect_Plane;
  CLHEP::HepRandomEngine* rndEngineStr;
 private:

  unsigned int RawtoSymb(uint32_t thedet);
  void readFile(std::string fileName, double* dataFile); 
  void LoadDeadSector(std::string inputFileNameDeadSect,std::vector<unsigned int> &VectDeadSect_Plane,std::vector<unsigned int> &VectDeadSect_Sector);
  bool IsStripInDeadSector(unsigned int RVal,unsigned int symbdetid);
  double EqThrFromEffi(double effi_measured);


  std::vector<int> bins_;
  std::vector<int> sigmaExtraNoise_;
  std::vector<int> simpleThreshold_;
  std::vector<double> capaNoiseFactorStrip_;
  double capaNoiseSigma[2][NUMBER_OF_CHANNELS];
  double threshold[2][NUMBER_OF_CHANNELS];
  bool UseCFGInfo;
  bool UseVFATs;

  bool SetVfatEfficiency;
  std::string inputFileNameEffi;
  std::string inputFileNameCorrupt;
  std::string inputFileNameDeadSect;
  std::string inputFileNameDeadChannels;
  std::string inputFileNameNoiseCorrelation;
  std::string inputFileNameNoisyChannels;

  
  boost::shared_ptr<TFile> MCfileNoise;  
  THnSparseD* EvtVfat_Strip_WhenCompletelyOn;
  
};

#endif
