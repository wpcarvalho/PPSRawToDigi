#ifndef T2_Pad_Electronics_Sim_h
#define T2_Pad_Electronics_Sim_h

/**
 * Class to simulate the electronics of the T2 pad
 * detector.
 *
 * Author: Erik Brücken / University of Helsinki
 * Email:  brucken@cc.helsinki.fi
 * Updateed:   2008-05-27
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"
#include <DataFormats/T2DigiVfat/interface/T2DigiVfat.h>
//#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/ESHandle.h"
#include "SimTotem/T2Digitizer/interface/VFatEfficiency.h"
#include <vector>
#include <map>
#include "CLHEP/Random/RandFlat.h"

#define NUMBER_OF_PCHANNELS 120

class T2DetectorHit;
class T2PadDigi;
class T2Geometry;

class T2PadElectronicsSim : public boost::noncopyable {

 public:

  T2PadElectronicsSim(const edm::ParameterSet & parameterSet, const edm::EventSetup& iSetup);

  void simulate(std::map<int, int*> & chargeMap);  

  void fillDigis(T2PadDigiCollection & padDigis,/* T2PadDigiCollection & deadPadActivated,*/
		 std::map<int, int*> & chargeMap);

  std::map<int, std::map<unsigned int, T2DigiVfat> > PadVFats;
  //boost::shared_ptr<std::map<int, std::map<unsigned int, T2DigiVfat> > > PadVFats;
  boost::shared_ptr<VFatEfficiency> theVFatEfficiency; 

  //Vector containing planes affected by shorts and correspinding sector;
  std::vector<unsigned int> VectDeadSect_Sector;
  std::vector<unsigned int> VectDeadSect_Plane;

 private:

  unsigned int RawtoSymb(uint32_t thedet);
  void readFile(std::string fileName, double* dataFile); 
  void LoadDeadSector(std::string inputFileNameDeadSect,std::vector<unsigned int> &VectDeadSect_Plane,std::vector<unsigned int> &VectDeadSect_Sector);
  bool IsPadInDeadSector(unsigned int RVal,unsigned int symbdetid);
  double EqThrFromEffi(double effi_measured);



  std::vector<int> bins_;
  std::vector<int> sigmaExtraNoise_;
  std::vector<int> simpleThreshold_;  
  std::vector<double> capaNoiseFactorPad_;

  bool SetVfatEfficiency;
  std::string inputFileNameEffi;
  std::string inputFileNameCorrupt;
  std::string inputFileNameDeadSect;
  std::string inputFileNameDeadChannels;
  std::string inputFileNameNoiseCorrelation;
  std::string inputFileNameNoisyChannels;

  double capaNoiseSigma[13][NUMBER_OF_PCHANNELS];
  double threshold[13][NUMBER_OF_PCHANNELS];
  bool UseCFGInfo;
  bool UseVFATs;
};

#endif
