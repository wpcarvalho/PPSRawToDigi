/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kaspar (jan.kaspar@gmail.com) 
*    
* $$RCSfile: RPTrackCounter.cc,v $: $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidateCollection.h"
//#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
//#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
//#include "DataFormats/Common/interface/DetSetVector.h"
//#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "DataFormats/TotemRPDataTypes/interface/RPRecoHit.h"
//#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecoElasticEvent.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
//#include "HepMC/GenEvent.h"

//#include "CLHEP/Vector/ThreeVector.h"
//#include "CLHEP/Vector/Rotation.h"

//#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"
//#include "Geometry/TotemRPGeometryBuilder/interface/DetGeomDesc.h"
//#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
//#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"

#include <map>
#include <set> 

/**
 *\brief  **/
class RPTrackCounter : public edm::EDAnalyzer
{
  public:
    RPTrackCounter(const edm::ParameterSet &ps);
    ~RPTrackCounter() {}

  private:
    edm::InputTag RPFittedTrackCollectionLabel;
    unsigned int verbosity;

    std::map<unsigned int, unsigned long> rpCount; 
    std::map<unsigned int, unsigned long> unitCount; 
    std::map<unsigned int, unsigned long> stationCount; 
    std::map< std::set<unsigned int>, unsigned long> confCount;

    virtual void analyze(const edm::Event &e, const edm::EventSetup &es);
    virtual void endJob();
};

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

RPTrackCounter::RPTrackCounter(const edm::ParameterSet &ps) : 
  verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0))
{
	RPFittedTrackCollectionLabel = ps.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");

}

//----------------------------------------------------------------------------------------------------

void RPTrackCounter::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  edm::Handle< RPFittedTrackCollection > fitTrCol;
  event.getByLabel(RPFittedTrackCollectionLabel, fitTrCol);

  set<unsigned int> rpSet, unitSet, stationSet;
  for (std::map<RPId, RPFittedTrack>::const_iterator it = fitTrCol->begin(); it != fitTrCol->end(); ++it) {
    if (!it->second.IsValid()) continue;

    unsigned int rpId = it->first;
    unsigned int stId = it->first / 10;
    unsigned int rpNum = it->first % 10;
    unsigned int unitId = stId * 10;
    if (rpNum > 2) unitId++;

    rpSet.insert(rpId);
    unitSet.insert(unitId);
    stationSet.insert(stId);
  }

  for (set<unsigned int>::iterator it = rpSet.begin(); it != rpSet.end(); ++it) rpCount[*it]++;
  for (set<unsigned int>::iterator it = unitSet.begin(); it != unitSet.end(); ++it) unitCount[*it]++;
  for (set<unsigned int>::iterator it = stationSet.begin(); it != stationSet.end(); ++it) stationCount[*it]++;
  confCount[rpSet]++;
}

//----------------------------------------------------------------------------------------------------

void RPTrackCounter::endJob()
{
  printf(">> RPTrackCounter::endJob\n");

  for (map< set<unsigned int>, unsigned long>::iterator it = confCount.begin(); it != confCount.end(); ++it) {
    printf("\t%6lu in RPs ", it->second);
    for (set<unsigned int>::iterator rit = it->first.begin(); rit != it->first.end(); ++rit) {
      if (rit != it->first.begin()) printf(", ");
      printf("%u", *rit);
    }
    printf("\n");
  }
  printf("\n");

  for (map<unsigned int, unsigned long>::iterator it = rpCount.begin(); it != rpCount.end(); ++it)
    printf("\t%6lu in RP %u\n", it->second, it->first);
  printf("\n");

  for (map<unsigned int, unsigned long>::iterator it = unitCount.begin(); it != unitCount.end(); ++it)
    printf("\t%6lu in unit %u\n", it->second, it->first);
  printf("\n");

  for (map<unsigned int, unsigned long>::iterator it = stationCount.begin(); it != stationCount.end(); ++it)
    printf("\t%6lu in station %u\n", it->second, it->first);
}



DEFINE_FWK_MODULE(RPTrackCounter);

