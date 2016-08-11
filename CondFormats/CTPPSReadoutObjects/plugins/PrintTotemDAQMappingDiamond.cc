/****************************************************************************
 * Seyed Mohsen Etesami
 ****************************************************************************/

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CondFormats/DataRecord/interface/DiamondReadoutRcd.h"
#include "CondFormats/CTPPSReadoutObjects/interface/TotemDAQMappingDiamond.h"
#include "CondFormats/CTPPSReadoutObjects/interface/DiamondAnalysisMask.h"

//----------------------------------------------------------------------------------------------------

/**
 *\brief Prints the DAQ mapping loaded by DAQMappingSourceXML.
 **/
class PrintTotemDAQMappingDiamond : public edm::one::EDAnalyzer<>
{
public:
  PrintTotemDAQMappingDiamond(const edm::ParameterSet &ps) {}
  ~PrintTotemDAQMappingDiamond() {}

private:
  virtual void analyze(const edm::Event &e, const edm::EventSetup &es) override;
};

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

void PrintTotemDAQMappingDiamond::analyze(const edm::Event &e, const edm::EventSetup &es)
{
  // get mapping
  ESHandle<TotemDAQMappingDiamond> mapping;
  es.get<DiamondReadoutRcd>().get(mapping);

  // get analysis mask to mask channels
  ESHandle<DiamondAnalysisMask> analysisMask;
  es.get<DiamondReadoutRcd>().get(analysisMask);

  for (const auto &p : mapping->VFATMapping)
    {
      cout << p.first << " -> " << p.second << endl;
    }
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(PrintTotemDAQMappingDiamond);
