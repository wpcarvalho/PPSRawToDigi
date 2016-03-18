/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemL1Trigger/interface/RPCCBits.h"
#include "DataFormats/TotemRawData/interface/TotemRawEvent.h"
#include "DataFormats/TotemRawData/interface/TotemRawToDigiStatus.h"

#include "CondFormats/DataRecord/interface/TotemReadoutRcd.h"
#include "CondFormats/TotemReadoutObjects/interface/TotemDAQMapping.h"
#include "CondFormats/TotemReadoutObjects/interface/TotemAnalysisMask.h"

#include <string>

//----------------------------------------------------------------------------------------------------

class TotemRawToDigi : public edm::EDProducer 
{
  public:
    explicit TotemRawToDigi(const edm::ParameterSet&);
    ~TotemRawToDigi();

    virtual void produce(edm::Event&, const edm::EventSetup&) override;

  private:
    edm::InputTag inputTag_;

    /// product labels
    std::string rpDataProductLabel;
    std::string rpCCProductLabel;
    std::string conversionStatusLabel;
};

//----------------------------------------------------------------------------------------------------

using namespace edm;
using namespace std;

//----------------------------------------------------------------------------------------------------

TotemRawToDigi::TotemRawToDigi(const edm::ParameterSet &conf):
  inputTag_((char const *)"rawDataCollector")

/*
  verbosity(conf.getUntrackedParameter<unsigned int>("verbosity", 0)),
  printErrorSummary(conf.getUntrackedParameter<unsigned int>("printErrorSummary", 1)),
  printUnknownFrameSummary(conf.getUntrackedParameter<unsigned int>("printUnknownFrameSummary", 1)),

  testFootprint(conf.getParameter<unsigned int>("testFootprint")),
  testCRC(conf.getParameter<unsigned int>("testCRC")),
  testID(conf.getParameter<unsigned int>("testID")),
  testECRaw(conf.getParameter<unsigned int>("testECRaw")),
  testECDAQ(conf.getParameter<unsigned int>("testECDAQ")),
  testECMostFrequent(conf.getParameter<unsigned int>("testECMostFrequent")),
  testBCMostFrequent(conf.getParameter<unsigned int>("testBCMostFrequent")),
  
  EC_min(conf.getUntrackedParameter<unsigned int>("EC_min", 10)),
  BC_min(conf.getUntrackedParameter<unsigned int>("BC_min", 10)),
  
  EC_fraction(conf.getUntrackedParameter<double>("EC_fraction", 0.6)),
  BC_fraction(conf.getUntrackedParameter<double>("BC_fraction", 0.6))
*/

{
  produces<TotemRawEvent>();

  // RP data
  rpDataProductLabel = conf.getUntrackedParameter<std::string>("rpDataProductLabel", "rpDataOutput");
  produces< edm::DetSetVector<RPStripDigi> > (rpDataProductLabel);

  // RP CC
  rpCCProductLabel = conf.getUntrackedParameter<std::string>("rpCCProductLabel", "rpCCOutput");
  produces < std::vector <RPCCBits> > (rpCCProductLabel);

  // status
  conversionStatusLabel = conf.getUntrackedParameter<std::string>("conversionStatusLabel", "conversionStatus");
  produces <TotemRawToDigiStatus>(conversionStatusLabel);
}

//----------------------------------------------------------------------------------------------------

TotemRawToDigi::~TotemRawToDigi()
{
}

//----------------------------------------------------------------------------------------------------

void TotemRawToDigi::produce(edm::Event& event, const edm::EventSetup &es)
{
  // get DAQ mapping
  ESHandle<TotemDAQMapping> mapping;
  es.get<TotemReadoutRcd>().get(mapping);

  // get analysis mask to mask channels
  ESHandle<TotemAnalysisMask> analysisMask;
  es.get<TotemReadoutRcd>().get(analysisMask);

  // raw data handle
  edm::Handle<FEDRawDataCollection> rawData;
  event.getByLabel(inputTag_, rawData);

  // book output products
  auto_ptr<TotemRawEvent> totemRawEvent(new TotemRawEvent);

  auto_ptr< DetSetVector<RPStripDigi> > rpDataOutput(new edm::DetSetVector<RPStripDigi>);  
  auto_ptr< vector<RPCCBits> > rpCCOutput(new std::vector<RPCCBits>);
  auto_ptr< TotemRawToDigiStatus > conversionStatus(new TotemRawToDigiStatus());

  // uptodate example
  /*
  const FEDRawData &fedData = rawdata->FEDData(ScalersRaw::SCALERS_FED_ID);
  unsigned short int length = fedData.size();
  
  const ScalersEventRecordRaw_v6 *raw = (struct ScalersEventRecordRaw_v6 *)fedData.data();
  */

  // commit products to event
  event.put(totemRawEvent);
  event.put(rpDataOutput, rpDataProductLabel);
  event.put(rpCCOutput, rpCCProductLabel);
  event.put(conversionStatus, conversionStatusLabel);
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(TotemRawToDigi);
