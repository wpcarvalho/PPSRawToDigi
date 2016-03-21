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

#include "EventFilter/TotemRawToDigi/interface/SimpleVFATFrameCollection.h"
#include "EventFilter/TotemRawToDigi/interface/RawDataUnpacker.h"
#include "EventFilter/TotemRawToDigi/interface/RawToDigiConverter.h"

#include <string>

//----------------------------------------------------------------------------------------------------

class TotemRawToDigi : public edm::EDProducer 
{
  public:
    explicit TotemRawToDigi(const edm::ParameterSet&);
    ~TotemRawToDigi();

    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob();

  private:
    edm::InputTag inputTag_;

    /// product labels
    std::string rpDataProductLabel;
    std::string rpCCProductLabel;
    std::string conversionStatusLabel;

    RawDataUnpacker rawDataUnpacker;
    RawToDigiConverter rawToDigiConverter;
};

//----------------------------------------------------------------------------------------------------

using namespace edm;
using namespace std;

//----------------------------------------------------------------------------------------------------

TotemRawToDigi::TotemRawToDigi(const edm::ParameterSet &conf):
  inputTag_((char const *)"rawDataCollector"),
  rawDataUnpacker(conf.getParameterSet("RawUnpacking")),
  rawToDigiConverter(conf.getParameterSet("RawToDigi"))
{
  produces<TotemRawEvent>();

  // RP data
  rpDataProductLabel = conf.getUntrackedParameter<std::string>("rpDataProductLabel", "");
  produces< edm::DetSetVector<RPStripDigi> > (rpDataProductLabel);

  // RP CC
  rpCCProductLabel = conf.getUntrackedParameter<std::string>("rpCCProductLabel", "");
  produces < std::vector <RPCCBits> > (rpCCProductLabel);

  // status
  conversionStatusLabel = conf.getUntrackedParameter<std::string>("conversionStatusLabel", "");
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

  // step 1: raw-data unpacking
  SimpleVFATFrameCollection vfatCollection;

  vector<int> fedIds = {}; // TODO

  for (const auto &fedId : fedIds)
  {
    rawDataUnpacker.Run(fedId, rawData->FEDData(fedId), vfatCollection, *totemRawEvent);
  }

  // step 2: raw to digi
  rawToDigiConverter.Run(vfatCollection, *mapping, *analysisMask,
    *rpDataOutput, *rpCCOutput, *conversionStatus);
  
  // commit products to event
  event.put(totemRawEvent);
  event.put(rpDataOutput, rpDataProductLabel);
  event.put(rpCCOutput, rpCCProductLabel);
  event.put(conversionStatus, conversionStatusLabel);
}

//----------------------------------------------------------------------------------------------------

void TotemRawToDigi::endJob()
{
  rawToDigiConverter.PrintSummaries();
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(TotemRawToDigi);
