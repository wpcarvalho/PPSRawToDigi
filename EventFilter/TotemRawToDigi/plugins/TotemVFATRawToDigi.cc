/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/TotemDigi/interface/TotemRPDigi.h"
#include "DataFormats/TotemDigi/interface/TotemVFATStatus.h"

#include "CondFormats/DataRecord/interface/TotemReadoutRcd.h"
#include "CondFormats/TotemReadoutObjects/interface/TotemDAQMapping.h"
#include "CondFormats/TotemReadoutObjects/interface/TotemAnalysisMask.h"

#include "EventFilter/TotemRawToDigi/interface/SimpleVFATFrameCollection.h"
#include "EventFilter/TotemRawToDigi/interface/RawDataUnpacker.h"
#include "EventFilter/TotemRawToDigi/interface/RawToDigiConverter.h"

#include <string>

//----------------------------------------------------------------------------------------------------

class TotemVFATRawToDigi : public edm::one::EDProducer<>
{
  public:
    explicit TotemVFATRawToDigi(const edm::ParameterSet&);
    ~TotemVFATRawToDigi();

    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob();

  private:
    std::string subSystem;

    // TODO: default values should be stored/read from
    //    DataFormats/FEDRawData/interface/FEDNumbering.h
  /* Hints from Michele
      DEVICE      ID     DETECTOR          OLD ID
      Trigger     577    LONEG             0x29c
      RX 1        578    5-6 210m FAR      0x1a1
      RX 2        579    5-6 210m NEAR     0x1a2
      RX 3        580    4-5 210m FAR      0x1a9
      RX 4        581    4-5 210m NEAR     0x1aa  
  */
    std::vector<unsigned int> fedIds;

    edm::EDGetTokenT<FEDRawDataCollection> fedDataToken;

    RawDataUnpacker rawDataUnpacker;
    RawToDigiConverter rawToDigiConverter;

    template <typename DigiType>
    void run(edm::Event&, const edm::EventSetup&);
};

//----------------------------------------------------------------------------------------------------

using namespace edm;
using namespace std;

//----------------------------------------------------------------------------------------------------

TotemVFATRawToDigi::TotemVFATRawToDigi(const edm::ParameterSet &conf):
  subSystem(conf.getParameter<string>("subSystem")),
  fedIds(conf.getParameter< vector<unsigned int> >("fedIds")),
  rawDataUnpacker(conf.getParameterSet("RawUnpacking")),
  rawToDigiConverter(conf.getParameterSet("RawToDigi"))
{
  fedDataToken = consumes<FEDRawDataCollection>(conf.getParameter<edm::InputTag>("rawDataTag"));

  // validate chosen subSystem
  if (subSystem != "RP")
    throw cms::Exception("TotemVFATRawToDigi::TotemVFATRawToDigi") << "Unknown sub-system string " << subSystem << "." << endl;

  // digi
  if (subSystem == "RP")
    produces< DetSetVector<TotemRPDigi> >(subSystem);

  // conversion status
  produces< DetSetVector<TotemVFATStatus> >(subSystem);
}

//----------------------------------------------------------------------------------------------------

TotemVFATRawToDigi::~TotemVFATRawToDigi()
{
}

//----------------------------------------------------------------------------------------------------

void TotemVFATRawToDigi::produce(edm::Event& event, const edm::EventSetup &es)
{
  if (subSystem == "RP")
    run< DetSetVector<TotemRPDigi> >(event, es);
}

//----------------------------------------------------------------------------------------------------

template <typename DigiType>
void TotemVFATRawToDigi::run(edm::Event& event, const edm::EventSetup &es)
{
  // get DAQ mapping
  ESHandle<TotemDAQMapping> mapping;
  es.get<TotemReadoutRcd>().get(mapping);

  // get analysis mask to mask channels
  ESHandle<TotemAnalysisMask> analysisMask;
  es.get<TotemReadoutRcd>().get(analysisMask);

  // raw data handle
  edm::Handle<FEDRawDataCollection> rawData;
  event.getByToken(fedDataToken, rawData);

  // book output products
  auto_ptr< DigiType > digi(new DigiType);  
  auto_ptr< DetSetVector<TotemVFATStatus> > conversionStatus(new DetSetVector<TotemVFATStatus>);

  // raw-data unpacking
  SimpleVFATFrameCollection vfatCollection;
  for (const auto &fedId : fedIds)
    rawDataUnpacker.Run(fedId, rawData->FEDData(fedId), vfatCollection);

  // raw-to-digi conversion
  rawToDigiConverter.Run(vfatCollection, *mapping, *analysisMask, *digi, *conversionStatus);

  // commit products to event
  event.put(digi, subSystem);
  event.put(conversionStatus, subSystem);
}

//----------------------------------------------------------------------------------------------------

void TotemVFATRawToDigi::endJob()
{
  rawToDigiConverter.PrintSummaries();
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(TotemVFATRawToDigi);
