
/****************************************************************************
* Seyed Mohsen Etesami
****************************************************************************/

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"

#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/CTPPSDigi/interface/DiamondDigi.h"
#include "DataFormats/CTPPSDigi/interface/DiamondVFATStatus.h"
#include "DataFormats/CTPPSDigi/interface/DiamondFEDInfo.h"

#include "CondFormats/DataRecord/interface/DiamondReadoutRcd.h"

#include "CondFormats/CTPPSReadoutObjects/interface/TotemDAQMappingDiamond.h"
#include "CondFormats/CTPPSReadoutObjects/interface/DiamondAnalysisMask.h"

#include "EventFilter/CTPPSRawToDigi/interface/DiamondVFATFrameCollection.h"
#include "EventFilter/CTPPSRawToDigi/interface/DiamondRawDataUnpacker.h"
#include "EventFilter/CTPPSRawToDigi/interface/DiamondRawToDigiConverter.h"

#include <string>

//----------------------------------------------------------------------------------------------------

class DiamondVFATRawToDigi : public edm::stream::EDProducer<>
{
  public:
    explicit DiamondVFATRawToDigi(const edm::ParameterSet&);
    ~DiamondVFATRawToDigi();

    virtual void produce(edm::Event&, const edm::EventSetup&) override;

  private:
    std::string subSystem;

    std::vector<unsigned int> fedIds;

    edm::EDGetTokenT<FEDRawDataCollection> fedDataToken;

    DiamondRawDataUnpacker rawDataUnpacker;
    DiamondRawToDigiConverter rawToDigiConverter;

    template <typename DigiType>
    void run(edm::Event&, const edm::EventSetup&);
};

//----------------------------------------------------------------------------------------------------

using namespace edm;
using namespace std;

//----------------------------------------------------------------------------------------------------

DiamondVFATRawToDigi::DiamondVFATRawToDigi(const edm::ParameterSet &conf):
  subSystem(conf.getParameter<string>("subSystem")),
  fedIds(conf.getParameter< vector<unsigned int> >("fedIds")),
  rawDataUnpacker(conf.getParameterSet("RawUnpacking")),
  rawToDigiConverter(conf.getParameterSet("RawToDigi"))
{
  fedDataToken = consumes<FEDRawDataCollection>(conf.getParameter<edm::InputTag>("rawDataTag"));

  // validate chosen subSystem
  if (subSystem != "RP")
    throw cms::Exception("DiamondVFATRawToDigi::DiaomndVFATRawToDigi") << "Unknown sub-system string " << subSystem << "." << endl;

  // FED (OptoRx) headers and footers
  produces< vector<DiamondFEDInfo> >(subSystem);

  // digi
    if (subSystem == "RP")
    produces< DetSetVector<DiamondDigi> >(subSystem);

  // set default IDs
  if (fedIds.empty())
  {
    if (subSystem == "RP")
    {
      for (int id = FEDNumbering::MINCTPPSDiamondsFEDID; id <= FEDNumbering::MAXCTPPSDiamondsFEDID; ++id)
        fedIds.push_back(id);
    }
  }

  // conversion status
  produces< DetSetVector<DiamondVFATStatus> >(subSystem);
}

//----------------------------------------------------------------------------------------------------

DiamondVFATRawToDigi::~DiamondVFATRawToDigi()
{
}

//----------------------------------------------------------------------------------------------------

void DiamondVFATRawToDigi::produce(edm::Event& event, const edm::EventSetup &es)
{
  if (subSystem == "RP")
    run< DetSetVector<DiamondDigi> >(event, es);
}

//----------------------------------------------------------------------------------------------------

template <typename DigiType>
void DiamondVFATRawToDigi::run(edm::Event& event, const edm::EventSetup &es)
{
  // get DAQ mapping
  ESHandle<TotemDAQMappingDiamond> mapping;
  es.get<DiamondReadoutRcd>().get(mapping);

  // get analysis mask to mask channels
  ESHandle<DiamondAnalysisMask> analysisMask;
 // es.get<DiamondReadoutRcd>().get(analysisMask);
  // raw data handle
  edm::Handle<FEDRawDataCollection> rawData;
  event.getByToken(fedDataToken, rawData);

  // book output products
  vector<DiamondFEDInfo> fedInfo;
  DigiType digi;
  DetSetVector<DiamondVFATStatus> conversionStatus;

  // raw-data unpacking
  DiamondVFATFrameCollection vfatCollection;
  for (const auto &fedId : fedIds)
  {
    const FEDRawData &data = rawData->FEDData(fedId);
    if (data.size() > 0)
      rawDataUnpacker.Run(fedId, data, fedInfo, vfatCollection);

  }

  // raw-to-digi conversion
  rawToDigiConverter.Run(vfatCollection, *mapping, *analysisMask, digi, conversionStatus);

  // commit products to event
  event.put(make_unique<vector<DiamondFEDInfo>>(fedInfo), subSystem);
  event.put(make_unique<DigiType>(digi), subSystem);
  event.put(make_unique<DetSetVector<DiamondVFATStatus>>(conversionStatus), subSystem);

}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(DiamondVFATRawToDigi);
