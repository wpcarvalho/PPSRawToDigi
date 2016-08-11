/****************************************************************************
 *  Seyed Mohsen Etesami
 ****************************************************************************/

#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "DataFormats/CTPPSDigi/interface/DiamondFEDInfo.h"

#include "EventFilter/CTPPSRawToDigi/interface/DiamondVFATFrameCollection.h"
#include "EventFilter/CTPPSRawToDigi/interface/DiamondRawDataUnpacker.h"

#include <string>

//----------------------------------------------------------------------------------------------------

class DiamondVFATFrameAnalyzer : public edm::global::EDAnalyzer<>
{ 
public:
  explicit DiamondVFATFrameAnalyzer(const edm::ParameterSet&);
  ~DiamondVFATFrameAnalyzer();

  virtual void analyze(edm::StreamID, const edm::Event &, const edm::EventSetup &) const override;

private:

  std::vector<unsigned int> fedIds;

  edm::EDGetTokenT<FEDRawDataCollection> fedDataToken;

  DiamondRawDataUnpacker rawDataUnpacker;

  template <typename DigiType>
  void run(edm::Event&, const edm::EventSetup&);
};

//----------------------------------------------------------------------------------------------------

using namespace edm;
using namespace std;

//----------------------------------------------------------------------------------------------------

DiamondVFATFrameAnalyzer::DiamondVFATFrameAnalyzer(const edm::ParameterSet &conf):
  fedIds(conf.getParameter< vector<unsigned int> >("fedIds")),
  rawDataUnpacker(conf.getParameterSet("RawUnpacking"))
{
  fedDataToken = consumes<FEDRawDataCollection>(conf.getParameter<edm::InputTag>("rawDataTag"));
}

//----------------------------------------------------------------------------------------------------

DiamondVFATFrameAnalyzer::~DiamondVFATFrameAnalyzer()
{
}

//----------------------------------------------------------------------------------------------------

void DiamondVFATFrameAnalyzer::analyze(edm::StreamID, const edm::Event& event, const edm::EventSetup &) const
{
  // raw data handle
  edm::Handle<FEDRawDataCollection> rawData;
  event.getByToken(fedDataToken, rawData);

  // raw-data unpacking
  vector<DiamondFEDInfo> fedInfo;
  DiamondVFATFrameCollection vfatCollection;
  for (const auto &fedId : fedIds)
    {
      const FEDRawData &data = rawData->FEDData(fedId);
      if (data.size() > 0)
	rawDataUnpacker.Run(fedId, data, fedInfo, vfatCollection);
    }

  // print VFAT frames
  cout << endl << "----------------------------------------------------------------------------------------------------" << endl;
  cout << event.id() << endl;

  for (DiamondVFATInterface::Iterator fr(&vfatCollection); !fr.IsEnd(); fr.Next())
    {
      cout << fr.Position() << " > ";
      fr.Data()->Print();
    }
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(DiamondVFATFrameAnalyzer);
