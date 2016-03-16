/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

// TODO: clean header files

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "DataFormats/Scalers/interface/L1AcceptBunchCrossing.h"
#include "DataFormats/Scalers/interface/L1TriggerScalers.h"
#include "DataFormats/Scalers/interface/Level1TriggerScalers.h"
#include "DataFormats/Scalers/interface/Level1TriggerRates.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/Scalers/interface/BeamSpotOnline.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/Scalers/interface/ScalersRaw.h"

//----------------------------------------------------------------------------------------------------

class TotemRawToDigi : public edm::EDProducer 
{
  public:
    explicit TotemRawToDigi(const edm::ParameterSet&);
    ~TotemRawToDigi();

    virtual void produce(edm::Event&, const edm::EventSetup&) override;

  private:
    edm::InputTag inputTag_;
};

//----------------------------------------------------------------------------------------------------

TotemRawToDigi::TotemRawToDigi(const edm::ParameterSet& iConfig):
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

  produces<L1AcceptBunchCrossingCollection>();
}

//----------------------------------------------------------------------------------------------------

TotemRawToDigi::~TotemRawToDigi()
{
}

//----------------------------------------------------------------------------------------------------

// Method called to produce the data 
void TotemRawToDigi::produce(edm::Event& event, const edm::EventSetup& /*eventSetup*/)
{
  using namespace edm;

  // Get a handle to the FED data collection
  edm::Handle<FEDRawDataCollection> rawdata;
  event.getByLabel(inputTag_, rawdata);

  //std::auto_ptr<DcsStatusCollection> pDcsStatus(new DcsStatusCollection());

  //const FEDRawData &fedData = rawdata->FEDData(ScalersRaw::SCALERS_FED_ID);
  //unsigned short int length = fedData.size();
  
  //const ScalersEventRecordRaw_v6 *raw = (struct ScalersEventRecordRaw_v6 *)fedData.data();


  //iEvent.put(pDcsStatus);
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(TotemRawToDigi);
