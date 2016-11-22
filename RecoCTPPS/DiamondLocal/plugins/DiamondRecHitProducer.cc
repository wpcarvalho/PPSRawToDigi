/****************************************************************************
*
* Plugin module to produce CTPPS diamond timing detector reconstructed hits 
* from digi hits.
*
* Based on RecoCTPPS/TotemRPLocal/plugins/TotemRPClusterProducer.cc 
*
* Author:
*   Wagner Carvalho (wcarvalh@cern.ch)
*
****************************************************************************/

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/CTPPSDigi/interface/DiamondDigi.h"

#include "RecoCTPPS/DiamondLocal/interface/DiamondRecHitProducerAlgorithm.h"
 
//----------------------------------------------------------------------------------------------------

/**
 * Produces RecHit from Digi.
 **/
class DiamondRecHitProducer : public edm::stream::EDProducer<>
{
  public:
  
    explicit DiamondRecHitProducer(const edm::ParameterSet& conf);
  
    virtual ~DiamondRecHitProducer() {}
  
    virtual void produce(edm::Event& e, const edm::EventSetup& c) override;
  
  private:
    edm::ParameterSet conf_;
    int verbosity_;
    std::string subSystem_;
    edm::InputTag digiInputTag_;
    edm::EDGetTokenT<edm::DetSetVector<DiamondDigi>> digiInputTagToken_;

    DiamondRecHitProducerAlgorithm algorithm_;

    void run(const edm::DetSetVector<DiamondDigi> &input, edm::DetSetVector<DiamondRecHit> &output);
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

DiamondRecHitProducer::DiamondRecHitProducer(edm::ParameterSet const& conf) :
  conf_(conf), algorithm_(conf)
{
  verbosity_ = conf.getParameter<int>("verbosity");
  
  subSystem_ = conf.getParameter<string>("subSystem");
  
  edm::LogInfo("CTPPS") << "DiamondRecHitProducer constructor with 'verbosity = " << verbosity_ 
                        << "' and 'subSystem = " << subSystem_ << "'" << std::endl;

  if(subSystem_ != "RP") throw cms::Exception("DiamondRecHitProducer::DiamondRecHitProducer") 
     << "Unknown subsystem string " << subSystem_ << "." << endl;
  
  digiInputTag_ = conf.getParameter<edm::InputTag>("tagDigi");
  digiInputTagToken_ = consumes<edm::DetSetVector<DiamondDigi> >(digiInputTag_);

  produces< edm::DetSetVector<DiamondRecHit> > (subSystem_);
}

//----------------------------------------------------------------------------------------------------
 
void DiamondRecHitProducer::produce(edm::Event& e, const edm::EventSetup& es)
{
  // get input
  edm::Handle< edm::DetSetVector<DiamondDigi> > input;
  e.getByToken(digiInputTagToken_, input);

  // prepare output
  DetSetVector<DiamondRecHit> output;
  
  // run reco
  if (input->size()) {
    run(*input, output);
  } else {
    edm::LogInfo("CTPPS") << "Input collection DetSetVector<DiamondDigi> is empty " << std::endl;
  }

  // save output to event
  e.put(make_unique<DetSetVector<DiamondRecHit>>(output),subSystem_);
}

//----------------------------------------------------------------------------------------------------

void DiamondRecHitProducer::run(const edm::DetSetVector<DiamondDigi>& input, edm::DetSetVector<DiamondRecHit> &output)
{
  for (const auto &ds_digi : input)
  {
    //  NEED TO IMPLEMENT THIS !!!
    edm::DetSet<DiamondRecHit> &ds_rechit = output.find_or_insert(ds_digi.id);

    //  NEED TO IMPLEMENT THIS !!!
    algorithm_.buildRecHit(ds_digi.id, ds_digi.data, ds_rechit.data);
  }
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(DiamondRecHitProducer);
