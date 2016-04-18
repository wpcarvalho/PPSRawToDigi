/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
* 	Hubert Niewiadomski
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/CTPPSReco/interface/TotemRPCluster.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"

#include "RecoLocalCTPPS/TotemRP/interface/TotemRPRecHitProducerAlgorithm.h"
 
//----------------------------------------------------------------------------------------------------

class TotemRPRecHitProducer : public edm::one::EDProducer<>
{
  public:
  
    explicit TotemRPRecHitProducer(const edm::ParameterSet& conf);
  
    virtual ~TotemRPRecHitProducer();
  
    virtual void beginJob();
  
    virtual void produce(edm::Event& e, const edm::EventSetup& c);
  
  private:
    const edm::ParameterSet conf_;
    int verbosity_;

    TotemRPRecHitProducerAlgorithm algorithm_;

    edm::InputTag tagCluster_;
    edm::EDGetTokenT<edm::DetSetVector<TotemRPCluster>> tokenCluster_;
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

TotemRPRecHitProducer::TotemRPRecHitProducer(const edm::ParameterSet& conf) :
  conf_(conf), algorithm_(conf)
{
  verbosity_ = conf.getParameter<int>("Verbosity");

  tagCluster_ = conf.getParameter<edm::InputTag>("tagCluster");
  tokenCluster_ = consumes<edm::DetSetVector<TotemRPCluster> >(tagCluster_);

  produces<edm::DetSetVector<TotemRPRecHit>>();
}

//----------------------------------------------------------------------------------------------------
 
TotemRPRecHitProducer::~TotemRPRecHitProducer()
{
}
 
//----------------------------------------------------------------------------------------------------

void TotemRPRecHitProducer::beginJob()
{
}

//----------------------------------------------------------------------------------------------------
 
void TotemRPRecHitProducer::produce(edm::Event& e, const edm::EventSetup& es)
{
  // get input
  edm::Handle< edm::DetSetVector<TotemRPCluster> > input;
  e.getByToken(tokenCluster_, input);
 
  // prepare output
  DetSetVector<TotemRPRecHit> output;

  // build reco hits
  for (auto &ids : *input)
  {
    unsigned int rpId = ids.detId();

    DetSet<TotemRPRecHit> ods(rpId);
    algorithm_.BuildRecoHits(ids, ods);
       
    // insert the DetSet<TotemRPCluster> in the  DetSetVec<TotemRPCluster> only if there is at least one digi
    if (ods.data.size())
      output.insert(ods);
  }
   
  // save output
  e.put(make_unique<DetSetVector<TotemRPRecHit>>(output));
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(TotemRPRecHitProducer);
