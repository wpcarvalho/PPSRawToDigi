/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*   Hubert Niewiadomski
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
#include "DataFormats/TotemDigi/interface/TotemRPDigi.h"

#include "RecoLocalCTPPS/TotemRP/interface/TotemRPClusterProducerAlgorithm.h"
 
#include <iostream>

//----------------------------------------------------------------------------------------------------

/**
 * Merges neighbouring active TOTEM RP strips into clusters.
 **/
class TotemRPClusterProducer : public edm::one::EDProducer<>
{
  public:
  
    explicit TotemRPClusterProducer(const edm::ParameterSet& conf);
  
    virtual ~TotemRPClusterProducer();
  
    virtual void beginJob();
  
    virtual void produce(edm::Event& e, const edm::EventSetup& c);
  
  private:
    edm::ParameterSet conf_;
    int verbosity_;
    edm::InputTag digiInputTag_;
    edm::EDGetTokenT<edm::DetSetVector<TotemRPDigi>> digiInputTagToken_;

    TotemRPClusterProducerAlgorithm algorithm_;

    void run(const edm::DetSetVector<TotemRPDigi> &input, edm::DetSetVector<TotemRPCluster> &output);
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

TotemRPClusterProducer::TotemRPClusterProducer(edm::ParameterSet const& conf) :
  conf_(conf), algorithm_(conf)
{
  verbosity_ = conf.getParameter<int>("Verbosity");

  digiInputTag_ = conf.getParameter<edm::InputTag>("DigiLabel");
  
  digiInputTagToken_ = consumes<edm::DetSetVector<TotemRPDigi> >(digiInputTag_);

  produces< edm::DetSetVector<TotemRPCluster> > ();
}

//----------------------------------------------------------------------------------------------------
 
TotemRPClusterProducer::~TotemRPClusterProducer()
{
}

//----------------------------------------------------------------------------------------------------
 
void TotemRPClusterProducer::beginJob()
{
}

//----------------------------------------------------------------------------------------------------
 
void TotemRPClusterProducer::produce(edm::Event& e, const edm::EventSetup& es)
{
  // get input
  edm::Handle< edm::DetSetVector<TotemRPDigi> > input;
  e.getByToken(digiInputTagToken_, input);

  // prepare output
  std::auto_ptr< edm::DetSetVector<TotemRPCluster> > output(new edm::DetSetVector<TotemRPCluster>() );
  
  // run clusterisation
  if (verbosity_)
    std::cout << " Reading " << input->size() << " of TotemRPDigi" << std::endl;
 
  if (input->size())
    run(*input, *output);
   
  if (verbosity_)
    std::cout << " Saving " << output->size() << " of TotemRPCluster" << std::endl;

  // save output to event
  e.put(output);
}

//----------------------------------------------------------------------------------------------------

void TotemRPClusterProducer::run(const edm::DetSetVector<TotemRPDigi>& input, edm::DetSetVector<TotemRPCluster> &output)
{
  for (const auto &ds_digi : input)
  {
    edm::DetSet<TotemRPCluster> &ds_cluster = output.find_or_insert(ds_digi.id);

    algorithm_.BuildClusters(ds_digi.id, ds_digi.data, ds_cluster.data);
  }
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(TotemRPClusterProducer);
