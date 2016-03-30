/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*   Hubert Niewiadomski
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "FWCore/Framework/interface/MakerMacros.h"
 
#include "RecoTotemRP/RPClusterizer/interface/RPClusterizer.h"

//----------------------------------------------------------------------------------------------------

RPClusterizer::RPClusterizer(edm::ParameterSet const& conf) :
  conf_(conf), RPClusterizerAlgorithm_(conf)
{
  verbosity_ = conf.getParameter<int>("Verbosity");

  digiInputTag_ = conf.getParameter<edm::InputTag>("DigiLabel");
  
  digiInputTagToken_ = consumes<edm::DetSetVector<TotemRPDigi> >(digiInputTag_);

  produces< edm::DetSetVector<TotemRPCluster> > ();
}

//----------------------------------------------------------------------------------------------------
 
RPClusterizer::~RPClusterizer()
{
}

//----------------------------------------------------------------------------------------------------
 
void RPClusterizer::beginJob()
{
}

//----------------------------------------------------------------------------------------------------
 
void RPClusterizer::produce(edm::Event& e, const edm::EventSetup& es)
{
  // get input
  edm::Handle< edm::DetSetVector<TotemRPDigi> > input;
 
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

void RPClusterizer::run(const edm::DetSetVector<TotemRPDigi>& input, edm::DetSetVector<TotemRPCluster> &output)
{
  for (const auto &ds_digi : input)
  {
    edm::DetSet<TotemRPCluster> &ds_cluster = output.find_or_insert(ds_digi.id);

    RPClusterizerAlgorithm_.BuildClusters(ds_digi.id, ds_digi.data, ds_cluster.data);
  }
}

DEFINE_FWK_MODULE(RPClusterizer);
