// File: RPRecoHitProducer.cc
// Description:  see RPRecoHitProducer.h
// Author:  Hubert Niewiadomski, CERN
// Creation Date:  OGU Aug. 1 2005 Initial version.
//
//--------------------------------------------


#include "RecoTotemRP/RPRecoHitProducer/interface/RPRecoHitProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

RPRecoHitProducer::RPRecoHitProducer(const edm::ParameterSet& conf) :
  conf_(conf), RPRecoHitProducerAlgorithm_(conf)
{
  edm::LogInfo("RPRecoHitProducer") << "[RPRecoHitProducer::RPRecoHitProducer] Constructing object...";
  verbosity_ = conf.getParameter<int>("Verbosity");
//  cluster_producer_ = conf.getParameter<std::string>("ClusterProducer");
  if(!conf.exists("ClusterLabel")){
	  std::cout<<"expecting ClusterLabel parameter for type: TotemRPCluster"<<endl;
  }
  	  cluster_label_ = conf.getParameter<edm::InputTag>("ClusterLabel");
//  rec_hit_label_ = conf.getParameter<std::string>("RecHitLabel");
//  produces< edm::DetSetVector<TotemRPRecHit> > (rec_hit_label_);
  produces< edm::DetSetVector<TotemRPRecHit> > ();
  cluster_label_Token_ = consumes<edm::DetSetVector<TotemRPCluster> >(cluster_label_);


}
 
// Virtual destructor needed.
RPRecoHitProducer::~RPRecoHitProducer()
{
  edm::LogInfo("RPRecoHitProducer") << "[RPRecoHitProducer::~RPRecoHitProducer] Destructing object...";
}
 
//Get at the beginning
void RPRecoHitProducer::beginJob()
{
  if(verbosity_)
  {
    edm::LogInfo("RPRecoHitProducer") << "[RPRecoHitProducer::beginJob]";
  }
  RPRecoHitProducerAlgorithm_.SetClusterSigmasService( &the_cluster_sigma_ );
}
 
// Functions that gets called by framework every event
void RPRecoHitProducer::produce(edm::Event& e, const edm::EventSetup& es)
{
//  es.get<RPClusterSigmaServiceRecord>().get(cluster_sigmas_);
//  RPRecoHitProducerAlgorithm_.SetClusterSigmasService( &(*cluster_sigmas_) );
  
  // Step B: Get Inputs
  edm::Handle< edm::DetSetVector<TotemRPCluster> > input;
 
  // Step C: produce output product
  std::vector< edm::DetSet<TotemRPRecHit> > vRPRecoHits;
  vRPRecoHits.reserve(240);
  
 // e.getByLabel(cluster_label_, input);  //FIXME: fix this label
//  e.getByType(input);
    e.getByToken(cluster_label_Token_, input);
 

 
  if(input->size())
    run(*input,vRPRecoHits);
   
  // Step D: create and fill output collection
  std::auto_ptr< edm::DetSetVector<TotemRPRecHit> > output(new edm::DetSetVector<TotemRPRecHit>(vRPRecoHits) );
 
  // Step D: write output to file
  e.put(output);
}


void RPRecoHitProducer::run(const edm::DetSetVector<TotemRPCluster>& input,std::vector<edm::DetSet<TotemRPRecHit> > & output)
{
  int number_detunits = 0;
  int total_hits_no = 0;

  //loop on all detset inside the input collection
  edm::DetSetVector<TotemRPCluster>::const_iterator DSViter=input.begin();
  for (; DSViter!=input.end();DSViter++)
  {
    ++number_detunits;
    if(verbosity_)
      LogDebug("RPRecoHitProducer")  << "[RPRecoHitProducer::run] DetID " << DSViter->id;
 
      edm::DetSet<TotemRPRecHit> ssc(DSViter->id);
      RPRecoHitProducerAlgorithm_.BuildRecoHits(*DSViter, ssc);
      total_hits_no += ssc.data.size();
       
      if (ssc.data.size())
        output.push_back(ssc);  // insert the DetSet<TotemRPCluster> in the  DetSetVec<TotemRPCluster> only if there is at least a digi
  }
  if(verbosity_)
    LogDebug("RPRecoHitProducer") << "[RPRecoHitProducer] generating " << total_hits_no << " RPRecoHits in " << number_detunits << " DetUnits.";
}

DEFINE_FWK_MODULE(RPRecoHitProducer);
