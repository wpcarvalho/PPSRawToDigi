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
    void run(const edm::DetSetVector<TotemRPCluster>& input,std::vector<edm::DetSet<TotemRPRecHit> > & output);
    const edm::ParameterSet conf_;
    int verbosity_;
    TotemRPRecHitProducerAlgorithm algorithm_;

    // TODO: clean
    //std::string cluster_producer_;
    //std::string cluster_label_;
    //std::string rec_hit_label_;
    edm::InputTag cluster_label_;
    edm::EDGetTokenT<edm::DetSetVector<TotemRPCluster> >cluster_label_Token_;
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

TotemRPRecHitProducer::TotemRPRecHitProducer(const edm::ParameterSet& conf) :
  conf_(conf), algorithm_(conf)
{
  edm::LogInfo("TotemRPRecHitProducer") << "[TotemRPRecHitProducer::TotemRPRecHitProducer] Constructing object...";
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

//----------------------------------------------------------------------------------------------------
 
// Virtual destructor needed.
TotemRPRecHitProducer::~TotemRPRecHitProducer()
{
}
 
//----------------------------------------------------------------------------------------------------

void TotemRPRecHitProducer::beginJob()
{
}

//----------------------------------------------------------------------------------------------------
 
// Functions that gets called by framework every event
void TotemRPRecHitProducer::produce(edm::Event& e, const edm::EventSetup& es)
{
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

//----------------------------------------------------------------------------------------------------

void TotemRPRecHitProducer::run(const edm::DetSetVector<TotemRPCluster>& input,std::vector<edm::DetSet<TotemRPRecHit> > & output)
{
  int number_detunits = 0;
  int total_hits_no = 0;

  //loop on all detset inside the input collection
  edm::DetSetVector<TotemRPCluster>::const_iterator DSViter=input.begin();
  for (; DSViter!=input.end();DSViter++)
  {
    ++number_detunits;
    if(verbosity_)
      LogDebug("TotemRPRecHitProducer")  << "[TotemRPRecHitProducer::run] DetID " << DSViter->id;
 
      edm::DetSet<TotemRPRecHit> ssc(DSViter->id);
      algorithm_.BuildRecoHits(*DSViter, ssc);
      total_hits_no += ssc.data.size();
       
      if (ssc.data.size())
        output.push_back(ssc);  // insert the DetSet<TotemRPCluster> in the  DetSetVec<TotemRPCluster> only if there is at least a digi
  }
  if(verbosity_)
    LogDebug("TotemRPRecHitProducer") << "[TotemRPRecHitProducer] generating " << total_hits_no << " RPRecoHits in " << number_detunits << " DetUnits.";
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(TotemRPRecHitProducer);
