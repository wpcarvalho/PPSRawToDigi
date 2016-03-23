#ifndef RecoTotemRP_RPRecoHitProducer_RPRecoHitProducerAlgorithm_h
#define RecoTotemRP_RPRecoHitProducer_RPRecoHitProducerAlgorithm_h

#include "Geometry/TotemRPDetTopology/interface/RPTopology.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/TotemRPReco/interface/TotemRPCluster.h"
#include "DataFormats/TotemRPReco/interface/TotemRPRecHit.h"
#include "RecoTotemRP/RPClusterSigmaService/interface/RPDetClusterSigmas.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class RPRecoHitProducerAlgorithm
{
  public:
    RPRecoHitProducerAlgorithm(const edm::ParameterSet& conf) : cluster_sigmas_(0) {}
    void SetClusterSigmasService(const RPDetClusterSigmas *cluster_sigmas) {cluster_sigmas_=cluster_sigmas;}
    void BuildRecoHits(const edm::DetSet<TotemRPCluster>& input, 
        edm::DetSet<TotemRPRecHit>& output);
  private:
    const RPDetClusterSigmas *cluster_sigmas_;
    RPTopology rp_topology_;
};

#endif  //RecoTotemRP_RPRecoHitProducer_RPRecoHitProducerAlgorithm_h
