/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
* 	Hubert Niewiadomski
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#ifndef RecoTotemRP_RPRecoHitProducer_TotemRPRecHitProducerAlgorithm_h
#define RecoTotemRP_RPRecoHitProducer_TotemRPRecHitProducerAlgorithm_h

#include "Geometry/VeryForwardRPTopology/interface/RPTopology.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/CTPPSReco/interface/TotemRPCluster.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"
#include "RecoTotemRP/RPClusterSigmaService/interface/RPDetClusterSigmas.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TotemRPRecHitProducerAlgorithm
{
  public:
    TotemRPRecHitProducerAlgorithm(const edm::ParameterSet& conf) : cluster_sigmas_(0) {}
    void SetClusterSigmasService(const RPDetClusterSigmas *cluster_sigmas) {cluster_sigmas_=cluster_sigmas;}
    void BuildRecoHits(const edm::DetSet<TotemRPCluster>& input, 
        edm::DetSet<TotemRPRecHit>& output);
  private:
    const RPDetClusterSigmas *cluster_sigmas_;
    RPTopology rp_topology_;
};

#endif  //RecoTotemRP_RPRecoHitProducer_TotemRPRecHitProducerAlgorithm_h
