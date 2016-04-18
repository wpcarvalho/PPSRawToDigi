/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
* 	Hubert Niewiadomski
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "RecoLocalCTPPS/TotemRP/interface/TotemRPRecHitProducerAlgorithm.h"

void TotemRPRecHitProducerAlgorithm::BuildRecoHits(const edm::DetSet<TotemRPCluster>& input, 
    edm::DetSet<TotemRPRecHit>& output)
{
  for (edm::DetSet<TotemRPCluster>::const_iterator it = input.begin(); it!=input.end(); ++it)
  {
    double sigma; // in mm
    if (cluster_sigmas_)
    {
      sigma = cluster_sigmas_->GetClusterSigma(it->DetId(), it->GetNumberOfStrips(), 
          it->CentreStripPos());
    } else {
      sigma = 0.0191;
    }
    output.push_back(TotemRPRecHit(it->DetId(), 
        rp_topology_.GetHitPositionInReadoutDirection(it->CentreStripPos()), sigma));
  }  
}
