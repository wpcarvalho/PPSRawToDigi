/****************************************************************************
*
* CTPPS diamond timing detector reconstruction 
*
* Based on RecoCTPPS/TotemRPLocal/interface/TotemRPClusterProducerAlgorithm.h 
*
* Author:
*   Wagner Carvalho (wcarvalh@cern.ch)
*
****************************************************************************/

#ifndef RecoCTPPS_DiamondLocal_DiamondRecHitProducerAlgorithm
#define RecoCTPPS_DiamondLocal_DiamondRecHitProducerAlgorithm

#include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "DataFormats/CTPPSDigi/interface/DiamondDigi.h"
#include "DataFormats/CTPPSReco/interface/DiamondRecHit.h"

class DiamondRecHitProducerAlgorithm
{
  public:
    
    DiamondRecHitProducerAlgorithm(const edm::ParameterSet& param);

    ~DiamondRecHitProducerAlgorithm();
    
    void buildRecHit(unsigned int detId, const std::vector<DiamondDigi> &digi, std::vector<DiamondRecHit> &rechit);
    
  private:
    
    const edm::ParameterSet &param_;

};

#endif
