/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*   Hubert Niewiadomski
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#ifndef RecoTotemRP_RPClusterizer_RPClusterizerAlgorithm
#define RecoTotemRP_RPClusterizer_RPClusterizerAlgorithm

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "DataFormats/TotemDigi/interface/TotemRPDigi.h"
#include "DataFormats/TotemRPReco/interface/TotemRPCluster.h"

#include "TVector3.h"

#include <vector>
#include <set>

class RPClusterizerAlgorithm
{
  public:
    RPClusterizerAlgorithm(const edm::ParameterSet& param);

    ~RPClusterizerAlgorithm();
    
    int BuildClusters(unsigned int detId, const std::vector<TotemRPDigi> &digi, std::vector<TotemRPCluster> &clusters);
    
  private:
    typedef std::set<TotemRPDigi> TotemRPDigiSet;

    void Reset();

    TotemRPDigiSet strip_digi_set_;  ///< input digi set, strip by strip

    const edm::ParameterSet &param_;

    int verbosity_;
};

#endif
