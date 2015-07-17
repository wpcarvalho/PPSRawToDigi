/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#ifndef RPStationMultiTrackFinderFitter_IClusterBuilder
#define RPStationMultiTrackFinderFitter_IClusterBuilder

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatternsCollection.h"

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/Cluster.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/GeometryHelper.h"

#include <vector>

using namespace std;

namespace RPStationMultiTrackFinderFitter {

/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief Interface for \ref Cluster builders
 */
class IClusterBuilder
{
  public:
    explicit IClusterBuilder(const edm::ParameterSet &ps, GeometryHelper &geometry) :
        verbosity(ps.getUntrackedParameter<unsigned int>("verbosity")),
        geometry(geometry) {}

    virtual ~IClusterBuilder() {}

    /// should build the collection of \ref Cluster "Clusters" based on U/V patterns
    virtual void build(
            vector<Cluster> &clusters,
            const RPRecognizedPatternsCollection &rpPatColl) = 0;

  protected:
    const unsigned int verbosity;

    /// geometry based on which the clusters are found
    GeometryHelper &geometry;
};

} // namespace RPStationMultiTrackFinderFitter

#endif // #ifndef RPStationMultiTrackFinderFitter_IClusterBuilder
