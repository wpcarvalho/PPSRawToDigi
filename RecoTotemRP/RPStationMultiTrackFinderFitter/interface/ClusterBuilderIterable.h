/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#ifndef RPStationMultiTrackFinderFitter_ClusterBuilderIterable
#define RPStationMultiTrackFinderFitter_ClusterBuilderIterable

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatternsCollection.h"

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/GeometryHelper.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/IClusterBuilder.h"


namespace RPStationMultiTrackFinderFitter {

/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief Creates and merges \ref Cluster "Clusters" in sequential manner
 */
class ClusterBuilderIterable : public IClusterBuilder
{
  public:
    explicit ClusterBuilderIterable(const edm::ParameterSet &ps, GeometryHelper &geometry);

    /**
     * Iterating over all possible track candidates (combinations of U and V patterns
     * from two RPs) it creates and merges them with \ref Cluster "Clusters".
     *
     * \deprecated
     * The process by itself is dependent on the order of analysis. It is possible to
     * produce multiple similar \ref Cluster "Clusters".
     */
    virtual void build(
            vector<Cluster> &clusters,
            const RPRecognizedPatternsCollection &rpPatColl);
  private:
    /// cut offs to reject apparently wrong track candidates
    const double ax_cut, ay_cut;

    /// the projection plane position
    const double z_h;

    /// threshold of distance between track candidate and \ref Cluster
    /// for them to be merged together
    const double cluster_merge_si_cut;
};

} // namespace RPStationMultiTrackFinderFitter

#endif // #ifndef RPStationMultiTrackFinderFitter_ClusterBuilderIterable
