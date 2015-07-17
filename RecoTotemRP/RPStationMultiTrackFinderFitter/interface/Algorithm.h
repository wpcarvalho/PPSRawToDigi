/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#ifndef RecoTotemRP_RPStationMultiTrackFinderFitter_Algorithm_h
#define RecoTotemRP_RPStationMultiTrackFinderFitter_Algorithm_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPStationTrackFitCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatternsCollection.h"

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/OverlapsRemoval.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/Cluster.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/GeometryHelper.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/IClusterBuilder.h"

#include "TMatrixD.h"
#include "TVectorD.h"

#include <vector>
#include <map>
#include <set>
#include <string>
#include <cstdio>
#include <cmath>
#include <utility>

using namespace std;

namespace RPStationMultiTrackFinderFitter {

/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief It is responsible for running of the stages of reconstruction
 *
 * The algorithm is expecting an input in form of U/V patterns being
 * reconstructed by 1-RP reconstruction algorithm of type
 * \ref RPRecognizedPatternsCollection and returns the reconstructed tracks
 * in object of type \ref RPStationTrackFitCollection.
 *
 * \author Jan Kaspar <jan.kaspar@gmail.com>
 * \author Jakub Sawicki <jakub.kuba.sawicki@gmail.com>
 */
class Algorithm
{
  public:
    explicit Algorithm(
        const edm::ParameterSet& ps,
        GeometryHelper &geometry,
        IClusterBuilder &clusterBuilder);

    virtual ~Algorithm();

    /**
     * \brief Performs the tracks reconstruction
     *
     * U/V patterns are fed into the object implementing \ref IClusterBuilder
     * in order to obtain the set of track candidate \ref Cluster "Clusters".
     * These are filtered based on how many RPs are included per Cluster and
     * distances from the U/V patterns. Due to the fact that resulting Clusters
     * might not be compatible with each other their population is later subjected
     * to the \ref OverlapsRemoval algorithm. What remains is being returned.
     */
    virtual void reconstruct(
            /// [out] collection of Clusters returned from the algorithm
            vector<Cluster> *output,
            /// [in] U/V patterns collection
            const RPRecognizedPatternsCollection *input);

protected:
    const unsigned int verbosity;
    GeometryHelper &geometry;
    IClusterBuilder &clusterBuilder;
    OverlapsRemoval overlapsRemoval;

    const double filter_dist_cut; ///< minimal squared sigma distance from U/V pattern
                                  ///< for a Cluster to be rejected

    const double track_pattern_assoc_cut; ///< maximal squared sigma distance from U/V
                                          ///< pattern for a Cluster to be considered
                                          ///< for association

    const double z_h; ///< the plane at which a and b track parameters are projected
    const unsigned int minRPsPerStation; ///< minimal number of RPs associated with
                                         ///< a Cluster not to be rejected
                                         //
    bool verbose(unsigned int level);

    void verboseEnable(bool enable);
    bool verboseEnabled;


    /// Defines the readout direction
    enum Direction { U, V };

    /* SELECTION */

    /// Selection of \ref Cluster "Clusters" based on the number of RPs involved
    void selectClusters(vector<Cluster> &out, const vector<Cluster> &in);

    /* FILTERING */

    /// Filtering out \ref Cluster "Clusters" which happen to diverge from U/V patterns by more
    /// than \ref filter_dist_cut
    void filterClusters(
            vector<Cluster> &out,
            const vector<Cluster> &in,
            const RPRecognizedPatternsCollection &rpPatColl);

    /* OVERLAPS REMOVAL */

    /// Prepare the input for a \ref OverlapsRemoval algorithm and run it
    void removeOverlaps(
            vector<Cluster> &out,
            const vector<Cluster> &in,
            const RPRecognizedPatternsCollection &rpPatColl);

    // Cluster / RP

    /// \returns if the cluster is worthy of selection
    /// The choice is made based on the RPs accounting for the cluster.
    bool selectCluster(const Cluster &cluster);

    /// Checks if the track should leave a pattern in the RP
    bool trackHitsRP(Track track, RPId rpID);

    // Cluster / patterns

    /// Verifies if the track associated with cluster matches patterns in all RPs
    bool clusterMatchesPatterns(
            const Cluster &cluster,
            const RPRecognizedPatternsCollection &rpPatColl);

    /// Associates cluster with patterns in each plane
    vector<unsigned int> associateClusterWithHits(
            const Cluster &cluster,
            const RPRecognizedPatternsCollection &rpPatColl);

    /// Associates cluster with a pattern in a single plane
    pair<unsigned int, double> assocClusterWithPattern(
            const Cluster &cluster,
            RPId rpID,
            const vector<RPRecognizedPatterns::Line> &vPatterns,
            Direction direction);

    /// Associates U/V patterns with \ref Cluster "Clusters" based on distance
    /// between them.
    void associateClustersWithHits(
            vector<vector<unsigned int> > &cluHits,
            const vector<Cluster> &clusters,
            const RPRecognizedPatternsCollection &rpPatColl);

};

} // namespace RPStationMultiTrackFinderFitter

#endif
