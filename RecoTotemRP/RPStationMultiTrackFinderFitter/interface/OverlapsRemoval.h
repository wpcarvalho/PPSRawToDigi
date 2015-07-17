/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jakub Sawicki (jakub.sawicki@cern.ch, jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#ifndef RPStationMultiTrackFinderFitter_OverlapsRemoval
#define RPStationMultiTrackFinderFitter_OverlapsRemoval

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatternsCollection.h"

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/Cluster.h"

#include <vector>
#include <map>
#include <utility>

namespace RPStationMultiTrackFinderFitter {

/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief Removes non-compatible \ref Cluster "Clusters" (tracks)
 *
 * It performs the search of the best combinations of Clusters
 * within the full problem space. If more than single solution
 * is found the algorithm returns empty result not being
 * able to distinguish between the solutions.
 *
 * The algorithm is recursive.
 *
 * Solution is considered good if:
 * * each U/V pattern is associated with at least one chosen track
 * * number of overlaps (in terms of U/V patterns) between chosen tracks
 *   is acceptable
 */
class OverlapsRemoval
{
  public:
    explicit OverlapsRemoval(const edm::ParameterSet &ps);
    ~OverlapsRemoval();

    /// initiates the recursion
    void run(
            vector<Cluster> &out, ///< [out] output Clusters
            const vector<Cluster> &in, ///< [in] input Clusters
            const vector<unsigned int> &hitCount, ///< [in] number of U/V patterns in each plane (U and V in each RP)
            const vector<vector<unsigned int> > &cluHits ///< [in] for each Cluster a vector of associated U/V patterns
            );
  private:
    const unsigned int verbosity;

    /// Cluster -> U/V plane -> U/V pattern
    const vector<vector<unsigned int> > *cluHits;

    /// U/V plane -> number of hits
    const vector<unsigned int> *hitCountOriginal;

    /// U/V plane -> number of hits that may be satisfied by cluHits
    vector<unsigned int> hitCountPossible;

    /// keeps a vector of solutions
    vector<vector<bool> > v_chosen;

    /// Main step of the recursive search
    void goDown(
            vector<bool> &chosen, ///< [in] currently chosen set of tracks
            unsigned int pos      ///< [in] current position in the chosen vector
            );

    /// Determines the quality of the solution
    void feasibilityCheck(
            bool *feasible, ///< [out] solution is not contradictory
            bool *complete, ///< [out] solution is good
            const vector<bool> &chosen ///< [in]
            );

    /// if a U/V pattern cannot be satisfied with available tracks then ignore
    /// it during further analysis
    void calculatePossibleHitCount();
};

} // namespace RPStationMultiTrackFinderFitter

#endif // RPStationMultiTrackFinderFitter_OverlapsRemoval
