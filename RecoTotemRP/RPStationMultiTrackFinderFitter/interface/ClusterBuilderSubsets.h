/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#ifndef RPStationMultiTrackFinderFitter_ClusterBuilderSubsets
#define RPStationMultiTrackFinderFitter_ClusterBuilderSubsets

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatternsCollection.h"
#include "DataFormats/TotemRPDataTypes/interface/RPRecoHit.h"

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/GeometryHelper.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/IClusterBuilder.h"

#include "TMatrixD.h"
#include "TVectorD.h"

#include <iterator>
#include <set>

namespace RPStationMultiTrackFinderFitter {

/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief Uses graph analysis to determine the best \ref Cluster choice
 */
class ClusterBuilderSubsets : public IClusterBuilder
{
  public:
    explicit ClusterBuilderSubsets(const edm::ParameterSet &ps, GeometryHelper &geometry);

    /**
     * \brief Manages the algorithm steps
     *
     * The graph is first created using populateGraph then
     * the \ref Cluster "Clusters" are created from graph subsets
     * returned from findLargestNonInclusiveCliques.
     */
    virtual void build(
            vector<Cluster> &clusters,
            const RPRecognizedPatternsCollection &rpPatColl);
  private:
    /// cut offs to reject apparently wrong track candidates
    const double ax_mean, ax_cut;
    const double ay_mean, ay_cut;

    /// the projection plane position
    const double z_h;

    /// assumed U/V pattern uncertainty
    const double uv_unc;

    /// threshold of distance between track candidates to create
    /// an edge between them
    const double edge_create_th;

    /// minimal number of RPs involved to consider candidate interesting
    const unsigned int minRPsPerStation;

    typedef unsigned int CandidateId;
    typedef RPRecognizedPatterns::Line Line;

    /**
     * \brief Edge representation in the graph of track candidates
     * Each edge is undirected, so Edge(a, b) == Edge(b, a).
     */
    struct Edge
    {
        CandidateId first;
        CandidateId second;

        Edge(CandidateId first, CandidateId second) : first(first), second(second) {}

        bool operator<(const Edge &other) const
        {
            unsigned int max_me = max(first, second),
                         min_me = min(first, second),
                         max_other = max(other.first, other.second),
                         min_other = min(other.first, other.second);

            if (min_me == min_other)
            {
                return max_me < max_other;
            } else {
                return min_me < min_other;
            }
        }

        bool operator==(const Edge &other) const
        {
            unsigned int max_me = max(first, second),
                         min_me = min(first, second),
                         max_other = max(other.first, other.second),
                         min_other = min(other.first, other.second);

            return (max_me == max_other && min_me == min_other);
        }
    };

    /**
     * \brief A track candidate
     * Consists of th keeping the parameters of the track (ax bx ay by)
     * and its uncertainty matrices.
     */
    struct Candidate
    {
        TVectorD th;
        TMatrixD V_th, V_th_i;
        RPId rpId1, rpId2;
        set<RPRecoHit> hits;
        Candidate() : th(4), V_th(4, 4), V_th_i(4, 4) {}

        /// \returns the squared distance between the \ref Candidate and the other one
        double dist(const Candidate &other);
    };

    /**
     * \brief Creates the graph of \ref Candidate "Candidates" with compatible ones connected
     * It iterates over all possible two RP combinations.
     * Uses buildCandidatesFrom2RPs to generate all candidates.
     */
    void populateGraph(
            vector<Candidate> &candidates, ///< [out] vertices
            set<Edge> &edges, ///< [out] edges between compatible vertices
            const RPRecognizedPatternsCollection &rpPatColl ///< [in]
            );

    /**
     * \brief Generates all sensible track candidates from combination of 2 RPs
     * Iterates over all combinations of U/V patterns to create Candidates.
     * Once each candidate is constructed it is validated.
     */
    vector<Candidate> buildCandidatesFrom2RPs(
        RPId rpID1, RPId rpID2,
        const RPRecognizedPatterns &rpPat1, const RPRecognizedPatterns &rpPat2);

    /**
     * \brief Creates a track candidate using data from U/V patterns
     * A Candidate is created. It crosses through the specified U/V patterns.
     */
    Candidate buildCandidate(RPId rpID1, RPId rpID2,
            const Line &u1, const Line &v1,
            const Line &u2, const Line &v2);

    /**
     * \brief Checks if the Candidate is to be considered further
     * If Candidate is too inclined or not passing through both specified
     * RPs (due to the edge cut) returns false.
     */
    bool candidateOk(RPId rpID1, RPId rpID2, const Candidate &candidate);

    /**
     * \brief Generates the global coordinates of a Candidate hit with specified RP
     */
    CLHEP::Hep3Vector coordsAtRP(RPId rpID, const Candidate &candidate);

    /**
     * \brief Finds the biggest set of largest non-inclusive cliques
     * Searches the graph for cliques which are mutually non-inclusive.
     * A search of a full problem space is performed. Due to this fact
     * the algorithm becomes highly time consuming for number of RPs
     * active over 4 (large cliques are hindering the performance).
     */
    void findLargestNonInclusiveCliques(
            set<set<CandidateId> > &subsets, ///< [out]
            const vector<Candidate> &candidates, ///< [in]
            const set<Edge> &edges ///< [in]
            );

    /// Checks if all members of cand_v are connected to cand_id.
    /// Cand_id is considered not to be connected to itself.
    /// \returns if cand_id connected to all members of cand_v
    bool allCandidateIdsConnectedTo(const vector<CandidateId> &cand_v,
            CandidateId cand_id, const set<Edge> &edges);

    /// Converts the subsets into \ref Cluster "Clusters"
    void mergeSubsetsIntoClusters(vector<Cluster> &clusters, ///< [out]
            const set<set<CandidateId> > &subsets, ///< [in]
            const vector<Candidate> &candidates ///< [in]
            );
};

} // namespace RPStationMultiTrackFinderFitter

#endif // #ifndef RPStationMultiTrackFinderFitter_ClusterBuilderSubsets
