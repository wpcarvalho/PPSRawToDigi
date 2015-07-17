#ifndef RecoTotemRP_RPStationMultiTrackFinderFitter_Fitter_h
#define RecoTotemRP_RPStationMultiTrackFinderFitter_Fitter_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPStationTrackFitCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatternsCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidateCollection.h"

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/Cluster.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/GeometryHelper.h"

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
 * \brief Refines the cluster fit
 *
 * It performs a least squared fitting to obtain the best
 * possible track fit.
 *
 * \author Jan Kaspar <jan.kaspar@gmail.com>
 * \author Jakub Sawicki <jakub.kuba.sawicki@gmail.com>
 */
class Fitter {
public:
    explicit Fitter(
        const edm::ParameterSet& ps,
        GeometryHelper &geometry);
    virtual ~Fitter();

    /**
     * \brief Performs the tracks' fitting
     */
    virtual void fitCollection(
            /// [out] collection of refined tracks
            RPStationTrackFitCollection *rpTrColl,
            /// [in] collection of Clusters
            const vector<Cluster> *input);

    /**
     * \brief Perform single track fitting
     */
    virtual void fitCollection(
            /// [out] collection of refined tracks
            RPStationTrackFitCollection *rpTrColl,
            /// [in] track candidate collection
            const RPTrackCandidateCollection *input);

    /**
     * \brief Performs the track fitting
     *
     * Using least squares method it determines the best possible track fit.
     */
    virtual void fit(
            /// [out] fitted track
            RPStationTrackFit &out,
            /// [in] the hits to be fitted
            const vector<RPRecoHit> &hits);

private:
    const unsigned int verbosity;
    GeometryHelper &geometry;
    const bool enableFitting;
    const double z_h;

    /**
     * \brief Reduces Cluster to track fit without actual fitting
     */
    void clusterToTrack(
            /// [out] fitted track
            RPStationTrackFit &out,
            /// [in] cluster to be reduced to track
            const Cluster &cluster);
};

} // namespace RPStationMultiTrackFinderFitter

#endif
