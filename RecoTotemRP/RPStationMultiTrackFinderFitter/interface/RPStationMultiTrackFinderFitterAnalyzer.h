/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jakub Sawicki (jakub.sawicki@cern.ch, jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPStationTrackFitCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatternsCollection.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/GeometryHelper.h"

#include "TFile.h"
#include "TH2D.h"

#include <cstdio>

using namespace std;

namespace RPStationMultiTrackFinderFitter {

/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief Compares the output of the reconstruction with simulation
 *
 * \author Jakub Sawicki <jakub.kuba.sawicki@gmail.com>
 */
class RPStationMultiTrackFinderFitterAnalyzer : public edm::EDAnalyzer {
public:

    /**
     * \brief Opens/Creates the files, initializes empty structures
     */
    RPStationMultiTrackFinderFitterAnalyzer(const edm::ParameterSet &ps);
    /**
     * \brief Saves the summary, closes all the open file handles
     */
    virtual ~RPStationMultiTrackFinderFitterAnalyzer();
    /**
     * \brief Deals with the CMSSW and calls the internal methods doing the analysis
     */
    virtual void analyze(const edm::Event &e, const edm::EventSetup &es);

private:
    GeometryHelper geometry;
    edm::ESWatcher<RealGeometryRecord> watcherRealGeometry;

    const unsigned int verbosity;

    /// name of the module which is doing the simulation
    const string simTracksProducer;
    
    /// name of the module which is doing the reconstruction
    const string recoTracksProducer;
    
    /// name of the module which is doing the 1-RP reconstruction
    const string oneRPPatternsProducer;

    /// name of a file where statistics will be saved, if empty no statistics are saved
    const string outputFileName;
    
    /// name of a file where histograms will be kept, if empty no histograms are saved
    const string rootFileName;
    
    /// upper bound on the count of missing tracks in a histogram
    const unsigned int missing_max;
    
    /// upper bound on the count of fake tracks in a histogram
    const unsigned int fake_max;

    /// one sigma distance in each direction
    const double si_de_ax; // urad
    const double si_de_ay; // urad
    const double si_de_bx; // mm
    const double si_de_by; // mm
    /// limit of squared distance between simulated and reconstructed tracks to be
    /// considered compatible in units of sigma
    const double acceptance;


    FILE *f_out; ///< pointer to output file
    const bool f_out_enable; ///< whether to enable output to file

    const bool tf_out_enable; ///< whether to enable output to ROOT file
    TH2D *missing_fake; ///< histogram of missing vs fake cases

    TH1D *h_min_pat_dist_all;         ///< histogram of minimal distance between U/V patterns
    TH1D *h_min_pat_dist_fail;        ///< histogram of minimal distance between U/V patterns in failing cases
    TH1D *h_min_pat_dist_fake_comb;   ///< histogram of minimal distance between U/V patterns for events where the nubmer of reasonable fits is larger than the number of true tracks

    TH1D *track_distances; ///< histogram of distances between matching tracks
    TH2D *dist_x;
    TH2D *dist_y;
    TH2D *fail_dist_x;
    TH2D *fail_dist_y;
    TH2D *dist_u;
    TH2D *dist_v;
    TH2D *fail_dist_u;
    TH2D *fail_dist_v;

    TH1D *h_num_reas_fits;

    /// keeping the statistics
    struct {
        unsigned int total,       ///< total number of qualified events
                     failed,      ///< number of cases where reconstruction failed
                     fake,        ///< number of cases where fake tracks were detected
                     missing,     ///< number of cases where missing tracks were detected
                     fakemissing; ///< number of cases where both fake and missing tracks were detected
    } stats;

    struct missing_fake_pair {
        unsigned int missing, fake;
    };

    /**
     * \brief Compare set of simulated and reconstructed tracks with each other
     *
     * For each simulated track it looks for a corresponding reconstructed track
     * while keeping the statistics about missing and fake. It also deals with
     * a track being assigned to multiple simulated tracks, to count it only once.
     *
     * To determine the distance between the tracks an euclidean metric is used
     * where distances are normalized by si_de_[ab][xy] properties. Threshold
     * of \ref acceptance is used.
     */
    virtual missing_fake_pair AnalyzeTracks(
            const edm::Event &event,                       ///< [in] the event object
            const vector<RPStationTrackFit> &tracks_simu,  ///< [in] the vector of simulated tracks
            const vector<RPStationTrackFit> &tracks_reco); ///< [in] tracks_reco the vector of reconstructed tracks

    /**
     * \brief Generates statistics of distances between U/V patterns
     *
     * Looks for the smallest distance between the U/V patterns provided and
     * saves the output in the h_min_pat_dist_all histogram.
     *
     * If reconstruction failed (based on miss_fake) the point is also inserted
     * into h_min_pat_dist_fail.
     */
    virtual void AnalyzePatterns(
            const RPRecognizedPatternsCollection &patterns, ///< [in] collection of U/V patterns
            const missing_fake_pair &miss_fake, ///< [in] number of failures
            bool fake_comb);

    /**
     * \brief Generate statistics of difference between simulated tracks
     */
    virtual void AnalyzeSimuTracks(
        const vector<RPStationTrackFit> &tracks_simu,
        const missing_fake_pair &miss_fake);
    
    // TODO
    unsigned int AnalyzeReasonableInputFits(const RPRecognizedPatternsCollection &patterns);

    /**
     * \brief Outputs the overall statistics to the specified files
     */
    virtual void WriteSummary();
};

} // namespace RPStationMultiTrackFinderFitter

