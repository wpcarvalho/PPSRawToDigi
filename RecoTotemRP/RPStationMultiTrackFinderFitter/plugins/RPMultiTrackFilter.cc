/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors:
*  Jan Kašpar (jan.kaspar@gmail.com)
*  Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
*
* $$Id: RPMultiTrackFilter.cc 8163 2013-09-16 10:10:04Z psikora $: $
* $Revision: 8163 $
* $Date: 2013-09-16 12:10:04 +0200 (pon, 16 wrz 2013) $
*
****************************************************************************/

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatternsCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPStationTrackFitCollection.h"


namespace RPStationMultiTrackFinderFitter {

/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief Filters tracks based on information about previous stages of reconstruction
 *
 * Available filters:
 * * (generated) track count
 * * reconstructible track count
 * * reconstructed track count
 * * minimal distance between U/V patterns
 * * minimal distance between tracks
 */
class RPMultiTrackFilter : public edm::EDFilter {
public:
    RPMultiTrackFilter(const edm::ParameterSet &);

protected:
    edm::InputTag generatorLabel,
                  oneRPRecoLabel,
                  stationRecoLabel;

    /// Keeps the information if the filter is active and its limits
    template<typename T>
    struct Criterion {
      bool active;
      T min, max;
    };

    Criterion<unsigned int> trackCount,
                            reconstructableTrackCount,
                            recoTrackCount;

    Criterion<double> minDistBetweenPatterns,
                      minDistBetweenTracks;


    /// Excludes events not meeting all the active criterions
    virtual bool filter(edm::Event&, const edm::EventSetup &);

    /// Converts the ParameterSet into a proper criterion
    template<typename T>
    Criterion<T> PSToCriterion(const edm::ParameterSet &ps);

    /// Checks if the value is inside [min;max)
    template<typename T>
    inline bool inside(T value, T min, T max) {
        return min <= value && value < max;
    }
};


using namespace std;


//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

RPMultiTrackFilter::RPMultiTrackFilter(const edm::ParameterSet &ps) :
    generatorLabel(ps.getParameter<edm::InputTag>("generatorLabel")),
    oneRPRecoLabel(ps.getParameter<edm::InputTag>("oneRPRecoLabel")),
    stationRecoLabel(ps.getParameter<edm::InputTag>("stationRecoLabel")),
    trackCount(
        PSToCriterion<unsigned int>(
            ps.getParameterSet("trackCount"))),
    reconstructableTrackCount(
        PSToCriterion<unsigned int>(
            ps.getParameterSet("reconstructableTrackCount"))),
    recoTrackCount(
        PSToCriterion<unsigned int>(
            ps.getParameterSet("recoTrackCount"))),
    minDistBetweenPatterns(
        PSToCriterion<double>(
            ps.getParameterSet("minDistBetweenPatterns"))),
    minDistBetweenTracks(
        PSToCriterion<double>(
            ps.getParameterSet("minDistBetweenTracks"))) { }

//----------------------------------------------------------------------------------------------------

template<typename T>
RPMultiTrackFilter::Criterion<T> RPMultiTrackFilter::PSToCriterion(const edm::ParameterSet &ps) {
    Criterion<T> crit;
    crit.active = ps.getParameter<bool>("active");
    crit.min    = ps.getParameter<T>("min");
    crit.max    = ps.getParameter<T>("max");
    return crit;
}

bool RPMultiTrackFilter::filter(edm::Event &event, const edm::EventSetup &es)
{
    if (trackCount.active) {
        edm::Handle<RPStationTrackFitCollection> stationSimu;
        event.getByLabel(generatorLabel, stationSimu);

        unsigned int total = 0;
        for (RPStationTrackFitCollection::const_iterator armIt = stationSimu->begin();
                armIt != stationSimu->end(); armIt++) {
            total += (unsigned int) armIt->second.size();
        }

        if (!inside<unsigned int>(total, trackCount.min, trackCount.max)) {
            return false;
        }
    }

    if (reconstructableTrackCount.active) {
        edm::Handle<RPStationTrackFitCollection> stationSimu;
        event.getByLabel(generatorLabel, stationSimu);

        unsigned int total = 0;
        for (RPStationTrackFitCollection::const_iterator armIt = stationSimu->begin();
                armIt != stationSimu->end(); armIt++) {
            for (vector<RPStationTrackFit>::const_iterator trackIt = armIt->second.begin();
                    trackIt != armIt->second.end(); trackIt++) {
                if (trackIt->valid) {
                    total++;
                }
            }
        }

        if (!inside<unsigned int>(total, reconstructableTrackCount.min, reconstructableTrackCount.max)) {
            return false;
        }
    }

    if (recoTrackCount.active) {
        edm::Handle<RPStationTrackFitCollection> stationReco;
        event.getByLabel(stationRecoLabel, stationReco);

        unsigned int total = 0;
        for (RPStationTrackFitCollection::const_iterator armIt = stationReco->begin();
                armIt != stationReco->end(); armIt++) {
            total += (unsigned int) armIt->second.size();
        }

        if (!inside<unsigned int>(total, recoTrackCount.min, recoTrackCount.max)) {
            return false;
        }
    }

    if (minDistBetweenPatterns.active) {
        edm::Handle<RPRecognizedPatternsCollection> patterns;
        event.getByLabel(oneRPRecoLabel, patterns);

        double minDist = 1e100;

        for (RPRecognizedPatternsCollection::const_iterator rpIt = patterns->begin();
                rpIt != patterns->end(); rpIt++) {
            if (rpIt->second.uLines.size() > 1) {
                for (vector<RPRecognizedPatterns::Line>::const_iterator lineIt = rpIt->second.uLines.begin();
                        lineIt != rpIt->second.uLines.end(); lineIt++) {
                    for (vector<RPRecognizedPatterns::Line>::const_iterator lineIt2 = lineIt + 1;
                            lineIt2 != rpIt->second.uLines.end(); lineIt2++) {
                        double dist = fabs(lineIt->b - lineIt2->b);
                        if (dist < minDist) {
                            minDist = dist;
                        }
                    }
                }
            }

            if (rpIt->second.vLines.size() > 1) {
                for (vector<RPRecognizedPatterns::Line>::const_iterator lineIt = rpIt->second.vLines.begin();
                        lineIt != rpIt->second.vLines.end(); lineIt++) {
                    for (vector<RPRecognizedPatterns::Line>::const_iterator lineIt2 = lineIt + 1;
                            lineIt2 != rpIt->second.vLines.end(); lineIt2++) {
                        double dist = fabs(lineIt->b - lineIt2->b);
                        if (dist < minDist) {
                            minDist = dist;
                        }
                    }
                }
            }
        }

        if (!inside(minDist, minDistBetweenPatterns.min, minDistBetweenPatterns.max)) {
            return false;
        }
    }

    if (minDistBetweenTracks.active) {
        edm::Handle<RPStationTrackFitCollection> stationSimu;
        event.getByLabel(generatorLabel, stationSimu);

        double minDist = 1e100;

        for (RPStationTrackFitCollection::const_iterator armIt = stationSimu->begin();
                armIt != stationSimu->end(); armIt++) {
            for (vector<RPStationTrackFit>::const_iterator track1It = armIt->second.begin();
                    track1It != armIt->second.end(); track1It++) {
                for (vector<RPStationTrackFit>::const_iterator track2It = track1It + 1;
                        track2It != armIt->second.end(); track2It++) {
                    const RPStationTrackFit &t_r = *track1It,
                                            &t_s = *track2It;

                    double si_de_ax = 3.,
                           si_de_ay = 3.,
                           si_de_bx = 0.02,
                           si_de_by = 0.02;

                    double de_ax = (t_r.ax - t_s.ax)*1E6;
                    double de_bx = (t_r.bx - t_s.bx - t_s.ax*(t_r.z0 - t_s.z0));
                    double de_ay = (t_r.ay - t_s.ay)*1E6;
                    double de_by = (t_r.by - t_s.by - t_s.ay*(t_r.z0 - t_s.z0));

                    double dist_de_ax = de_ax / si_de_ax;
                    double dist_de_ay = de_ay / si_de_ay;
                    double dist_de_bx = de_bx / si_de_bx;
                    double dist_de_by = de_by / si_de_by;

                    double dist_sq = dist_de_ax*dist_de_ax + dist_de_ay*dist_de_ay + dist_de_bx*dist_de_bx + dist_de_by*dist_de_by;

                    if (dist_sq < minDist) {
                        minDist = dist_sq;
                    }
                }
            }
        }

        if (!inside(minDist, minDistBetweenTracks.min, minDistBetweenTracks.max)) {
            return false;
        }
    }

    return true;
}


DEFINE_FWK_MODULE(RPMultiTrackFilter);

} // namespace RPStationMultiTrackFinderFitter
