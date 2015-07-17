/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jan Kašpar (jan.kaspar@gmail.com)
 *   Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#ifndef RecoTotemRP_RPStationMultiTrackFinderFitter_RPStationMultiTrackFinderFitter_h
#define RecoTotemRP_RPStationMultiTrackFinderFitter_RPStationMultiTrackFinderFitter_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/Algorithm.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/Fitter.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/GeometryHelper.h"

#include <string>

/**
 * \defgroup RPStationMultiTrackFinderFitter
 * \brief Multi-track station reconstruction algorithm
 */
namespace RPStationMultiTrackFinderFitter {

/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief Provides the interface for CMSSW edm::EDProducer
 *
 * \author Jakub Sawicki <jakub.kuba.sawicki@gmail.com>
 */
class RPStationMultiTrackFinderFitter : public edm::EDProducer
{
  public:
    explicit RPStationMultiTrackFinderFitter(const edm::ParameterSet&);
    ~RPStationMultiTrackFinderFitter();
    
    /**
     * \brief Takes care of communicating with CMSSW
     * Retrieves the input data from edm::Event, initializes the output data
     * structures, runs the algorithm and returns the results to the event.
     */
    virtual void produce(edm::Event&, const edm::EventSetup&) override;

  private:
    unsigned int verbosity;
    /// the geometry helper used
    GeometryHelper geometry;

    /// cluster builder to be used
    IClusterBuilder *clusterBuilder;

    /// reconstruction algorithm used
    Algorithm recoAlgorithm;

    /// track fitter
    Fitter fitter;

    /// ID of a module performing 1-RP pattern reconstruction
    string RPTrackCandCollProducer;

    /// watcher for the changes in the geometry
    edm::ESWatcher<RealGeometryRecord> watcherRealGeometry;
};

} // namespace RPStationMultiTrackFinderFitter

#endif
