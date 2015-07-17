/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jan Kašpar (jan.kaspar@gmail.com)
 *   Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/


#include <memory>

#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/RPStationFitter.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/Algorithm.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/ClusterBuilderIterable.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/ClusterBuilderSubsets.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPStationTrackFitCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatternsCollection.h"

namespace RPStationMultiTrackFinderFitter {

/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief Fits the sensor hits from RPTrackCandidateCollection
 */
RPStationFitter::RPStationFitter(const edm::ParameterSet& ps) :
    verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),
    geometry(ps),
    fitter(ps, geometry),
    RPTrackCandCollProducer(ps.getParameter<string>("RPTrackCandCollProducer"))
{
    // register products
    produces<RPStationTrackFitCollection>();
}

RPStationFitter::~RPStationFitter() { }

//
// member functions
//

void RPStationFitter::produce(edm::Event& event, const edm::EventSetup& es) {
    using namespace edm;

    if (verbosity > 0) {
        printf(">> RPStationFitter::produce (%u:%u)\n", event.id().run(), event.id().event());
    }

    // get geometry
    if (watcherRealGeometry.check(es)) {
        ESHandle<TotemRPGeometry> newGeometry;
        es.get<RealGeometryRecord>().get(newGeometry);
        geometry.setGeometry(&(*newGeometry));
    }

    // get input
    Handle<RPTrackCandidateCollection> input;
    event.getByLabel(RPTrackCandCollProducer, input);

    //prepare output
    auto_ptr<RPStationTrackFitCollection> output(new RPStationTrackFitCollection());

    fitter.fitCollection(output.get(), &*input);

    // save output
    event.put(output);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RPStationFitter);
} // namespace RPStationMultiTrackFinderFitter

