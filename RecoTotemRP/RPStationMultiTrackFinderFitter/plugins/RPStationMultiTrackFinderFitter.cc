/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jan Kašpar (jan.kaspar@gmail.com)
 *   Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

/**\class RPStationMultiTrackFinderFitter RPStationMultiTrackFinderFitter.cc RecoTotemRP/RPStationMultiTrackFinderFitter/plugins/RPStationMultiTrackFinderFitter.cc

 Description: It performs multi track reconstruction of events

 Implementation:

 Iterate through all the combinations of pots and create track candidates
 using the Hough transform.

 Match track candidates with hits in U and V planes and use a simple algorithm
 to reject tracks that create too many conflicts with other tracks.

 The additional fitting is performed to improve the track fit.
*/


#include <memory>

#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/RPStationMultiTrackFinderFitter.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/Algorithm.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/ClusterBuilderSubsets.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPStationTrackFitCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatternsCollection.h"

namespace RPStationMultiTrackFinderFitter {

RPStationMultiTrackFinderFitter::RPStationMultiTrackFinderFitter(const edm::ParameterSet& ps) :
    verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),
    geometry(ps),
    clusterBuilder(new ClusterBuilderSubsets(ps, geometry)),
    recoAlgorithm(ps, geometry, *clusterBuilder),
    fitter(ps, geometry),
    RPTrackCandCollProducer(ps.getParameter<string>("RPTrackCandCollProducer"))
{
    // register products
    produces<RPStationTrackFitCollection>();
}

//----------------------------------------------------------------------------------------------------

RPStationMultiTrackFinderFitter::~RPStationMultiTrackFinderFitter()
{
    delete clusterBuilder;
}

//----------------------------------------------------------------------------------------------------

void RPStationMultiTrackFinderFitter::produce(edm::Event& event, const edm::EventSetup& es)
{
    using namespace edm;

    if (verbosity > 1)
        printf("\n\n----------------------------------------------------------------------------------------------------\n");

    if (verbosity > 0)
        printf(">> RPStationMultiTrackFinderFitter::produce (%u:%u)\n", event.id().run(), event.id().event());

    // get geometry
    if (watcherRealGeometry.check(es))
    {
        ESHandle<TotemRPGeometry> newGeometry;
        es.get<RealGeometryRecord>().get(newGeometry);
        geometry.setGeometry(&(*newGeometry));
    }

    // get input
    Handle<RPRecognizedPatternsCollection> input;
    event.getByLabel(RPTrackCandCollProducer, input);

    // prepare output
    auto_ptr<RPStationTrackFitCollection> output(new RPStationTrackFitCollection());

    // recognize linear patterns
    vector<Cluster> clusters;
    recoAlgorithm.reconstruct(&clusters, input.product());

    // refit the recognized patterns
    fitter.fitCollection(output.get(), &clusters);

    // save output
    event.put(output);
}

DEFINE_FWK_MODULE(RPStationMultiTrackFinderFitter);

} // namespace RPStationMultiTrackFinderFitter
