/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "DataFormats/TotemRPDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"

#include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"

#include "RecoTotemRP/RPNonParallelTrackCandidateFinder/interface/FastLineRecognition.h"

/**
 * \brief Class to recognize straight line tracks, based on optimized Hough trasform.
 *
 * The search is perfomed in global U,V coordinates (wrt. beam). In this way (some of)
 * the alignment corrections can be taken into account.
**/
class RPNonParallelTrackCandidateFinder : public edm::EDProducer
{
  public:
    RPNonParallelTrackCandidateFinder(const edm::ParameterSet& conf);
    virtual ~RPNonParallelTrackCandidateFinder();
    virtual void produce(edm::Event& e, const edm::EventSetup& c);
  
  private:
    edm::InputTag detSetVectorTotemRPRecHitLabel;
    edm::EDGetTokenT<edm::DetSetVector<TotemRPRecHit> > detSetVectorTotemRPRecHitToken;

    unsigned int verbosity;

    /// minimal required number of active planes per projection to even start track recognition
    unsigned char minPlanesPerProjectionToSearch;

    /// minimal required number of active planes per projection to mark track candidate as fittable
    unsigned char minPlanesPerProjectionToFit;

    /// above this limit, planes are considered noisy
    unsigned int maxHitsPerPlaneToSearch;

    /// the line recognition algorithm
    FastLineRecognition *lrcgn;

    /// minimal weight of (Hough) cluster to accept it as candidate
    double threshold;

    /// maximal angle (in any projection) to mark candidate as fittable - controls track parallelity
    double max_a_toFit;

    /// whether to allow to combine most significant U and V pattern, in case there several of them
    bool allowAmbiguousCombination;

    /// exceptional settings, per RP and per projection
    std::vector<edm::ParameterSet> exceptionalSettings;

    edm::ESWatcher<VeryForwardRealGeometryRecord> geometryWatcher;

    /// executes line recognition in a projection
    void RecognizeAndSelect(TotemRPUVPattern::ProjectionType proj, double z0, double threshold,
      unsigned int planes_required,
      const std::vector<const TotemRPRecHit *> &hits, edm::DetSet<TotemRPUVPattern> &patterns);
};

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

RPNonParallelTrackCandidateFinder::RPNonParallelTrackCandidateFinder(const edm::ParameterSet& conf) :
  detSetVectorTotemRPRecHitLabel(conf.getParameter<edm::InputTag>("DetSetVectorTotemRPRecHitLabel")),
  verbosity(conf.getUntrackedParameter<unsigned int>("verbosity", 0)),
  minPlanesPerProjectionToSearch(conf.getParameter<unsigned int>("minPlanesPerProjectionToSearch")),
  minPlanesPerProjectionToFit(conf.getParameter<unsigned int>("minPlanesPerProjectionToFit")),
  maxHitsPerPlaneToSearch(conf.getParameter<unsigned int>("maxHitsPerPlaneToSearch")),
  lrcgn(new FastLineRecognition(conf.getParameter<double>("clusterSize_a"), conf.getParameter<double>("clusterSize_b"))),
  threshold(conf.getParameter<double>("threshold")),
  max_a_toFit(conf.getParameter<double>("max_a_toFit")),
  allowAmbiguousCombination(conf.getParameter<bool>("allowAmbiguousCombination")),
  exceptionalSettings(conf.getParameter< vector<ParameterSet> >("exceptionalSettings"))
{
  detSetVectorTotemRPRecHitToken = consumes<edm::DetSetVector<TotemRPRecHit> >(detSetVectorTotemRPRecHitLabel);

  produces<DetSetVector<TotemRPUVPattern>>();
}

//----------------------------------------------------------------------------------------------------

RPNonParallelTrackCandidateFinder::~RPNonParallelTrackCandidateFinder()
{
  delete lrcgn;
}

//----------------------------------------------------------------------------------------------------

void RPNonParallelTrackCandidateFinder::RecognizeAndSelect(TotemRPUVPattern::ProjectionType proj,
    double z0, double threshold_loc, unsigned int planes_required,
    const vector<const TotemRPRecHit *> &hits, DetSet<TotemRPUVPattern> &patterns)
{
  // run recognition
  lrcgn->GetPatterns(hits, z0, threshold_loc, patterns);
  
  // set pattern properties
  for (auto &p : patterns)
  {
    p.setProjection(proj);

    p.setFittable(true);

    set<unsigned int> planes;
    for (const auto &h : p.getHits())
      planes.insert(h.DetId() % 10);

    if (planes.size() < planes_required)
      p.setFittable(false);
    
    if (fabs(p.getA()) > max_a_toFit)
      p.setFittable(false);
  }

  if (verbosity > 5)
  {
    for (const auto &p : patterns)
      printf("\t\t\ta = %.3f, b = %.3f, w = %.3f\n", p.getA(), p.getB(), p.getW());
  }
}

//----------------------------------------------------------------------------------------------------

void RPNonParallelTrackCandidateFinder::produce(edm::Event& event, const edm::EventSetup& es)
{
  if (verbosity > 5)
    printf(">> RPNonParallelTrackCandidateFinder::produce (%u:%llu)\n", event.id().run(), event.id().event());

  // geometry
  ESHandle<TotemRPGeometry> geometry;
  es.get<VeryForwardRealGeometryRecord>().get(geometry);
  if (geometryWatcher.check(es))
    lrcgn->ResetGeometry(geometry.product());
  
  // get input
  edm::Handle< edm::DetSetVector<TotemRPRecHit> > input;
  event.getByToken(detSetVectorTotemRPRecHitToken, input);

  // prepare output
  DetSetVector<TotemRPUVPattern> patternsVector;
  
  // map hits to RP ids
  map< unsigned int, pair< vector<const TotemRPRecHit*>, vector<const TotemRPRecHit*> > > hitMap;
  map< unsigned int, pair< map<unsigned char, unsigned int>, map<unsigned char, unsigned int> > > planeOccupancyMap;

  for (DetSetVector<TotemRPRecHit>::const_iterator dit = input->begin(); dit != input->end(); ++dit)
  {
    unsigned int detId = TotemRPDetId::rawToDecId(dit->detId());
    unsigned int RPId = TotemRPDetId::rpOfDet(detId);
    unsigned int plane = detId % 10;
    bool uDir = TotemRPDetId::isStripsCoordinateUDirection(detId);

    for (DetSet<TotemRPRecHit>::const_iterator hit = dit->begin(); hit != dit->end(); ++hit)
    {
      if (uDir)
      {
        hitMap[RPId].first.push_back(& (*hit));
        planeOccupancyMap[RPId].first[plane]++;
      } else {
        hitMap[RPId].second.push_back(& (*hit));
        planeOccupancyMap[RPId].second[plane]++;
      }
    }
  }

  // track recognition pot by pot
  for (map< unsigned int, pair< vector<const TotemRPRecHit*>, vector<const TotemRPRecHit*> > >::iterator it = hitMap.begin();
      it != hitMap.end(); ++it)
  {
    unsigned int rpId = it->first;

    // merge default and exceptional settings (if available)
    unsigned int minPlanesPerProjectionToFit_U = minPlanesPerProjectionToFit;
    unsigned int minPlanesPerProjectionToFit_V = minPlanesPerProjectionToFit;
    double threshold_U = threshold;
    double threshold_V = threshold;

    for (auto ps : exceptionalSettings)
    {
      unsigned int setId = ps.getParameter<unsigned int>("rpId");
      if (setId == rpId)
      {
        minPlanesPerProjectionToFit_U = ps.getParameter<unsigned int>("minPlanesPerProjectionToFit_U");
        minPlanesPerProjectionToFit_V = ps.getParameter<unsigned int>("minPlanesPerProjectionToFit_V");
        threshold_U = ps.getParameter<double>("threshold_U");
        threshold_V = ps.getParameter<double>("threshold_V");
      }
    }

    const map<unsigned char, unsigned int> &uColl = planeOccupancyMap[rpId].first;
    const map<unsigned char, unsigned int> &vColl = planeOccupancyMap[rpId].second;

    if (verbosity > 5)
    {
      printf("\tRP %u\n", rpId);
      printf("\t\tall hits: u = %lu, v = %lu\n", it->second.first.size(), it->second.second.size());
      printf("\t\tall planes: u = %lu, v = %lu\n", uColl.size(), vColl.size());
    }

    // count planes with clean data (no showers, noise, ...)
    unsigned int uPlanes = 0, vPlanes = 0;
    for (map<unsigned char, unsigned int>::const_iterator pit = uColl.begin(); pit != uColl.end(); ++pit)
      if (pit->second <= maxHitsPerPlaneToSearch)
        uPlanes++;
    for (map<unsigned char, unsigned int>::const_iterator pit = vColl.begin(); pit != vColl.end(); ++pit)
      if (pit->second <= maxHitsPerPlaneToSearch)
        vPlanes++;

    if (verbosity > 5)
      printf("\t\tplanes with clean data: u = %u, v = %i\n", uPlanes, vPlanes);

    // discard RPs with too few reasonable planes
    if (uPlanes < minPlanesPerProjectionToSearch || vPlanes < minPlanesPerProjectionToSearch)
      continue;

    // prepare data containers
    DetSet<TotemRPUVPattern> &patterns = patternsVector.find_or_insert(rpId);

    // "typical" z0 for the RP
    double z0 = geometry->GetRPDevice(rpId)->translation().z();

    // u then v recognition
    if (verbosity > 5)
      printf("\t\tu recognition\n");
    RecognizeAndSelect(TotemRPUVPattern::projU, z0, threshold_U, minPlanesPerProjectionToFit_U, it->second.first, patterns);

    if (verbosity > 5)
      printf("\t\tv recognition\n");
    RecognizeAndSelect(TotemRPUVPattern::projU, z0, threshold_V, minPlanesPerProjectionToFit_V, it->second.second, patterns);
  }
 
  // save output
  event.put(make_unique<DetSetVector<TotemRPUVPattern>>(patternsVector));
}
 
//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(RPNonParallelTrackCandidateFinder);
