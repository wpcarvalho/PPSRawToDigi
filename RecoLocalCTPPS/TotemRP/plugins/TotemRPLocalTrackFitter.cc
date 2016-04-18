/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
* 	Hubert Niewiadomski
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"
#include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"

#include "RecoLocalCTPPS/TotemRP/interface/TotemRPLocalTrackFitterAlgorithm.h"

//----------------------------------------------------------------------------------------------------

/**
 *\brief TODO: Fits every track candidate in RPTrackCandidateCollection.
 **/
class TotemRPLocalTrackFitter : public edm::one::EDProducer<>
{
  public:
    explicit TotemRPLocalTrackFitter(const edm::ParameterSet& conf);
    virtual ~TotemRPLocalTrackFitter();
    virtual void beginJob();
    virtual void produce(edm::Event& e, const edm::EventSetup& c);

  private:
    int verbosity_;

    /// Selection of the pattern-recognition module.
    edm::InputTag tagUVPattern;

    edm::EDGetTokenT<edm::DetSetVector<TotemRPUVPattern>> patternCollectionToken;
    
    /// A watcher to detect geometry changes.
    edm::ESWatcher<VeryForwardRealGeometryRecord> geometryWatcher;
    
    /// The instance of the fitter module
    TotemRPLocalTrackFitterAlgorithm fitter_;
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

TotemRPLocalTrackFitter::TotemRPLocalTrackFitter(const edm::ParameterSet& conf)
   : verbosity_(conf.getParameter<int>("Verbosity")), fitter_(conf)
{
  tagUVPattern = conf.getParameter<edm::InputTag>("tagUVPattern");
  patternCollectionToken = consumes<DetSetVector<TotemRPUVPattern>>(tagUVPattern);

  produces<DetSetVector<TotemRPLocalTrack>>();
}

//----------------------------------------------------------------------------------------------------

TotemRPLocalTrackFitter::~TotemRPLocalTrackFitter()
{
}

//----------------------------------------------------------------------------------------------------

void TotemRPLocalTrackFitter::beginJob()
{
}

//----------------------------------------------------------------------------------------------------

void TotemRPLocalTrackFitter::produce(edm::Event& e, const edm::EventSetup& setup)
{
  // get geometry
  edm::ESHandle<TotemRPGeometry> geometry;
  setup.get<VeryForwardRealGeometryRecord>().get(geometry);

  if (geometryWatcher.check(setup))
    fitter_.Reset();
  
  // get input
  edm::Handle<DetSetVector<TotemRPUVPattern>> input;
  e.getByToken(patternCollectionToken, input);

  // run fit for each RP
  DetSetVector<TotemRPLocalTrack> output;
  
  for (const auto &rpv : *input)
  {
    det_id_type rpId =  rpv.detId();

    // is U-V association unique?
    unsigned int n_U=0, n_V=0;
    unsigned int idx_U=0, idx_V=0;
    for (unsigned int pi = 0; pi < rpv.size(); pi++)
    {
      const TotemRPUVPattern &pattern = rpv[pi];

      if (pattern.getFittable() == false)
        continue;

      switch (pattern.getProjection())
      {
        case TotemRPUVPattern::projU:
          n_U++;
          idx_U=pi;
          break;
     
        case TotemRPUVPattern::projV:
          n_V++;
          idx_V=pi;
          break;

        default:
          break;
      }   
    }

    if (n_U != 1 || n_V != 1)
    {
      if (verbosity_)
        printf(">> TotemRPLocalTrackFitter::produce > Impossible to combine U and V patterns in RP %u (n_U=%u, n_V=%u).\n",
          rpId, n_U, n_V);

      continue;
    }

    // combine U and V hits
    DetSetVector<TotemRPRecHit> hits;
    for (auto &ids : rpv[idx_U].getHits())
    {
      auto &ods = hits.find_or_insert(ids.detId());
      for (auto &h : ids)
        ods.push_back(h);
    }

    for (auto &ids : rpv[idx_V].getHits())
    {
      auto &ods = hits.find_or_insert(ids.detId());
      for (auto &h : ids)
        ods.push_back(h);
    }

    // run fit
    double z0 = geometry->GetRPGlobalTranslation(rpId).z();

    TotemRPLocalTrack track;
    fitter_.FitTrack(hits, z0, *geometry, track);
    
    if (track.isValid())
    {
      DetSet<TotemRPLocalTrack> ds = output.find_or_insert(rpId);
      ds.push_back(track);
    }
  }

  // save results
  e.put(make_unique<DetSetVector<TotemRPLocalTrack>>(output));
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(TotemRPLocalTrackFitter);
