/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
* 	Hubert Niewiadomski
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

// TODO: move the header file here
#include "RecoTotemRP/RPTrackCandidateCollectionFitter/interface/RPTrackCandidateCollectionFitter.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"

// TODO: needed ?
#include <iostream>

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

RPTrackCandidateCollectionFitter::RPTrackCandidateCollectionFitter(const edm::ParameterSet& conf)
   : verbosity_(conf.getParameter<int>("Verbosity")), fitter_(conf)
{
  patternCollectionLabel = conf.getParameter<edm::InputTag>("RPTrackCandidateCollectionLabel");
  patternCollectionToken = consumes<DetSetVector<TotemRPUVPattern>>(patternCollectionLabel);

  produces<DetSetVector<TotemRPLocalTrack>> ();
}

//----------------------------------------------------------------------------------------------------

RPTrackCandidateCollectionFitter::~RPTrackCandidateCollectionFitter()
{
}

//----------------------------------------------------------------------------------------------------

void RPTrackCandidateCollectionFitter::beginJob()
{
}

//----------------------------------------------------------------------------------------------------

void RPTrackCandidateCollectionFitter::produce(edm::Event& e, const edm::EventSetup& setup)
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
        printf(">> RPTrackCandidateCollectionFitter::produce > Impossible to combine U and V patterns in RP %u (n_U=%u, n_V=%u).\n",
          rpId, n_U, n_V);

      continue;
    }

    // build hit collection
    vector<const TotemRPRecHit *> hits;
    for (auto &h : rpv[idx_U].getHits())
      hits.push_back(&h);
    for (auto &h : rpv[idx_V].getHits())
      hits.push_back(&h);

    // run fit
    double z0 = geometry->GetRPGlobalTranslation(rpId).z();

    TotemRPLocalTrack track;
    fitter_.FitTrack(hits, z0, *geometry, track);
    
    if (track.IsValid())
    {
      DetSet<TotemRPLocalTrack> ds = output.find_or_insert(rpId);
      ds.push_back(track);
    }
  }

  /*
  RPTrackCandidateCollection::const_iterator in_it;
  for(in_it=input.begin(); in_it!=input.end(); ++in_it)
  {
    unsigned int rpId = in_it->first;

    double z0 = geometry->GetRPGlobalTranslation(rpId).z();

    TotemRPLocalTrack fitted_track;
    the_track_candidate_fitter_.FitTrack(in_it->second, z0, fitted_track, *geometry);

    if (fitted_track.IsValid())
      output[rpId] = fitted_track;
  }
  */

  // save results
  e.put(make_unique<DetSetVector<TotemRPLocalTrack>>(output));
}

//----------------------------------------------------------------------------------------------------
DEFINE_FWK_MODULE(RPTrackCandidateCollectionFitter);
