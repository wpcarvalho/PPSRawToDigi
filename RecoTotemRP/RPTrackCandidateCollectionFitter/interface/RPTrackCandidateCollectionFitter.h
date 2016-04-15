/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
* 	Hubert Niewiadomski
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#ifndef RecoTotemRP_RPTrackCandidateCollectionFitter_RPTrackCandidateCollectionFitter_h
#define RecoTotemRP_RPTrackCandidateCollectionFitter_RPTrackCandidateCollectionFitter_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"
#include "RecoTotemRP/RPTrackCandidateFitter/interface/RPTrackCandidateFitter.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"

#include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"


/**
 *\brief TODO: Fits every track candidate in RPTrackCandidateCollection.
 **/
class RPTrackCandidateCollectionFitter : public edm::EDProducer
{
  public:
    explicit RPTrackCandidateCollectionFitter(const edm::ParameterSet& conf);
    virtual ~RPTrackCandidateCollectionFitter();
    virtual void beginJob();
    virtual void produce(edm::Event& e, const edm::EventSetup& c);

  private:
    int verbosity_;

    /// Selection of the pattern-recognition module.
    edm::InputTag patternCollectionLabel;

    edm::EDGetTokenT<edm::DetSetVector<TotemRPUVPattern>> patternCollectionToken;
    
    /// A watcher to detect geometry changes.
    edm::ESWatcher<VeryForwardRealGeometryRecord> geometryWatcher;
    
    /// The instance of the fitter module
    RPTrackCandidateFitter fitter_;
};

#endif
