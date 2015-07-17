import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Reconstruction_cff import *

# Ferenc Tracking --------------------------------------------------------------------
from RecoPixelVertexing.PixelLowPtUtilities.MinBiasTracking_cff import *
#     This is only for MC, if Rechits are not present
from RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi import *
#myTTRHBuilderWithoutAngle4PixelTriplets.ComputeCoarseLocalPositionFromDisk = True

# Merging mbTracks & generalTracks
import RecoTracker.FinalTrackSelectors.simpleTrackListMerger_cfi
generalPlusMinBiasTracks = RecoTracker.FinalTrackSelectors.simpleTrackListMerger_cfi.simpleTrackListMerger.clone(
    TrackProducer1 = 'generalTracks',
    TrackProducer2 = 'allTracks',
    promoteTrackQuality = True,
    newQuality = 'goodIterative'
    #promoteTrackQuality = False
    )

# Filtering to remove loopers
#load('myProducers/V0TrackFilter/v0trackfilter_cfg')
#generalPlusMinBiasFilteredTracks.pxCut = cms.double(0.04)
#generalPlusMinBiasFilteredTracks.pyCut = cms.double(0.04)
#generalPlusMinBiasFilteredTracks.pzCut = cms.double(0.025)

# Ferenc vertex on Tracks ------------------------------------------------------------

import UserCode.FerencSiklerVertexing.NewVertexProducer_cfi

generalVertices = UserCode.FerencSiklerVertexing.NewVertexProducer_cfi.newVertices.clone()
generalVertices.TrackCollection = 'generalTracks'
generalVertices.PtMin = cms.double(0.0)

allVertices = UserCode.FerencSiklerVertexing.NewVertexProducer_cfi.newVertices.clone()
allVertices.TrackCollection = 'allTracks'
allVertices.PtMin = cms.double(0.0)

mergedVertices = UserCode.FerencSiklerVertexing.NewVertexProducer_cfi.newVertices.clone()
mergedVertices.TrackCollection = 'generalPlusMinBiasTracks'
mergedVertices.PtMin = cms.double(0.0)

# Standard vertexing on merged tracks ------------------------------------------------------------

from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *
offlinePrimaryVerticesWithMBTracks = offlinePrimaryVertices.clone(
    TrackLabel = cms.InputTag('generalPlusMinBiasTracks')
)
#offlinePrimaryVerticesWithMBTracks.PVSelParameters.maxDistanceToBeam = 2  # ! changed in 5.2.4
offlinePrimaryVerticesWithMBTracks.vertexCollections.maxDistanceToBeam = 2  
offlinePrimaryVerticesWithMBTracks.TkFilterParameters.maxNormalizedChi2 = 20
offlinePrimaryVerticesWithMBTracks.TkFilterParameters.maxD0Significance = 100
offlinePrimaryVerticesWithMBTracks.TkFilterParameters.minPixelLayersWithHits = 2
offlinePrimaryVerticesWithMBTracks.TkFilterParameters.minSiliconLayersWithHits = 5
#offlinePrimaryVerticesWithMBTracks.TkClusParameters.zSeparation = 1
#offlinePrimaryVerticesWithMBTracks.TkClusParameters.TkGapClusParameters.zSeparation = 1
#offlinePrimaryVerticesWithMBTracks.TkGapClusParameters.zSeparation = 1

redoSiHits = cms.Sequence(
                         siPixelRecHits *
                         siStripMatchedRecHits
                       )


fulltracking = cms.Sequence(minBiasTracking *
                          generalPlusMinBiasTracks *
#                        generalPlusMinBiasFilteredTracks *

                         #Vertices
                         allVertices *
                         generalVertices *
                         mergedVertices *
                         offlinePrimaryVerticesWithMBTracks
)
