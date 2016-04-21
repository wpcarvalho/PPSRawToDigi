import FWCore.ParameterSet.Config as cms

# clusterization
from RecoLocalCTPPS.TotemRP.totemRPClusterProducer_cfi import *

# reco hit production
from RecoLocalCTPPS.TotemRP.totemRPRecHitProducer_cfi import *

# non-parallel pattern recognition
from RecoLocalCTPPS.TotemRP.totemRPUVPatternFinder_cfi import *

# local track fitting
from RecoLocalCTPPS.TotemRP.totemRPLocalTrackFitter_cfi import *

totemRPLocalReconstruction = cms.Sequence(
    totemRPClusterProducer *
    totemRPRecHitProducer *
    totemRPUVPatternFinder *
    totemRPLocalTrackFitter
)
