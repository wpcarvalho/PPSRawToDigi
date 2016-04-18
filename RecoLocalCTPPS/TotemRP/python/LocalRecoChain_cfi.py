import FWCore.ParameterSet.Config as cms

# clusterization
from RecoLocalCTPPS.TotemRP.TotemRPClusterProducer_cfi import *

# reco hit production
from RecoLocalCTPPS.TotemRP.TotemRPRecHitProducer_cfi import *

# non-parallel pattern recognition
from RecoLocalCTPPS.TotemRP.TotemRPUVPatternFinder_cfi import *

# local track fitting
from RecoLocalCTPPS.TotemRP.TotemRPLocalTrackFitter_cfi import *

TotemRPLocalReconstruction = cms.Sequence(
    TotemRPClusterProducer *
    TotemRPRecHitProducer *
    TotemRPUVPatternFinder *
    TotemRPLocalTrackFitter
)
