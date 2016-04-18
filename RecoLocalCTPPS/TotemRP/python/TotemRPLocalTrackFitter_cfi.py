import FWCore.ParameterSet.Config as cms

TotemRPLocalTrackFitter = cms.EDProducer("TotemRPLocalTrackFitter",
    verbosity = cms.int32(0),

    tagUVPattern = cms.InputTag("TotemRPUVPatternFinder")
)
