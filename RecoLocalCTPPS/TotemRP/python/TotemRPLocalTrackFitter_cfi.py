import FWCore.ParameterSet.Config as cms

TotemRPLocalTrackFitter = cms.EDProducer("TotemRPLocalTrackFitter",
    Verbosity = cms.int32(0),

    tagUVPattern = cms.InputTag("TotemRPUVPatternFinder")
)
