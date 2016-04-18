import FWCore.ParameterSet.Config as cms

TotemRPRecHitProducer = cms.EDProducer("TotemRPRecHitProducer",
    Verbosity = cms.int32(0),
    tagCluster = cms.InputTag("TotemRPClusterProducer")
)
