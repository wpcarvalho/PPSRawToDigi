import FWCore.ParameterSet.Config as cms

TotemRPClusterProducer = cms.EDProducer("TotemRPClusterProducer",
    Verbosity = cms.int32(0),
    tagDigi = cms.InputTag("TotemRPRawToDigi", "RP")
)
