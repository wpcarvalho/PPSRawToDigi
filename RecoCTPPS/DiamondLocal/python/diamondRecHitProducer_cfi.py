import FWCore.ParameterSet.Config as cms

diamondRecHitProducer = cms.EDProducer("DiamondRecHitProducer",
    verbosity = cms.int32(1),
    tagDigi = cms.InputTag("diamondRPRawToDigi", "RP"),
    subSystem = cms.string("RP")
)
