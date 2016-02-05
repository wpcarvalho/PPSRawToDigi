import FWCore.ParameterSet.Config as cms

RPDataDigiProducer = cms.EDProducer("RPDataDigiProducer",
    verbosity = cms.untracked.uint32(0),
    stopOnMissing = cms.untracked.bool(False)
)
