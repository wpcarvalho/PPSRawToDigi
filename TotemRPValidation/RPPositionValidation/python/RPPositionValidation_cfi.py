import FWCore.ParameterSet.Config as cms

RPPosVal = cms.EDAnalyzer("RPPositionValidation",
    Verbosity = cms.untracked.uint32(2),
)
