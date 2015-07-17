import FWCore.ParameterSet.Config as cms

SmearingValidation = cms.EDAnalyzer("SmearingValidation",
    verbosity = cms.untracked.uint32(1),
    generatorLabel = cms.string('generator'),
    originalLabel = cms.string('original'),
    outputFile = cms.string('SmearingValidation.root'),
    thetaLimit = cms.double(0.001)
)
