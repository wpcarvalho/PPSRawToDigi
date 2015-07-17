import FWCore.ParameterSet.Config as cms

GeneratorValidation = cms.EDAnalyzer("GeneratorValidation",
    verbosity = cms.untracked.uint32(1),
    thetaLimit = cms.double(0.001),   # rad
    E_nom = cms.double(3500),         # GeV
    outputFile = cms.string('GeneratorValidation.root'),

    HepMCProductLabel = cms.InputTag("generator")
)
