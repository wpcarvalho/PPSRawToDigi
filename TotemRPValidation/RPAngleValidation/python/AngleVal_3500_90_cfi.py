import FWCore.ParameterSet.Config as cms

RPAngleVal = cms.EDAnalyzer("RPAngleValidation",
    Verbosity = cms.int32(0),

    SmearedHepMCProductLabel = cms.string(''),
    SmearedHepMCModuleName = cms.string('generator'),

    OriginalHepMCProductLabel = cms.string('original'),
    OriginalHepMCModuleName = cms.string('SmearingGenerator'),

    HistogramFileName = cms.string('RPAngleValidationHistsBeta90Energy3500GeV.root'),
    RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit")
)


