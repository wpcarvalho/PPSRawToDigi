import FWCore.ParameterSet.Config as cms

RPRecTracksVal = cms.EDAnalyzer("RPReconstructedTracksValidation",

    Verbosity = cms.int32(0),

    StripAlignmentResolutionDegradation = cms.double(1.7),

    SmearedHepMCProductLabel = cms.string(''),
    SmearedHepMCModuleName = cms.string('generator'),

    OriginalHepMCProductLabel = cms.string('original'),
    OriginalHepMCModuleName = cms.string('SmearingGenerator'),
    SDValidation = cms.bool(False),
    ProductLabelSimu = cms.string('RPCCSimuBits'),
    ModulLabelSimu = cms.string('RPCC'),

    HistogramFileName = cms.string('RPTracksValidationHistsBeta90.root'),

    ParameterizationFileName220Right = cms.string('Geometry/TotemRPOptics/data/parametrization_6500GeV_90_reco.root'),
    ParameterizationFileName220Left = cms.string('Geometry/TotemRPOptics/data/parametrization_6500GeV_90_reco.root'),
    ParameterizationFileName150Right = cms.string('Geometry/TotemRPOptics/data/parametrization_6500GeV_90_reco.root'),
    ParameterizationFileName150Left = cms.string('Geometry/TotemRPOptics/data/parametrization_6500GeV_90_reco.root'),

    ParameterizationNamePrefix220Right = cms.string('ip5_to_station_220'),
    ParameterizationNamePrefix220Left = cms.string('ip5_to_station_220'),
    ParameterizationNamePrefix150Right = cms.string('ip5_to_station_150'),
    ParameterizationNamePrefix150Left = cms.string('ip5_to_station_150'),

    RightBeamPostfix = cms.string('lhcb1'),
    LeftBeamPostfix = cms.string('lhcb1'),

    RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit"),
    RPDigClusterSetLabel = cms.InputTag("RPClustProd"),
    RPDetTriggerSetLabel = cms.InputTag("RPSiDetDigitizer")
)


