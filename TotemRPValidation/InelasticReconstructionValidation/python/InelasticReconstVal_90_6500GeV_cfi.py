import FWCore.ParameterSet.Config as cms

RPInelProtRecVal = cms.EDAnalyzer("InelasticReconstructionValidation", RPReconstructedProtonCollectionLabel = cms.InputTag("RP220Reconst"),
    Verbosity = cms.int32(0),
    PrimaryProtonHepMCProductLabel = cms.string('original'),
    PrimaryProtonHepMCModuleName = cms.string('SmearingGenerator'),
    SmearedVertexHepMCProductLabel = cms.string(''),
    SmearedVertexHepMCModuleName = cms.string('generator'),
    HistogramFileName = cms.string('RPInelasticRecValHists.root'),
    RandomGenSeed = cms.int32(31415927),

    ParameterizationFileName220Right = cms.string('Geometry/TotemRPOptics/data/parametrization_6500GeV_90_reco.root'),
    ParameterizationFileName220Left = cms.string('Geometry/TotemRPOptics/data/parametrization_6500GeV_90_reco.root'),
    ParameterizationFileName150Right = cms.string('Geometry/TotemRPOptics/data/parametrization_6500GeV_90_reco.root'),
    ParameterizationFileName150Left = cms.string('Geometry/TotemRPOptics/data/parametrization_6500GeV_90_reco.root'),

    ParameterizationNamePrefix220Right = cms.string('ip5_to_station_220'),
    ParameterizationNamePrefix220Left = cms.string('ip5_to_station_220'),
    ParameterizationNamePrefix150Right = cms.string('ip5_to_station_150'),
    ParameterizationNamePrefix150Left = cms.string('ip5_to_station_150'),

    RightBeamPostfix = cms.string('lhcb1'),
    LeftBeamPostfix = cms.string('lhcb2'),
    SDValidation = cms.bool(False),
    Chi2UpperTailProbCut = cms.double(0.005)
)


