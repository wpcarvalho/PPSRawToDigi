import FWCore.ParameterSet.Config as cms

RP2StationInelastReconst = cms.EDProducer("RPPrimaryVertexInelasticReconstruction",

    Verbosity = cms.int32(0), RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit"),

    StripAlignmentResolutionDegradation = cms.double(1.7),

    ConstrainPrimaryVertex = cms.bool(True),
    ElasticScatteringReconstruction = cms.bool(False),

    ExternalPrimaryVertex = cms.bool(False),
    Station150MandatoryInReconstruction = cms.bool(True),
    Station220MandatoryInReconstruction = cms.bool(True),
    AnyStationInReconstruction = cms.bool(False),

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

    ComputeFullVarianceMatrix = cms.bool(True),

    RPMultipleScatteringSigma = cms.double(5.7e-07),    # rad

    ReconstructionPrecisionX = cms.double(0.18),    # mm
    ReconstructionPrecisionY = cms.double(0.07),    # mm
    ReconstructionPrecisionThetaX = cms.double(1e-05),  # rad
    ReconstructionPrecisionThetaY = cms.double(3e-07),  # rad
    ReconstructionPrecisionKsi = cms.double(0.007), # -1 .. 0

    InitMinX = cms.double(-0.6),
    InitMinY = cms.double(-0.6),
    InitMinThetaX = cms.double(-0.00045),   # reconstructed crossing angle is positive for both beams
    InitMinThetaY = cms.double(-0.00045),
    InitMinKsi = cms.double(-0.3),

    InitMaxX = cms.double(0.6),
    InitMaxY = cms.double(0.6),
    InitMaxThetaX = cms.double(0.00045),    # reconstructed crossing angle is positive for both beams
    InitMaxThetaY = cms.double(0.00045),
    InitMaxKsi = cms.double(0.0),

    RandomSearchProbability = cms.double(0.3),

    InitIterationsNumber = cms.int32(20),
    MaxChiSqNDFOfConvergedProton = cms.double(15.0),
    MaxChiSqOfConvergedInitialisation = cms.double(100.0),
    MaxAllowedReconstructedXi = cms.double(0.001),
    OutOfXiRangePenaltyFactor = cms.double(10000.0),
    InverseParamRandSeed = cms.int32(14142135),

    BeamProtTransportSetup = cms.PSet(),
    ExpectedRPResolution = cms.double(0.016)
)


