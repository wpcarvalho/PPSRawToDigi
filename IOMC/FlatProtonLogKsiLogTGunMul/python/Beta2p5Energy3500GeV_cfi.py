import FWCore.ParameterSet.Config as cms

source = cms.Source("FlatProtonLogKsiLogTGunMul",
    LogKsiDistShape = cms.untracked.bool(True),
    LogTDistShape = cms.untracked.bool(True),
    GenerateMirroredProtons = cms.untracked.bool(False),

    RightArm = cms.untracked.bool(True),
    LeftArm = cms.untracked.bool(True),
    NominalEnergy = cms.untracked.double(3500.0),
    MinT = cms.untracked.double(-0.0001),
    MaxT = cms.untracked.double(1.0),
    MinPhi = cms.untracked.double(-3.141592654),
    MaxPhi = cms.untracked.double(3.141592654),
    MinKsi = cms.untracked.double(-0.01),
    MaxKsi = cms.untracked.double(-0.15),
    Verbosity = cms.untracked.int32(0),
    MaximumIterationNumber = cms.untracked.int32(20)
)


