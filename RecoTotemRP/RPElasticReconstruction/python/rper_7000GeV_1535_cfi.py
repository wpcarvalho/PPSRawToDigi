import FWCore.ParameterSet.Config as cms

ElasticReconstruction = cms.EDProducer("RPElasticReconstruction",
    verbosity = cms.untracked.uint32(0),

    # 0 ... do not add at all, 1 ... add only to left and right fits, 2 ... add to all fits
    addIPHit = cms.uint32(0),

    # road size (in radians)
    roadsize_x = cms.double(5e-06),
    roadsize_y = cms.double(8e-06),

    # angular tolerance (in radians)
    angleTolerance_y = cms.double(1.5e-06),
    angleTolerance_x = cms.double(1.5e-06),

    # vertex position tolerance (in mm)
    vertexTolerance_y = cms.double(5.0),
    vertexTolerance_x = cms.double(2.0),

    edgeEffectCorrection = cms.PSet(
        enabled = cms.bool(False),
        si_bd = cms.double(0),
        th_L = cms.double(0)
    )  ,
    RPFittedTrackCollLabel = cms.InputTag("RPSingleTrackCandCollFit")
)
