import FWCore.ParameterSet.Config as cms

ElasticReconstruction = cms.EDProducer("RPElasticReconstruction",
    verbosity = cms.untracked.uint32(0),

    # 0 ... do not add at all, 1 ... add only to left and right fits, 2 ... add to all fits
    addIPHit = cms.uint32(0),

    # road size (in radians)    
    roadsize_x = cms.double(160E-6),        # 8 sigma
    roadsize_y = cms.double(160E-6),        # 8 sigma

    # angular tolerance (in radians)
    angleTolerance_y = cms.double(80E-6),   # 4 sigma
    angleTolerance_x = cms.double(80E-6),   # 4 sigma

    # vertex position tolerance (in mm)
    vertexTolerance_y = cms.double(0.14),   # 4 sigma
    vertexTolerance_x = cms.double(0.14),   # 4 sigma

    edgeEffectCorrection = cms.PSet(
        enabled = cms.bool(False),
        si_bd = cms.double(0),
        th_L = cms.double(0)
    )  ,
    RPFittedTrackCollLabel = cms.InputTag("RPSingleTrackCandCollFit")
)
