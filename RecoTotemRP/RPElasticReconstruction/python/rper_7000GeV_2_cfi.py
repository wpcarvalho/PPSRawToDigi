import FWCore.ParameterSet.Config as cms

ElasticReconstruction = cms.EDProducer("RPElasticReconstruction",
    verbosity = cms.untracked.uint32(0),
    
    # 0 ... do not add at all, 1 ... add only to left and right fits, 2 ... add to all fits
    addIPHit = cms.uint32(0),

    # road size (in radians)    
    roadsize_x = cms.double(0.001),
    roadsize_y = cms.double(0.0001),

    # angular tolerance (in radians)
    angleTolerance_y = cms.double(0.0001),
    angleTolerance_x = cms.double(0.0001),

    # vertex position tolerance (in mm)
    vertexTolerance_y = cms.double(0.4),
    vertexTolerance_x = cms.double(0.15),

    edgeEffectCorrection = cms.PSet(
        enabled = cms.bool(False),
        si_bd = cms.double(1.58e-05),
        th_L = cms.double(0.000182)
    )  ,
    RPFittedTrackCollLabel = cms.InputTag("RPSingleTrackCandCollFit")
)
