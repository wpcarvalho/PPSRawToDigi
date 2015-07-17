import FWCore.ParameterSet.Config as cms

T2TrackColl3B = cms.EDProducer("T2TrackProducer3",
    ChiRidThr = cms.double(2.0),    #  Chi2 Threshold
    RoadModuleLabel = cms.string('T2RoadColl'),
    RoadInstanceLabel= cms.string('T2RoadColl'),
    theT2Hit_to_Track_Map0_Module= cms.string('HitToRoadAssociator'),
    theT2Hit_to_Track_Map0_Instance= cms.string('HitToRoadAssociator'),
    TrackInstanceLabel= cms.string('T2TrackColl3'),
    HitLabel = cms.string('T2Hits'),
    MinHitFinalTrk = cms.uint32(3),
    MinHitFinalTrkCl1Hit = cms.uint32(3),    
    verbosity = cms.bool(False),                       
    forceRZfit=cms.bool(False),
    forceXYfit=cms.bool(True),
    MaxHitsInRoad=cms.uint32(40),    
    UseRefittingProcedure=cms.bool(True),
    RecoHitRErrorPerStripcount =cms.vdouble(0.09,0.134,0.109,0.109,0.118),
    DropWorstRecHitChisquareRThreshold = cms.double(0.01),
    StripFitting=cms.bool(False),                           
    VtxPositionEX = cms.double(10),
    VtxPositionEY = cms.double(10),
    VtxPositionEZ = cms.double(2000),
    VtxPositionX  = cms.double(0),
    VtxPositionY  = cms.double(0),
    VtxPositionZ  = cms.double(0),
    FitVertexPosition=cms.bool(False), 
    RemoveOutliers=cms.bool(False),
    GhostSuppression=cms.bool(False), 
    PickUpDisplacedHit=cms.bool(False),
    PickUpRadius=cms.double(0.)

)


