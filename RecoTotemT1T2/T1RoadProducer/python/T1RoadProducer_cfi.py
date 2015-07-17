import FWCore.ParameterSet.Config as cms

t1roads = cms.EDProducer("T1RoadProducer",
    Verbosity = cms.int32(0),
    DeltaEta = cms.double(.1),
    DeltaPhi = cms.double(.2),
    MinHits = cms.int32(3),
    MaxHits = cms.int32(30),
    Alignment = cms.bool(False),
    T1RecHit2DCollectionLabel = cms.InputTag("t1rechit")
)


