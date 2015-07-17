import FWCore.ParameterSet.Config as cms

t1roadsWOaplane = cms.EDProducer("T1RoadProducerWOaplane",
    Verbosity = cms.int32(0),
    DeltaEta = cms.double(.1),
    DeltaPhi = cms.double(.2),
    MinHits = cms.int32(3),
    MaxHits = cms.int32(30),
    Alignment = cms.bool(False),
#    AlignmentFile = cms.string("TotemAlignment/T1Alignment/test/T1Align.txt"),
Plane2Skip = cms.int32(0),
    T1RecHit2DCollectionLabel = cms.InputTag("t1rechit")
)


