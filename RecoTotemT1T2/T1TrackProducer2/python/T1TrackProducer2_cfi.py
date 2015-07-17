import FWCore.ParameterSet.Config as cms

t1tracks2 = cms.EDProducer("T1TrackProducer2",
    Verbosity = cms.int32(0),
    ChiRidThr = cms.double(40.0),
    T1RoadCollectionLabel = cms.InputTag("t1roads")
)


