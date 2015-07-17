import FWCore.ParameterSet.Config as cms

t1tracks4 = cms.EDProducer("T1TrackProducer4",
    Verbosity = cms.int32(0),
    ChiProbThr = cms.double(0.01),
    eVx = cms.double(10),
    eVy = cms.double(10),
    eVz = cms.double(500),
    T1RoadCollectionLabel = cms.InputTag("t1roads")
)


