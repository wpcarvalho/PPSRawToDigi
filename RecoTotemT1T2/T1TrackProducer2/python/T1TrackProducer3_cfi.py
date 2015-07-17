import FWCore.ParameterSet.Config as cms

t1tracks3 = cms.EDProducer("T1TrackProducer3",
    Verbosity = cms.int32(0),
    ChiProbThr = cms.double(0.01),
    T1RoadCollectionLabel = cms.InputTag("t1roads")
)


