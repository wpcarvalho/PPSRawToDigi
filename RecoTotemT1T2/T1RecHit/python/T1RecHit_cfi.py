import FWCore.ParameterSet.Config as cms

t1rechit = cms.EDProducer("T1RecHitProducer",
    Verbosity = cms.int32(0),
    ChiSquare = cms.double(2.0),
    T1ClusterCollectionLabel = cms.InputTag("t1cluster"),
    T1DigiWireCollectionLabel = cms.InputTag("T1Digis", "T1DigiWire"),
    T1RecHit2DCollectionLabel = cms.InputTag("t1rechit")
)


