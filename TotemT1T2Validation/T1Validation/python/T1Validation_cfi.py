import FWCore.ParameterSet.Config as cms

t1valid = cms.EDAnalyzer("T1Validation",
    Verbosity = cms.int32(0),
    SIM = cms.double(1.0),
    DIGI = cms.double(1.0),
    RECO = cms.double(1.0),
    TrackLabel = cms.string('t1tracks'),
    OutputFile = cms.string('T1Valid.root'),
    T1ClusterCollectionLabel = cms.InputTag("t1cluster"),
    T1RoadCollectionLabel = cms.InputTag("t1roads"),
    T1DigiWireCollectionLabel = cms.InputTag("T1Digis", "T1DigiWire"),
    T1DigiVfatCollectionLabel = cms.InputTag("T1Digis", "T1DigiVfat")
)


