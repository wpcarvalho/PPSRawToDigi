import FWCore.ParameterSet.Config as cms

t1cluster = cms.EDProducer("T1MakeCluster",
    Verbosity = cms.int32(0),
    Electronics = cms.string('VFAT'),
    WEIGHT1 = cms.double(50.0),
    WEIGHT2 = cms.double(50.0),
    ActivateDeadChannels = cms.bool(False),
#    T1DigiSantiardCollectionLabel = cms.InputTag("T1Digis", "T1DigiSantiard"),
    T1DigiVfatCollectionLabel = cms.InputTag("T1Digis", "T1DigiVfat")
)


