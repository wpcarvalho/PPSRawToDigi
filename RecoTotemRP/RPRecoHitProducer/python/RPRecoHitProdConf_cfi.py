import FWCore.ParameterSet.Config as cms

RPRecoHitProd = cms.EDProducer("RPRecoHitProducer",
    Verbosity = cms.int32(1),
    ClusterLabel = cms.InputTag("RPClustProd")
)


