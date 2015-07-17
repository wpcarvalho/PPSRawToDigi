import FWCore.ParameterSet.Config as cms

from RecoTotemRP.RPClusterSigmaService.ClusterSigmaServiceConf_cfi import *
RPHecoHitProd = cms.EDProducer("RPRecoHitProducer",
    Verbosity = cms.int32(0),
    ClusterLabel = cms.InputTag("RPClustProd")
)


