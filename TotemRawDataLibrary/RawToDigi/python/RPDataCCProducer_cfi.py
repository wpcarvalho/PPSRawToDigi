import FWCore.ParameterSet.Config as cms

RPDataCCProducer = cms.EDProducer("RPDataCCProducer",
    verbosity = cms.untracked.uint32(0),

    # default label of the raw product
    productLabelRaw = cms.string('RPCCRawBits')
)
