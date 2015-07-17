import FWCore.ParameterSet.Config as cms

SmearingGenerator = cms.EDProducer("SmearingGenerator",
    verbosity = cms.untracked.uint32(1),
    modifyLabel = cms.string('generator'),
    originalLabel = cms.string('original'),

    # other parameters are now set at Configuration/TotemOpticsConfiguration/data/OpticsConfig_xxx.cfi, 
    # don't forget to include one of these files
)
