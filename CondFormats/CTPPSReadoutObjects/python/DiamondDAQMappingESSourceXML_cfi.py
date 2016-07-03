import FWCore.ParameterSet.Config as cms

DiamondDAQMappingESSourceXML = cms.ESSource("DiamondDAQMappingESSourceXML",
  verbosity = cms.untracked.uint32(0),

  mappingFileNames = cms.untracked.vstring(),
  maskFileNames = cms.untracked.vstring()
)
