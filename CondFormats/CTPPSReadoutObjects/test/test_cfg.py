import FWCore.ParameterSet.Config as cms

process = cms.Process('test')

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# raw data source
process.source = cms.Source("EmptySource",
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# load a mapping
process.load("CondFormats.CTPPSReadoutObjects.DiamondDAQMappingESSourceXML_cfi")
process.DiamondDAQMappingESSourceXML.mappingFileNames.append("CondFormats/CTPPSReadoutObjects/xml/ctpps_timing_diamond_215_mapping.xml")
 
# print the mapping
process.printTotemDAQMappingDiamond = cms.EDAnalyzer("PrintTotemDAQMappingDiamond"
)

process.path = cms.Path(
  process.printTotemDAQMappingDiamond
)
