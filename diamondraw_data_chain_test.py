import FWCore.ParameterSet.Config as cms

process = cms.Process("DiamondRawDToDigiTest")
process.maxEvents = cms.untracked.PSet(
input = cms.untracked.int32(1000)
)
# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('DEBUG')

    )
)

# raw data source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/j/jkaspar/public/run273062_ls0001-2_stream.root')
)


# raw-to-digi conversion
process.load('CondFormats.CTPPSReadoutObjects.DiamondDAQMappingESSourceXML_cfi')
process.DiamondDAQMappingESSourceXML.mappingFileNames.append("CondFormats/CTPPSReadoutObjects/xml/ctpps_timing_diamond_215_mapping.xml")

process.load("EventFilter.TotemRawToDigi.totemTriggerRawToDigi_cfi")
process.totemTriggerRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")

process.load('EventFilter.CTPPSRawToDigi.diamondRPRawToDigi_cfi')
process.diamondRPRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")


process.dump = cms.EDAnalyzer("EventContentAnalyzer")

# ntuplizer
process.load("TotemAnalysis.TotemNtuplizer.TotemNtuplizer_cfi")
process.totemNtuplizer.outputFileName = "ntuple.root"

process.p = cms.Path(
    process.diamondRPRawToDigi 
)


# output configuration
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:./DiamondDigi.root"),
     outputCommands = cms.untracked.vstring(
    'keep DiamondFEDInfos_diamondRPRawToDigi_*_*',
    'keep DiamondDigiedmDetSetVector_diamondRPRawToDigi_*_*',
    'keep DiamondVFATStatusedmDetSetVector_diamondRPRawToDigi_*_*'
 )
)





process.outpath = cms.EndPath(process.output)
