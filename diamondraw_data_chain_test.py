import FWCore.ParameterSet.Config as cms

process = cms.Process("DiamondRawDToDigiTest")
process.maxEvents = cms.untracked.PSet(
input = cms.untracked.int32(500)
)
# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet( threshold = cms.untracked.string('DEBUG') )
)

# raw data source
process.source = cms.Source("NewEventStreamFileReader",
    #fileNames = cms.untracked.vstring('/store/t0streamer/Minidaq/A/000/276/395/run276395_ls0001_streamA_StorageManager.dat')
    #fileNames = cms.untracked.vstring('/store/t0streamer/Minidaq/A/000/278/884/run278884_ls0001_streamA_StorageManager.dat')
    fileNames = cms.untracked.vstring('/store/t0streamer/Minidaq/A/000/280/102/run280102_ls0045_streamA_StorageManager.dat') 
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
#process.load("TotemAnalysis.TotemNtuplizer.TotemNtuplizer_cfi")
#process.totemNtuplizer.outputFileName = "ntuple.root"

process.p = cms.Path(
    process.totemTriggerRawToDigi
    * process.diamondRPRawToDigi
#    * process.dump
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
