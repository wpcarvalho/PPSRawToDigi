import FWCore.ParameterSet.Config as cms

process = cms.Process('RECODQM')

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# load DQM framework
process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = "CTPPS"
process.dqmEnv.eventInfoFolder = "EventInfo"
process.dqmSaver.path = ""
process.dqmSaver.tag = "CTPPS"

## raw data source
#process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/j/jkaspar/public/run273062_ls0001-2_stream.root')
#)

# raw data source
process.source = cms.Source("NewEventStreamFileReader",
    fileNames = cms.untracked.vstring('/store/t0streamer/Minidaq/A/000/276/390/run276390_ls0001_streamA_StorageManager.dat')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# raw-to-digi conversion
process.load('CondFormats.CTPPSReadoutObjects.DiamondDAQMappingESSourceXML_cfi')
process.DiamondDAQMappingESSourceXML.mappingFileNames.append("CondFormats/CTPPSReadoutObjects/xml/ctpps_timing_diamond_215_mapping.xml")

process.load("EventFilter.TotemRawToDigi.totemRawToDigi_cff")

process.load('EventFilter.CTPPSRawToDigi.diamondRPRawToDigi_cfi')
process.diamondRPRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")
process.diamondRPRawToDigi.RawToDigi.verbosity = cms.untracked.uint32(0)
# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")

# CTPPS DQM modules
process.load("DQM.CTPPS.totemDQM_cff")

process.path = cms.Path(
  process.totemTriggerRawToDigi *
  (process.totemRPRawToDigi +
  process.diamondRPRawToDigi )*

  process.recoCTPPS *

  process.totemDQM
)

process.end_path = cms.EndPath(
  process.dqmEnv +
  process.dqmSaver
)

process.schedule = cms.Schedule(
  process.path,
  process.end_path
)
