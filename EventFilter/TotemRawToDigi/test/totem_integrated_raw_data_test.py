import FWCore.ParameterSet.Config as cms

process = cms.Process("TotemIntegratedRawDataTest")

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# raw data source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/j/jkaspar/public/run268608_ls0001_streamA_StorageManager.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# raw-to-digi conversion
process.load('CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi')
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/ctpps_210_mapping.xml")

process.load("EventFilter.TotemRawToDigi.TotemTriggerRawToDigi_cfi")
process.TotemTriggerRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")

process.load('EventFilter.TotemRawToDigi.TotemRPRawToDigi_cfi')
process.TotemRPRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")
process.TotemRPRawToDigi.fedIds = cms.vuint32(578, 579, 580) # TODO: remove
process.TotemRPRawToDigi.RawToDigi.printErrorSummary = 0  # TODO: back to 1
process.TotemRPRawToDigi.RawToDigi.printUnknownFrameSummary = 0  # TODO: back to 1

# execution configuration
process.p = cms.Path(
    process.TotemTriggerRawToDigi *
    process.TotemRPRawToDigi
)
