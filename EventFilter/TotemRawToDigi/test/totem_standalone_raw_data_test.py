import FWCore.ParameterSet.Config as cms

process = cms.Process("TotemStandaloneRawDataTest")

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# raw data source
process.load('EventFilter.TotemRawToDigi.TotemStandaloneRawDataSource_cfi')
process.source.verbosity = 10
process.source.printProgressFrequency = 0
process.source.fileNames.append('/afs/cern.ch/user/j/jkaspar/public/run_9987_EVB11_1.003.srs')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# raw-to-digi conversion
process.load('CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi')
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/totem_rp_210far_220_mapping.xml")

process.load("EventFilter.TotemRawToDigi.TotemTriggerRawToDigi_cfi")
process.TotemTriggerRawToDigi.rawDataTag = cms.InputTag("source")
process.TotemTriggerRawToDigi.fedId = 0x29c

process.load('EventFilter.TotemRawToDigi.TotemRPRawToDigi_cfi')
process.TotemRPRawToDigi.rawDataTag = cms.InputTag("source")
process.TotemRPRawToDigi.fedIds = cms.vuint32(0x1a1, 0x1a2, 0x1a9, 0x1aa)
process.TotemRPRawToDigi.RawToDigi.printErrorSummary = 1
process.TotemRPRawToDigi.RawToDigi.printUnknownFrameSummary = 1

process.p = cms.Path(
    process.TotemTriggerRawToDigi *
    process.TotemRPRawToDigi
)
