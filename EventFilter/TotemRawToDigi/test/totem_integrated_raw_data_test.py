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
# TODO: put here the CMS raw-data source module

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)  # TODO: number of events to process
)

# raw to digi conversion
process.load('CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi')
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/totem_rp_210far_220_mapping.xml")

process.load('EventFilter.TotemRawToDigi.TotemRawToDigi_cfi')
process.TotemRawToDigi.rawDataTag = cms.InputTag("source")  # TODO: this may need tuning
process.TotemRawToDigi.fedIds = cms.vuint32(577, 578, 579, 580, 581)  # TODO: this may need tuning
process.TotemRawToDigi.RawToDigi.printErrorSummary = 1
process.TotemRawToDigi.RawToDigi.printUnknownFrameSummary = 1

# execution configuration
process.p = cms.Path(
    process.TotemRawToDigi
)
