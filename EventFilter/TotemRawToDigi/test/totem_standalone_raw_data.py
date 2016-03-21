import FWCore.ParameterSet.Config as cms

process = cms.Process("TotemStandaloneRawDataTest")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(3)
)

process.load('EventFilter.TotemRawToDigi.TotemStandaloneRawDataSource_cfi')
process.source.verbosity = 10
process.source.printProgressFrequency = 0
# eos cp /eos/totem/data/rawdata/2015/run_9998_EVB11_1.000.srs .
process.source.fileNames.append('run_9998_EVB11_1.000.srs')

# raw to digi conversion
process.load('CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi')
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/totem_rp_210_mapping.xml")
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/totem_rp_220_mapping.xml")

process.load('EventFilter.TotemRawToDigi.TotemRawToDigi_cfi')
process.TotemRawToDigi.rawDataTag = cms.InputTag("source")
process.TotemRawToDigi.RawToDigi.printErrorSummary = 0
process.TotemRawToDigi.RawToDigi.printUnknownFrameSummary = 0

process.PrintTotemDAQMapping = cms.EDAnalyzer("PrintTotemDAQMapping")

process.p = cms.Path(
    #process.PrintTotemDAQMapping *
    process.TotemRawToDigi
)
