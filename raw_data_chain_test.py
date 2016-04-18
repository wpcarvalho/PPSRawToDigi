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
process.load('TotemRawData.Readers.TotemStandaloneRawDataSource_cfi')
process.source.verbosity = 10
process.source.printProgressFrequency = 0
process.source.fileNames.append('/afs/cern.ch/user/j/jkaspar/public/run_9987_EVB11_1.003.srs')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1) # TODO: to -1
)

# raw-to-digi conversion
process.load('CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi')
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/totem_rp_210far_220_mapping.xml")

process.load("EventFilter.TotemRawToDigi.TotemTriggerRawToDigi_cfi")
process.TotemTriggerRawToDigi.rawDataTag = cms.InputTag("source")
process.TotemTriggerRawToDigi.fedId = 0x29c

process.load('EventFilter.TotemRawToDigi.TotemRPRawToDigi_cfi')
process.TotemRPRawToDigi.rawDataTag = cms.InputTag("source")
process.TotemRPRawToDigi.fedIds = cms.vuint32(0x1a1, 0x1a2, 0x1a9, 0x1aa, 0x1b5, 0x1bd)
process.TotemRPRawToDigi.RawToDigi.printErrorSummary = 1
process.TotemRPRawToDigi.RawToDigi.printUnknownFrameSummary = 0

# geometry
process.load("Geometry.VeryForwardGeometry.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/VeryForwardData/data/RP_Garage/RP_Dist_Beam_Cent.xml")

process.load("Alignment.RPDataFormats.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.RealFiles = cms.vstring()

# clusterization
process.load("RecoLocalCTPPS.TotemRP.TotemRPClusterProducer_cfi")

# reco hit production
process.load("RecoLocalCTPPS.TotemRP.TotemRPRecHitProducer_cfi")

# non-parallel pattern recognition
process.load("RecoLocalCTPPS.TotemRP.TotemRPUVPatternFinder_cfi")

# local track fitting
process.load("RecoLocalCTPPS.TotemRP.TotemRPLocalTrackFitter_cfi")

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
    process.TotemTriggerRawToDigi *
    process.TotemRPRawToDigi *
    process.TotemRPClusterProducer *
    process.TotemRPRecHitProducer *
    process.TotemRPUVPatternFinder *
    process.TotemRPLocalTrackFitter
)

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:./out.root")
)

process.outpath = cms.EndPath(process.output)
