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
process.source.fileNames.append('sample.srs')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# raw to digi conversion
process.load('CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi')
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/totem_rp_210far_220_mapping.xml")

process.load('EventFilter.TotemRawToDigi.TotemRawToDigi_cfi')
process.TotemRawToDigi.rawDataTag = cms.InputTag("source")
process.TotemRawToDigi.RawToDigi.printErrorSummary = 1
process.TotemRawToDigi.RawToDigi.printUnknownFrameSummary = 1

# clusterization
process.load("RecoTotemRP.RPClusterSigmaService.ClusterSigmaServiceConf_cfi")
process.load("RecoTotemRP.RPClusterizer.RPClusterizationConf_cfi")
process.RPClustProd.DigiLabel = cms.InputTag("TotemRawToDigi")

# reco hit production
process.load("RecoTotemRP.RPRecoHitProducer.RPRecoHitProdConf_cfi")

# geometry
process.load("Geometry.VeryForwardGeometry.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/VeryForwardData/data/RP_Garage/RP_Dist_Beam_Cent.xml")

process.load("Alignment.RPDataFormats.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.RealFiles = cms.vstring()

# non-parallel pattern recognition
process.load("RecoTotemRP.RPNonParallelTrackCandidateFinder.RPNonParallelTrackCandidateFinder_cfi")
process.NonParallelTrackFinder.verbosity = 0
process.NonParallelTrackFinder.maxHitsPerPlaneToSearch = 5
process.NonParallelTrackFinder.minPlanesPerProjectionToSearch = 2
process.NonParallelTrackFinder.minPlanesPerProjectionToFit = 3
process.NonParallelTrackFinder.threshold = 2.99

# local track fitting
process.load("RecoTotemRP.RPTrackCandidateCollectionFitter.RPSingleTrackCandCollFitted_cfi")
process.RPSingleTrackCandCollFit.Verbosity = 0
process.RPSingleTrackCandCollFit.RPTrackCandCollProducer = 'NonParallelTrackFinder' # takes up the non-parallel pattern recognition

process.p = cms.Path(
    process.TotemRawToDigi *
    process.RPClustProd *
    process.RPRecoHitProd
#    process.NonParallelTrackFinder *
#    process.RPSingleTrackCandCollFit
)

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:./out.root")
)

process.outpath = cms.EndPath(process.output)
