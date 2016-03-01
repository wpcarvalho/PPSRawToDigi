import FWCore.ParameterSet.Config as cms

process = cms.Process("rpReconstruction")
process.maxEvents = cms.untracked.PSet(
input = cms.untracked.int32(-1)
)

process.load('TotemRawData.Readers.RawDataSource_cfi')
process.source.fileNames = cms.untracked.vstring()

process.source.fileNames.append('root://eostotem//eos/totem/data/rawdata/2015/run_9998_EVB15_2.127.srs')


process.load("Configuration.TotemCommon.LoggerMin_cfi")

process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/2015_10_18_fill4511/RP_Dist_Beam_Cent.xml")

process.load("TotemAlignment.RPDataFormats.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.RealFiles = cms.vstring(
'/afs/cern.ch/exp/totem/scratch/Release/raw_data_reco/rev11326/CMSSW_7_0_4/src/TotemAlignment/RPData/LHC/2015_10_18_fill4511/version2/sr+el/45.xml',
'/afs/cern.ch/exp/totem/scratch/Release/raw_data_reco/rev11326/CMSSW_7_0_4/src/TotemAlignment/RPData/LHC/2015_10_18_fill4511/version2/sr+el/56.xml'
)

process.load('TotemCondFormats.DAQInformation.DAQMappingSourceXML_cfi')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_220_210far.xml')

process.load('TotemRawData.RawToDigi.Raw2DigiProducer_cfi')
process.Raw2DigiProducer.verbosity = 0

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_6500GeV_90_50urad_cfi")

process.load("RecoTotemRP.RPClusterSigmaService.ClusterSigmaServiceConf_cfi")
process.load("RecoTotemRP.RPClusterizer.RPClusterizationConf_cfi")
process.RPClustProd.DigiLabel = cms.InputTag("Raw2DigiProducer","rpDataOutput")

process.load("RecoTotemRP.RPRecoHitProducer.RPRecoHitProdConf_cfi")

process.load("RecoTotemRP.RPNonParallelTrackCandidateFinder.RPNonParallelTrackCandidateFinder_cfi")
process.NonParallelTrackFinder.verbosity = 0
process.NonParallelTrackFinder.maxHitsPerPlaneToSearch = 5
process.NonParallelTrackFinder.minPlanesPerProjectionToSearch = 2
process.NonParallelTrackFinder.minPlanesPerProjectionToFit = 3
process.NonParallelTrackFinder.threshold = 2.99

process.load("RecoTotemRP.RPTrackCandidateCollectionFitter.RPSingleTrackCandCollFitted_cfi")
process.RPSingleTrackCandCollFit.Verbosity = 0
process.RPSingleTrackCandCollFit.RPTrackCandCollProducer = 'NonParallelTrackFinder'

process.load("RecoTotemRP.RPInelasticReconstruction.Rec_6500GeV_beta_90_50urad_cfi")
process.RP220Reconst.BeamProtTransportSetup = process.BeamProtTransportSetup

process.output = cms.OutputModule(
"PoolOutputModule",
fileName = cms.untracked.string("file:./totcsi_9998.152127.root"),
outputCommands = cms.untracked.vstring('keep *')
)
process.p = cms.Path(
process.Raw2DigiProducer*
process.RPClustProd*
process.RPHecoHitProd*
process.NonParallelTrackFinder*
process.RPSingleTrackCandCollFit*
process.RP220Reconst
)

process.outpath = cms.EndPath(process.output)