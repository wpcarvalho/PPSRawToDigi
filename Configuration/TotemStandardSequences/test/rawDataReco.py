import FWCore.ParameterSet.Config as cms

process = cms.Process("HubertBigTest")

process.maxEvents = cms.untracked.PSet(
input = cms.untracked.int32(-1)
)

process.load('TotemRawData.Readers.RawDataSource_cfi')
process.source.verbosity = 1
process.source.printProgressFrequency = 0

process.source.fileNames = cms.untracked.vstring()
# process.source.fileNames.append('/castor/cern.ch/totem/LHCRawData/2013/Physics/run_8781.000.vmea')
process.source.fileNames.append('/castor/cern.ch/user/k/kmielnik/run_8318.000.vmea')


process.source.skipEvents = cms.untracked.uint32(0)

process.load("Configuration.TotemCommon.LoggerMin_cfi")

process.load("Configuration/TotemOpticsConfiguration/OpticsConfig_4000GeV_90_cfi")

process.load("TotemAlignment.RPDataFormats.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.RealFiles = cms.vstring('/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAlignment/RPData/LHC/2012_07_12-13/sr+hsx/45_220.xml','/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAlignment/RPData/LHC/2012_07_12-13/sr+hsx/56_220.xml','/afs/cern.ch/exp/totem/scratch/jkaspar/software/offline/424/src/TotemAlignment/RPData/LHC/2012_07_12-13/el/corrections.xml')

process.load('TotemCondFormats/DAQInformation/DAQMappingSourceXML_cfi')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_220.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_147.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/t1_all_run1.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/t2_4quarters.xml')
process.DAQMappingSourceXML.maskFileNames.append('TotemCondFormats/DAQInformation/data/T1DeadNoisyChannelsListRun1.xml')

process.load('TotemRawData.RawToDigi.Raw2DigiProducer_cfi')
process.Raw2DigiProducer.verbosity = 10

process.load("RecoTotemRP.RPClusterSigmaService.ClusterSigmaServiceConf_cfi")
process.load("RecoTotemRP.RPClusterizer.RPClusterizationConf_cfi")

process.load("RecoTotemRP.RPRecoHitProducer.RPRecoHitProdConf_cfi")

process.load("Configuration/TotemCommon/geometryGlobal_real_cfi")
toberemoved = []
for xmlfile in process.XMLIdealGeometryESSource.geomXMLFiles:
    if xmlfile.endswith("RP_Dist_Beam_Cent.xml"):
       toberemoved.append(xmlfile)
for xmlfile in toberemoved:
    process.XMLIdealGeometryESSource.geomXMLFiles.remove(xmlfile)
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/2012_07_12-13/RP_Dist_Beam_Cent.xml")
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

process.load("RecoTotemRP.RPSingleCandidateTrackFinder.RPSingleTrackCandFindConf_cfi")
process.RPSinglTrackCandFind.Verbosity = 0
process.RPSinglTrackCandFind.RoadSize = 0.3
process.RPSinglTrackCandFind.MinHitsPerCoord = 3
process.RPSinglTrackCandFind.MaxHitsPerDetector = 5
process.RPSinglTrackCandFind.ReduceWeightsWithMultiplicity = False

process.load("RecoTotemRP.RPNonParallelTrackCandidateFinder.RPNonParallelTrackCandidateFinder_cfi")
process.NonParallelTrackFinder.verbosity = 0
process.NonParallelTrackFinder.maxHitsPerPlaneToSearch = 5
process.NonParallelTrackFinder.minPlanesPerProjectionToSearch = 3
process.NonParallelTrackFinder.minPlanesPerProjectionToFit = 3
process.NonParallelTrackFinder.threshold = 2.99

process.load("RecoTotemRP.RPTrackCandidateCollectionFitter.RPSingleTrackCandCollFitted_cfi")
process.RPSingleTrackCandCollFit.Verbosity = 0

process.RPSingleTrackCandCollFit.RPTrackCandCollProducer = 'NonParallelTrackFinder'

process.load("RecoTotemRP/RPInelasticReconstruction/Rec_4000GeV_beta_90p0_220_2Arm_cfi")
process.RP2202ArmReconst.BeamProtTransportSetup = process.BeamProtTransportSetup
process.RP2202ArmReconst.ExpectedRPResolution = 0.020
process.RP2202ArmReconst.Verbosity = 0

process.load("RecoTotemRP/RPInelasticReconstruction/Rec_4000GeV_beta_90p0_220_cfi")
process.RP220Reconst.BeamProtTransportSetup = process.BeamProtTransportSetup
process.RP220Reconst.ExpectedRPResolution = 0.020
process.RP220Reconst.Verbosity = 0
process.RP220Reconst.ElasticScatteringReconstruction = False

process.TriggerBits = cms.EDProducer("RPTriggerBitsProducer",
verbose = cms.bool(False)
)


process.load('L1TriggerTotem.CoincidenceChip.RPCoincidenceProducer_cfi')

process.load("RecoTotemRP.RPMulTrackCandidateCollectionFitter.RPMulTrackNonParallelRecoFitter_cfi")

process.RawDataDumper = cms.EDAnalyzer("RawDataDumper")
process.RawDataDumper.RawEventLabel = cms.InputTag("source")


process.Raw2DigiProducer.rpDataProductLabel = cms.untracked.string("")
process.Raw2DigiProducer.rpCCProductLabel = cms.untracked.string("")
process.TriggerBits.StripDigiLabel = cms.InputTag("Raw2DigiProducer")
process.RPCC.DetTriggerLabel = cms.InputTag("TriggerBits")
process.RPClustProd.DigiLabel  = cms.InputTag("Raw2DigiProducer")

process.p = cms.Path(
#    process.RawDataDumper
     process.Raw2DigiProducer
    *process.TriggerBits
    *process.RPCC
    *process.RPClustProd
    *process.RPHecoHitProd
    *process.RPSinglTrackCandFind
    *process.NonParallelTrackFinder
    *process.RPSingleTrackCandCollFit
    *process.RPMulTrackNonParallelCandCollFit
    *process.RP220Reconst
    *process.RP2202ArmReconst
)

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:file:./rawDataReco.root"),
    outputCommands = cms.untracked.vstring('keep *')
)

process.outpath = cms.EndPath(process.output)
