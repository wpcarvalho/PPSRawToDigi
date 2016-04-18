import FWCore.ParameterSet.Config as cms

process = cms.Process('RECODQM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# minimum of logs
MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')  #for MC

# load DQM frame work
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")

# raw data source
source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/j/jkaspar/public/run268608_ls0001_streamA_StorageManager.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

# raw-to-digi conversion
process.load('CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi')
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/totem_rp_210far_220_mapping.xml")

process.load('EventFilter.TotemRawToDigi.TotemRPRawToDigi_cfi')
process.TotemRPRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")
process.TotemRPRawToDigi.fedIds = cms.vuint32(577, 578, 579, 580)
process.TotemRPRawToDigi.RawToDigi.printErrorSummary = 0
process.TotemRPRawToDigi.RawToDigi.printUnknownFrameSummary = 0

# local RP reconstruction chain with standard settings
process.load("RecoLocalCTPPS.TotemRP.LocalRecoChain_cfi")

# TOTEM RP DQM module
process.TotemRPDQMSource = cms.EDAnalyzer("TotemRPDQMSource",
    tagStripDigi = cms.InputTag("TotemRawToDigi"),
	tagDigiCluster = cms.InputTag("RPClustProd"),
	tagRecoHit = cms.InputTag("RPRecoHitProd"),
	tagPatternColl = cms.InputTag("NonParallelTrackFinder"),
	tagTrackColl = cms.InputTag("RPSingleTrackCandCollFit"),
	tagTrackCandColl = cms.InputTag("NonParallelTrackFinder"),
	tagMultiTrackColl = cms.InputTag(""),

    buildCorrelationPlots = cms.untracked.bool(False),
    correlationPlotsLimit = cms.untracked.uint32(50),
	correlationPlotsFilter = cms.untracked.string("default=0,1")
)

# DQM output
process.DQMOutput = cms.OutputModule("DQMRootOutputModule",
  fileName = cms.untracked.string("OUT_step1.root")
)

# execution schedule
process.reco_step = cms.Path(
  process.TotemRPRawToDigi *
  process.TotemRPLocalReconstruction
)

process.dqm_offline_step = cms.Path(
  process.TotemRPDQMSource
)

process.dqm_output_step = cms.EndPath(
    process.DQMOutput
)

process.schedule = cms.Schedule(
    process.reco_step,
    process.dqm_offline_step,
    process.dqm_output_step
)
