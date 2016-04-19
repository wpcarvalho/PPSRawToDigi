import FWCore.ParameterSet.Config as cms

process = cms.Process('RECODQM')

# minimum of logs
MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')  #for MC

# load DQM frame work
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")

# raw data source
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/j/jkaspar/public/run268608_ls0001_streamA_StorageManager.root')
#)

process.load('TotemRawData.Readers.TotemStandaloneRawDataSource_cfi')
process.source.verbosity = 10
process.source.printProgressFrequency = 0
process.source.fileNames.append('/afs/cern.ch/user/j/jkaspar/public/run_9987_EVB11_1.003.srs')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

# raw-to-digi conversion
process.load('CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi')
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/totem_rp_210far_220_mapping.xml")

# process.load('EventFilter.TotemRawToDigi.TotemRPRawToDigi_cfi')
# process.TotemRPRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")
# process.TotemRPRawToDigi.fedIds = cms.vuint32(577, 578, 579, 580)
# process.TotemRPRawToDigi.RawToDigi.printErrorSummary = 0
# process.TotemRPRawToDigi.RawToDigi.printUnknownFrameSummary = 0

process.load("EventFilter.TotemRawToDigi.TotemTriggerRawToDigi_cfi")
process.TotemTriggerRawToDigi.rawDataTag = cms.InputTag("source")
process.TotemTriggerRawToDigi.fedId = 0x29c

process.load('EventFilter.TotemRawToDigi.TotemRPRawToDigi_cfi')
process.TotemRPRawToDigi.rawDataTag = cms.InputTag("source")
process.TotemRPRawToDigi.fedIds = cms.vuint32(0x1a1, 0x1a2, 0x1a9, 0x1aa, 0x1b5, 0x1bd)
process.TotemRPRawToDigi.RawToDigi.printErrorSummary = 0
process.TotemRPRawToDigi.RawToDigi.printUnknownFrameSummary = 0

# RP geometry
process.load("Geometry.VeryForwardGeometry.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/VeryForwardData/data/RP_Garage/RP_Dist_Beam_Cent.xml")

# local RP reconstruction chain with standard settings
process.load("RecoLocalCTPPS.TotemRP.LocalRecoChain_cfi")

# TOTEM DQM modules
process.load("DQM.Totem.TotemDAQTriggerDQMSource_cfi")
process.load("DQM.Totem.TotemRPDQMSource_cfi")

# DQM output
process.DQMOutput = cms.OutputModule("DQMRootOutputModule",
  fileName = cms.untracked.string("OUT_step1.root")
)

# execution schedule
process.reco_step = cms.Path(
  process.TotemTriggerRawToDigi *
  process.TotemRPRawToDigi *
  process.TotemRPLocalReconstruction
)

process.dqm_produce_step = cms.Path(
  process.TotemDAQTriggerDQMSource *
  process.TotemRPDQMSource
)

process.dqm_output_step = cms.EndPath(
    process.DQMOutput
)

process.schedule = cms.Schedule(
    process.reco_step,
    process.dqm_produce_step,
    process.dqm_output_step
)
