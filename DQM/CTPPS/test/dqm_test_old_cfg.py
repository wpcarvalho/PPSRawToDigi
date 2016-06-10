import FWCore.ParameterSet.Config as cms

process = cms.Process('RECODQM')

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')  #for MC

# load DQM frame work
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load("DQMServices.Components.DQMStoreStats_cfi")

process.dqmEnv.subSystemFolder = 'CTPPS'

process.dqmSaver = cms.EDAnalyzer("DQMFileSaverOnline",
  producer = cms.untracked.string("DQM"),
  tag = cms.untracked.string("CTPPS"),
  path = cms.untracked.string(".")
)

# raw data source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/j/jkaspar/public/run273062_ls0001-2_stream.root')
)

process.maxEvents = cms.untracked.PSet(
    # TODO: revert to -1
    #input = cms.untracked.int32(-1)
    input = cms.untracked.int32(10)
)

# raw-to-digi conversion
process.load("CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi")
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/ctpps_210_mapping.xml")

process.load("EventFilter.TotemRawToDigi.totemTriggerRawToDigi_cfi")
process.totemTriggerRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")

process.load("EventFilter.TotemRawToDigi.totemRPRawToDigi_cfi")
process.totemRPRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")

process.totemRawToDigi = cms.Sequence(
  #process.totemTriggerRawToDigi *
  process.totemRPRawToDigi
)

# RP geometry
process.load("Geometry.VeryForwardGeometry.geometryRP_cfi")

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.TotemRPLocal.totemRPLocalReconstruction_cff")

# TOTEM DQM modules
process.load("DQM.CTPPS.totemDAQTriggerDQMSource_cfi")
process.load("DQM.CTPPS.totemRPDQMSource_cfi")
process.load("DQM.CTPPS.totemRPDQMHarvester_cfi")

# ntuplizer
process.load("TotemAnalysis.TotemNtuplizer.TotemNtuplizer_cfi")
process.totemNtuplizer.outputFileName = "dqm_ntuple.root"

# output configuration
from RecoCTPPS.Configuration.RecoCTPPS_EventContent_cff import RecoCTPPSRECO
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:./dqm_reco.root"),
    outputCommands = RecoCTPPSRECO.outputCommands
)

# execution schedule
process.path = cms.Path(
  process.totemRawToDigi *
  process.totemRPLocalReconstruction *

  process.totemDAQTriggerDQMSource *
  process.totemRPDQMSource *
  process.totemRPDQMHarvester *
  
  process.dqmEnv *
  process.dqmSaver

  #process.totemNtuplizer
)

process.outpath = cms.EndPath(process.output)
