import FWCore.ParameterSet.Config as cms

process = cms.Process("DQM")

#----------------------------
#### Event Source
#----------------------------
# for live online DQM in P5
# TODO: uncomment
#process.load("DQM.Integration.config.inputsource_cfi")

# for testing in lxplus
# TODO: comment
process.load("DQM.Integration.config.fileinputsource_cfi")
process.source.fileNames = cms.untracked.vstring(
  "file:/afs/cern.ch/user/j/jkaspar/public/run268608_ls0001_streamA_StorageManager.root"
)

#----------------------------
#### DQM Environment
#----------------------------
process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = 'Totem'
process.dqmSaver.tag = 'Totem'

process.dqmSaver.path = "." # TODO: remove this line, only for testing

#-----------------------------
process.load("DQMServices.Components.DQMProvInfo_cfi")

# message logger
process.MessageLogger = cms.Service("MessageLogger",
  destinations = cms.untracked.vstring('cout'),
  cout = cms.untracked.PSet(threshold = cms.untracked.string('WARNING'))
)

# Global tag - Condition for P5 cluster
process.load("DQM.Integration.config.FrontierCondition_GT_cfi")

# raw-to-digi conversion
process.load('CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi')
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/ctpps_210_mapping.xml")

# TODO: uncomment if/once data contain LoneG frame
# process.load('EventFilter.TotemRawToDigi.TotemRPRawToDigi_cfi')
# process.totemRPRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")
# process.totemRPRawToDigi.fedIds = cms.vuint32(577, 578, 579, 580)
# process.totemRPRawToDigi.RawToDigi.printErrorSummary = 0
# process.totemRPRawToDigi.RawToDigi.printUnknownFrameSummary = 0

process.load("EventFilter.TotemRawToDigi.totemTriggerRawToDigi_cfi")
process.totemTriggerRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")
process.totemTriggerRawToDigi.fedId = 0x29c

process.load('EventFilter.TotemRawToDigi.totemRPRawToDigi_cfi')
process.totemRPRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")
process.totemRPRawToDigi.RawToDigi.testID = 0 # TODO: remove this line once using non-emulated data
process.totemRPRawToDigi.fedIds = cms.vuint32(578, 579, 580) # TODO: remove this line once all TOTEM FED are used

# RP geometry
process.load("Geometry.VeryForwardGeometry.geometryRP_cfi")
process.XMLIdealGeometryESSource_CTPPS.geomXMLFiles.append("Geometry/VeryForwardData/data/RP_Garage/RP_Dist_Beam_Cent.xml")

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.TotemRPLocal.totemRPLocalReconstruction_cff")

# DQM Modules
process.load("DQM.Totem.totemDAQTriggerDQMSource_cfi")
process.load("DQM.Totem.totemRPDQMSource_cfi")

# processing path
process.recoStep = cms.Sequence(
  process.totemTriggerRawToDigi *
  process.totemRPRawToDigi *
  process.totemRPLocalReconstruction
)

process.dqmModules = cms.Sequence(
  process.totemDAQTriggerDQMSource +
  process.totemRPDQMSource +
  process.dqmEnv +
  process.dqmSaver
)

process.path = cms.Path(
  process.recoStep *
  #  process.dqmProvInfo * # TODO: needed?
  process.dqmModules
)

process.schedule = cms.Schedule(process.path)

process.dqmProvInfo.runType = process.runType.getRunTypeName()

# Process customizations included here
from DQM.Integration.config.online_customizations_cfi import *
process = customise(process)
