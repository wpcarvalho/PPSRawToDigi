import FWCore.ParameterSet.Config as cms

process = cms.Process("DQM")

#----------------------------
#### Event Source
#----------------------------
# for live online DQM in P5
process.load("DQM.Integration.config.inputsource_cfi")

# for testing in lxplus
#process.load("DQM.Integration.config.fileinputsource_cfi")

#----------------------------
#### DQM Environment
#----------------------------
process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = 'Info'
process.dqmSaver.tag = 'Info'
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
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/totem_rp_210far_220_mapping.xml")

# process.load('EventFilter.TotemRawToDigi.TotemRPRawToDigi_cfi')
# process.totemRPRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")
# process.totemRPRawToDigi.fedIds = cms.vuint32(577, 578, 579, 580)
# process.totemRPRawToDigi.RawToDigi.printErrorSummary = 0
# process.totemRPRawToDigi.RawToDigi.printUnknownFrameSummary = 0

process.load("EventFilter.TotemRawToDigi.totemTriggerRawToDigi_cfi")
process.totemTriggerRawToDigi.rawDataTag = cms.InputTag("source")
process.totemTriggerRawToDigi.fedId = 0x29c

process.load('EventFilter.TotemRawToDigi.totemRPRawToDigi_cfi')
process.totemRPRawToDigi.rawDataTag = cms.InputTag("source")
process.totemRPRawToDigi.fedIds = cms.vuint32(0x1a1, 0x1a2, 0x1a9, 0x1aa, 0x1b5, 0x1bd)

# RP geometry
process.load("Geometry.VeryForwardGeometry.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/VeryForwardData/data/RP_Garage/RP_Dist_Beam_Cent.xml")

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

process.dqmModules = cms.Sequence(process.dqmEnv + process.dqmSaver)

process.path = cms.Path(
  process.recoStep *
  process.dqmProvInfo *
  process.dqmModules
)

process.schedule = cms.Schedule(process.path)

process.dqmProvInfo.runType = process.runType.getRunTypeName()

# Process customizations included here
from DQM.Integration.config.online_customizations_cfi import *
process = customise(process)
