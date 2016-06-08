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

# load DQM framework
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")

process.dqmEnv.subSystemFolder = 'CTPPS'

process.dqmSaver.convention = 'Offline'
process.dqmSaver.workflow = '/CTPPS/myTest/DQM'

# raw data source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/j/jkaspar/public/run273062_ls0001-2_stream.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# raw-to-digi conversion
process.load("CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi")
process.TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/ctpps_210_mapping.xml")

process.load("EventFilter.TotemRawToDigi.totemTriggerRawToDigi_cfi")
process.totemTriggerRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")

process.load("EventFilter.TotemRawToDigi.totemRPRawToDigi_cfi")
process.totemRPRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")

process.totemRawToDigi = cms.Sequence(
  process.totemTriggerRawToDigi *
  process.totemRPRawToDigi
)

# RP geometry
process.load("Geometry.VeryForwardGeometry.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/VeryForwardData/data/RP_Garage/RP_Dist_Beam_Cent.xml")

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.TotemRPLocal.totemRPLocalReconstruction_cff")

# CTPPS DQM modules
process.load("DQM.CTPPS.totemDAQTriggerDQMSource_cfi")
process.load("DQM.CTPPS.totemRPDQMSource_cfi")
process.load("DQM.CTPPS.totemRPDQMHarvester_cfi")

process.dqmoffline_step = cms.Path(
  process.totemRawToDigi *
  process.totemRPLocalReconstruction *

  process.totemDAQTriggerDQMSource *
  process.totemRPDQMSource *
  process.totemRPDQMHarvester *

  process.dqmSaver
)
