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

process.dqmEnv.subSystemFolder = 'Totem'
process.dqmSaver.workflow = '/Totem/Test/Workflow'

process.dqmSaver = cms.EDAnalyzer("DQMFileSaverOnline",
  producer = cms.untracked.string("DQM"),
  tag = cms.untracked.string("Totem"),
  path = cms.untracked.string("."),
)

# RP raw data and digi
process.load("DQM.Totem.test.standaloneDataInput_cff")
#process.load("DQM.Totem.test.emulatedDataInput_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# RP geometry
process.load("Geometry.VeryForwardGeometry.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/VeryForwardData/data/RP_Garage/RP_Dist_Beam_Cent.xml")

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.TotemRPLocal.totemRPLocalReconstruction_cff")

# TOTEM DQM modules
process.load("DQM.Totem.totemDAQTriggerDQMSource_cfi")
process.load("DQM.Totem.totemRPDQMSource_cfi")

# execution schedule
process.reco_totem = cms.Path(
  process.totemRawToDigi *
  process.totemRPLocalReconstruction
)

process.dqm_totem = cms.Path(
  process.totemDAQTriggerDQMSource *
  process.totemRPDQMSource
)

process.dqm_common = cms.Path(
    process.dqmEnv *
    process.dqmSaver
    #process.dqmStoreStats
)

process.schedule = cms.Schedule(
    process.reco_totem,
    process.dqm_totem,
    process.dqm_common
)
