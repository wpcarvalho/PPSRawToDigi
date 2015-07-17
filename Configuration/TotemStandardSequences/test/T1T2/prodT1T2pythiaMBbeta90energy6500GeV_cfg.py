import FWCore.ParameterSet.Config as cms

from Configuration.TotemStandardSequences.prodT1T2Default_cfg import *

process.setName_("prodT1T2pythiaMBbeta90energy6500GeV")

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(10)
)

# Specify the output filename
exec 'process.' + str(process.outpath) + '.fileName = cms.untracked.string("file:prodT1T2pythiaMBbeta90energy6500GeV.root")'

# Pythia source labeled as "process.generator" @ CMSSW_3_1_1
process.load("Configuration.TotemCommon.PythiaMB_cfi")
process.generator.comEnergy = cms.double(6500)

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_6500GeV_90_cfi")

process.p1 = cms.Path(process.generator*process.SmearingGenerator*process.g4SimHits*process.mix*process.T1Digis*process.t1cluster*process.t1rechit*process.t1roads*process.t1tracks2*process.T2Digis*process.T2MCl*process.T2Hits*process.T2RoadPadFinder*process.T2TrackColl3*process.T2CC)
