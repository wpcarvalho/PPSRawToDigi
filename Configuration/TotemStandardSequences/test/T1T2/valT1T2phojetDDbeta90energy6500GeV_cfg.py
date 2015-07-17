import FWCore.ParameterSet.Config as cms

from Configuration.TotemStandardSequences.valT1T2Default_cfg import *

process.setName_("valT1T2phojetDDbeta90energy6500GeV")

process.source.fileNames = cms.untracked.vstring('file:prodT1T2phojetDDbeta90energy6500GeV.root')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)     # -1 means to take all the events from the source as input
)

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_6500GeV_90_cfi")

process.T1Val.OutputFile = cms.string("file:valT1T2phojetDDbeta90energy6500GeV_T1Val.root")

process.T2RecoAn.HepMCProductLabel = cms.string("generator")
process.T2RecoAn.OutputFile = "file:valT1T2phojetDDbeta90energy6500GeV_T2RecoAn.root"
