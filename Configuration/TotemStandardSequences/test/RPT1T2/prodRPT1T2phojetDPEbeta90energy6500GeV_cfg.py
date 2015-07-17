import FWCore.ParameterSet.Config as cms

from Configuration.TotemStandardSequences.prodRPT1T2Default_cfg import *

process.setName_("prodRPT1T2phojetDPEbeta90energy6500GeV")

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(10)
)

# Specify the output filename
exec 'process.' + str(process.outpath) + '.fileName = cms.untracked.string("file:prodRPT1T2phojetDPEbeta90energy6500GeV.root")'

import IOMC.Phojet.Phojet_cfi
process.generator = IOMC.Phojet.Phojet_cfi.generator
process.generator.verbosity = 0
process.generator.bufferSize = 100
process.generator.cmsEnergy = 6500.0
process.generator.process = "DPE"

process.load("Configuration.TotemCommon.geometryGlobal_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90/RP_Dist_Beam_Cent.xml')

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_6500GeV_90_cfi")

process.g4SimHits.Generator.HepMCProductLabel = 'generator'
process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup

process.load("RecoTotemRP.RPInelasticReconstruction.Rec_6500GeV_beta_90_cfi")
process.RP220Reconst.BeamProtTransportSetup = process.BeamProtTransportSetup

#process.load("RecoTotemRP.RPInelasticReconstruction.Rec_6500GeV_beta_90_220_2Arm_cfi")
#process.RP2202ArmReconst.BeamProtTransportSetup = process.BeamProtTransportSetup

process.p1 *= cms.Sequence(process.RP220Reconst)
process.p1 *= cms.Sequence(process.RP2202ArmReconst)
