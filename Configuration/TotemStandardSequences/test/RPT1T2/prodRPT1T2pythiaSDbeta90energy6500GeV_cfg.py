import FWCore.ParameterSet.Config as cms

from Configuration.TotemStandardSequences.prodRPT1T2Default_cfg import *

process.setName_("prodRPT1T2pythiaSDbeta90energy6500GeV")

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(10)
)

# Specify the output filename
exec 'process.' + str(process.outpath) + '.fileName = cms.untracked.string("file:prodRPT1T2pythiaSDbeta90energy6500GeV.root")'

# Pythia source labeled as "process.generator" @ CMSSW_3_1_1
process.load("Configuration.TotemCommon.PythiaSD_cfi")
process.generator.comEnergy = cms.double(6500)

process.load("Configuration.TotemCommon.geometryGlobal_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90/RP_Dist_Beam_Cent.xml')

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_6500GeV_90_cfi")

process.g4SimHits.Generator.HepMCProductLabel = 'generator'
process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup

process.load("RecoTotemRP.RPInelasticReconstruction.Rec_6500GeV_beta_90_cfi")
process.RP220Reconst.BeamProtTransportSetup = process.BeamProtTransportSetup

process.p1 *= cms.Sequence(process.RP220Reconst)
