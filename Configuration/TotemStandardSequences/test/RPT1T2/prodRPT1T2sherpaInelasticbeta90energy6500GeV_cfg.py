import FWCore.ParameterSet.Config as cms

process = cms.Process("prodT1T2RPDefault")
process.setName_("prodRPT1T2sherpaInelasticbeta90energy6500GeV")

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(10)
)
# Configure the output module (save the result in a file)
process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *', 'drop *_*mix*_*_*',  'drop *_*_*TrackerHits*_*', 'drop *_*_*Muon*_*', 'drop *_*_*Ecal*_*', 'drop *_*_*Hcal*_*', 'drop *_*_*Calo*_*', 'drop *_*_*Castor*_*', 'drop *_*_*FP420SI_*', 'drop *_*_*ZDCHITS_*', 'drop *_*_*BSCHits_*', 'drop *_*_*ChamberHits_*', 'drop *_*_*FibreHits_*', 'drop *_*_*WedgeHits_*'),
    fileName = cms.untracked.string('file:prodT1T2RPDefault.root')
)
process.outpath = cms.EndPath(process.o1)

process.source = cms.Source("EmptySource")

# Configure if you want to detail or simple log information.
# LoggerMax -- detail log info output including: errors.log, warnings.log, infos.log, debugs.log
# LoggerMin -- simple log info output to the standard output (e.g. screen)
process.load("Configuration.TotemCommon.LoggerMin_cfi")




# Specify the output filename
exec 'process.' + str(process.outpath) + '.fileName = cms.untracked.string("file:prodRPT1T2sherpaInelasticbeta90energy6500GeV.root")'

# Pythia source labeled as "process.generator" @ CMSSW_3_1_1
process.load("Configuration.TotemCommon.SherpaInelastic_cfi")

# Sherpa needs a temporary directory, which is specified in
# Configuration.TotemCommon.SherpaInelastic_cfi.
# Create the directory, if it doesn't already exist.
tmpDir = process.generator.libDir._value
import os
if not os.path.exists(tmpDir):
    os.makedirs(tmpDir)

# G4 geometry
process.load("Configuration.TotemCommon.geometryGlobal_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90/RP_Dist_Beam_Cent.xml')

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_6500GeV_90_cfi")

# process.g4SimHits.Generator.HepMCProductLabel = 'generator'
# process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup
#
# process.load("RecoTotemRP.RPInelasticReconstruction.Rec_6500GeV_beta_90p0_220_cfi")
# process.RP220Reconst.BeamProtTransportSetup = process.BeamProtTransportSetup
#
# process.g4SimHits.Generator.HepMCProductLabel = 'generator'
# process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup

process.p1 = cms.Path(process.generator)

