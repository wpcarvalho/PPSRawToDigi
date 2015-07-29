import FWCore.ParameterSet.Config as cms

process = cms.Process("TDTestElastic")

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Configure the output module (save the result in a file)
process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('file:test.root')
)
process.outpath = cms.EndPath(process.o1)

# Configure if you want to detail or simple log information.
# LoggerMax -- detail log info output including: errors.log, warnings.log, infos.log, debugs.log
# LoggerMin -- simple log info output to the standard output (e.g. screen)
process.load("Configuration.TotemCommon.LoggerMax_cfi")

################## STEP 1 - process.generator
process.source = cms.Source("EmptySource")

# Use random number generator service
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# Monte Carlo gun - elastic specific
energy = "6500"
import IOMC.Elegent.ElegentSource_cfi
process.generator = IOMC.Elegent.ElegentSource_cfi.generator
process.generator.fileName = IOMC.Elegent.ElegentSource_cfi.ElegentDefaultFileName(energy)

################## STEP 2 process.SmearingGenerator

# declare optics parameters
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_6500GeV_90_cfi")

# Smearing
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")

################## STEP 3 process.g4SimHits

# Geometry - beta* specific
process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90/RP_Dist_Beam_Cent.xml')

# Magnetic Field, by default we have 3.8T
process.load("Configuration.StandardSequences.MagneticField_cff")

# G4 simulation & proton transport
process.load("Configuration.TotemCommon.g4SimHits_cfi")
process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup
process.g4SimHits.Generator.HepMCProductLabel = 'generator'    # The input source for G4 module is connected to "process.source".
process.g4SimHits.G4TrackingManagerVerbosity = cms.untracked.int32(3)
process.g4SimHits.UseMagneticField = cms.bool(False) # todo enable magnetic field

# Use particle table
process.load("SimGeneral.HepPDTESSource.pdt_cfi")

################## STEP 4 mix pdt_cfi

#process.load("Configuration.TotemStandardSequences/RP_Digi_and_TrackReconstruction_cfi")
#process.load("Configuration.TotemCommon/mixNoPU_cfi")


process.p1 = cms.Path(
	process.generator
	*process.SmearingGenerator
	*process.g4SimHits
)


