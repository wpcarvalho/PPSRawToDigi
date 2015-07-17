import FWCore.ParameterSet.Config as cms

process = cms.Process("SmearingGeneratorTest")

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(10)
)

# Configure the output module
process.o1 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:smearing_generator_output.root')
)

# process.load("Configuration.TotemCommon.LoggerMax_cfi")

# process.load("SimGeneral.HepPDTESSource.pdt_cfi")

##################### STEP 1
process.source = cms.Source("EmptySource")

##################### STEP 2 - process.generator

# Use random number generator service
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# Monte Carlo gun - elastic specific
energy = "7000"
import IOMC.Elegent.ElegentSource_cfi
process.generator = IOMC.Elegent.ElegentSource_cfi.generator
process.generator.t_min = '6E-2'
process.generator.t_max = '6E-1'
process.generator.fileName = IOMC.Elegent.ElegentSource_cfi.ElegentDefaultFileName(energy)

##################### STEP 3 - process.SmearingGenerator

process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")

# optics
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_1180GeV_11_cfi")


#    # G4 geometry
#    process.load("Configuration.TotemCommon.geometryRP_cfi")
#    process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_1180_Beta_11_220/RP_Dist_Beam_Cent.xml')

#    process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup

#    # Elastic reconstruction options
#    process.load("RecoTotemRP.RPElasticReconstruction.rper_1180GeV_11_cfi")

process.SmearingGenerator.verbosity = 7
process.p1 = cms.Path(
                    process.generator
                    *process.SmearingGenerator
                      )

process.outpath = cms.EndPath(process.o1)

# print process.dumpConfig()