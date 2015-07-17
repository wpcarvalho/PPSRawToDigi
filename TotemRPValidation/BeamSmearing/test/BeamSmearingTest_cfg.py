import FWCore.ParameterSet.Config as cms

process = cms.Process("example")

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Configure the output module (save the result in a file)
process.o1 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:output.root')
)

process.load("Configuration.TotemCommon.LoggerMax_cfi")

# Use random number generator service
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# Use particle table
process.load("SimGeneral.HepPDTESSource.pdt_cfi")

# declare optics parameters
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_90_cfi")

# Smearing
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")
process.SmearingGenerator.originalLabel = 'original'
process.SmearingGenerator.modifyLabel = 'generator'
process.SmearingGenerator.verbosity = 0

################## STEP 1
process.source = cms.Source("EmptySource")

# particle producer
process.load("IOMC.FlatProtonLogKsiLogTGun.Beta90_cfi")
# process.FlatProtonLogKsiLogTGun.Verbosity = 1

# smearing validation
# TODO it should be possible to specify label in the configuration file
# currently it works only with IOMC.FlatProtonLogKsiLogTGun
process.load("TotemRPValidation.BeamSmearing.SmearingValidation_cfi")

process.p1 = cms.Path(process.generator*process.SmearingGenerator*process.SmearingValidation)

process.outpath = cms.EndPath(process.o1)

