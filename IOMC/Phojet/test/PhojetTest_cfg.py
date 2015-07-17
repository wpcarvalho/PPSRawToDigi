import FWCore.ParameterSet.Config as cms

process = cms.Process("PhojetTest")

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

# Configure the output module (save the result in a file)
process.o1 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:phojet_output.root')
)

process.load("Configuration.TotemCommon.LoggerMax_cfi")

# Use random number generator service
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# Use particle table
process.load("SimGeneral.HepPDTESSource.pdt_cfi")

################## STEP 1
process.source = cms.Source("EmptySource")

################## STEP 2 - process.generator

# Use random number generator service
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# Monte Carlo gun
import IOMC.Phojet.Phojet_cfi
process.generator = IOMC.Phojet.Phojet_cfi.generator
process.generator.verbosity = 0
process.generator.bufferSize = 100
process.generator.cmsEnergy = 14000.
#process.generator.process = "SD"  #specified from the list  
process.generator.process = "0 0 0 0 1 1 0 1" #specified manually

process.p1 = cms.Path(process.generator)

process.outpath = cms.EndPath(process.o1)
