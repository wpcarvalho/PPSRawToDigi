import FWCore.ParameterSet.Config as cms

process = cms.Process("GenerateES")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# random seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")
process.RandomNumberGeneratorService.generator.initialSeed = 123
process.RandomNumberGeneratorService.SmearingGenerator.initialSeed = 234

# set number of events
process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)

# physics process generator - Single Diffraction
process.load("Configuration.TotemCommon.PythiaSD_cfi")
process.generator.comEnergy = 8000

# beam smearing
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")
process.SmearingGenerator.verbosity = 1

# load common settings (optics, geometry)
process.load("common_cfi")

# geometry printer
process.GeometryInfo = cms.EDAnalyzer("GeometryInfoModule")

# fast simulation + proton dumper
process.load("TotemAlignment.RPFastSimulation.RPFastFullSimulation200_cfi")
process.RPFastFullSimulation200.verbosity = 5
process.RPFastFullSimulation200.plotFileName = "plots_SD.root"
process.RPFastFullSimulation200.dumpTopRPProtons = True

process.p = cms.Path(
    #process.GeometryInfo *
    process.generator *
    process.SmearingGenerator *
    process.RPFastFullSimulation200
)
