import FWCore.ParameterSet.Config as cms

process = cms.Process("T2CCSimulation2")

process.load("Configuration.TotemCommon.LoggerMax_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:input.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.Products = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:output-t2simu.root')
)

process.load("L1TriggerTotem.CoincidenceChip.T2CoincidenceProducer_cfi")
process.T2CC.coincidenceChipConfig.V = 7  # number of planes in "V out of NP" block
process.T2CC.verbose=3

process.producer = cms.Path(process.T2CC)

process.end = cms.EndPath(process.Products)


