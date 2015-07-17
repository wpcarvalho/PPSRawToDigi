import FWCore.ParameterSet.Config as cms

process = cms.Process("RPCCSimulation")

process.load("Configuration.TotemCommon.LoggerMax_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:input2.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)
process.o1 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:output.root')
)

# RPCC module
process.load("L1TriggerTotem.CoincidenceChip.RPCoincidenceProducer_cfi")

process.p1 = cms.Path(process.RPCC)

process.outpath = cms.EndPath(process.o1) 
