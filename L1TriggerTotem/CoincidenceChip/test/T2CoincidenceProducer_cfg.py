import FWCore.ParameterSet.Config as cms

process = cms.Process("T2CCSimulation")

process.load("Configuration.TotemCommon.LoggerMax_cfi")

process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:input.root')
     fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/totem/offline/Reco/2010/T2/run_2136_NoClean_IntGlobAL24666_.root')

)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

process.Products = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:output-t2simu-10kEvts-r2136-quarterPF-q2-logPrintItAll.root')
)

process.load("L1TriggerTotem.CoincidenceChip.T2CoincidenceProducer_cfi")
process.T2CC.coincidenceChipConfig.useControlRegisters = False  # if True then the control numbers V,NP,OV... are not taken into account
process.T2CC.quarter = cms.uint32(2)
process.T2CC.coincidenceChipConfig.V = 5   # number of planes in "V out of NP" block
process.T2CC.coincidenceChipConfig.Z = 12   # number of planes in "Z out of 8/16" block -- change zEff too, in plugins/T2*.cc
process.T2CC.coincidenceChipConfig.useLogicWithWrongNP=1

process.producer = cms.Path(process.T2CC)

process.end = cms.EndPath(process.Products)


