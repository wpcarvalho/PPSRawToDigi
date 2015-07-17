import FWCore.ParameterSet.Config as cms

process = cms.Process("ValidationPlots")

process.load("Configuration.TotemCommon.LoggerMin_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:elastic_90_smeared_sim.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.ValPlots = cms.EDAnalyzer("ValidationPlots",
    verbosity = cms.untracked.uint32(0),
    outputFile = cms.string('plots-el-90.root'),
    type = cms.string('elastic'),
    smearedSource = cms.untracked.bool(True),
    beta = cms.untracked.uint32(90),
    rpList = cms.vuint32(100, 101, 102, 103, 104, 
        105, 120, 121, 122, 123, 
        124, 125, 0, 1, 2, 
        3, 4, 5, 20, 21, 
        22, 23, 24, 25)
)

process.p = cms.Path(process.ValPlots)

