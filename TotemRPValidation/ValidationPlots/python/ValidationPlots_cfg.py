import FWCore.ParameterSet.Config as cms

process = cms.Process("ValidationPlots")

process.load("Configuration.TotemCommon.LoggerMin_cfi")

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.ValPlots = cms.EDAnalyzer("ValidationPlots",
    verbosity = cms.untracked.uint32(2),
    hitsFileList = cms.vstring('inelastic_90_smeared_hit_1.root', 
        'inelastic_90_smeared_hit_2.root'),
    resFile = cms.string('inelastic_90_smeared_val.root'),
    outputFile = cms.string('temp-plots.root'),
    beta = cms.uint32(90),
    type = cms.string('inelastic')
)

process.p = cms.Path(process.ValPlots)

