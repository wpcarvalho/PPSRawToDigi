import FWCore.ParameterSet.Config as cms

process = cms.Process("ValidationPlots")

process.load("Configuration.TotemCommon.LoggerMin_cfi")

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.HitDists1 = cms.EDAnalyzer("HitDistributions",
    verbosity = cms.untracked.uint32(2),
    outputFile = cms.string('/tmp/1.root'),
    rpList = cms.vuint32(123, 124, 125)
)

process.HitDists2 = cms.EDAnalyzer("HitDistributions",
    verbosity = cms.untracked.uint32(2),
    outputFile = cms.string('/tmp/2.root'),
    rpList = cms.vuint32(120, 121, 122)
)

process.HitDists3 = cms.EDAnalyzer("HitDistributions",
    verbosity = cms.untracked.uint32(2),
    outputFile = cms.string('/tmp/3.root'),
    rpList = cms.vuint32(20, 21, 22)
)

process.HitDists4 = cms.EDAnalyzer("HitDistributions",
    verbosity = cms.untracked.uint32(2),
    outputFile = cms.string('/tmp/4.root'),
    rpList = cms.vuint32(23, 24, 25)
)

process.ValPlots = cms.EDAnalyzer("ValidationPlots",
    verbosity = cms.untracked.uint32(0),
    hitsFileList = cms.vstring('/tmp/1.root', 
        '/tmp/2.root', 
        '/tmp/3.root', 
        '/tmp/4.root'),
    resFile = cms.string('/tmp/DPE_right_220_no_vertex_90_smeared_hist.root'),
    outputFile = cms.string('/tmp/dpe.root'),
    beta = cms.uint32(90),
    type = cms.string('inelastic')
)

process.p = cms.Path(process.ValPlots)

