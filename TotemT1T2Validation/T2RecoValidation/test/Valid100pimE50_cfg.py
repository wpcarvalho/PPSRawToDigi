import FWCore.ParameterSet.Config as cms

process = cms.Process("T2Reco")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:100pim_bp_E50_phiFlat_etaFlat_B.root')
)

# keep the logging output to a nice level
process.MessageLogger = cms.Service("MessageLogger")

process.analyzerReco = cms.EDAnalyzer("T2RecoAnalyzer",
    trkZRange = cms.double(2800.0),
    chicut = cms.double(1.0),
    e50range1m = cms.double(5.35),
    e50range1M = cms.double(6.13),
    e50range2m = cms.double(6.13),
    e50range2M = cms.double(6.26),

    e10range1m = cms.double(5.4),
    e10range1M = cms.double(6.18),
    typeDzdistr = cms.double(1.0),
    numsigmaZ = cms.double(3.0),
    energy = cms.double(50.0),      # choose 10. or 50. GeV
    OutputFile = cms.untracked.string('T2Reco.root')
)

process.p = cms.Path(process.analyzerReco)

