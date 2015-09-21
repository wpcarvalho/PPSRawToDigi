import FWCore.ParameterSet.Config as cms

RPHitDists = cms.EDAnalyzer("HitDistributions",
    verbosity = cms.untracked.uint32(0),
    outputFile = cms.string('RPHitDists.root'),
    validTracksOnly = cms.untracked.uint32(1),
    rpList = cms.vuint32(100, 101, 102, 103, 104, 
        105, 120, 121, 122, 123, 
        124, 125, 0, 1, 2, 
        3, 4, 5, 20, 21, 
        22, 23, 24, 25),
    RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit")
)


