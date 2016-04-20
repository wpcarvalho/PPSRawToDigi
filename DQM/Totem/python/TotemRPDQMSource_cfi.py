import FWCore.ParameterSet.Config as cms

TotemRPDQMSource = cms.EDAnalyzer("TotemRPDQMSource",
    tagStatus = cms.InputTag("TotemRPRawToDigi", "RP"),
    tagDigi = cms.InputTag("TotemRPRawToDigi", "RP"),
	tagCluster = cms.InputTag("TotemRPClusterProducer"),
	tagRecHit = cms.InputTag("TotemRPRecHitProducer"),
	tagUVPattern = cms.InputTag("TotemRPUVPatternFinder"),
	tagLocalTrack = cms.InputTag("TotemRPLocalTrackFitter"),

    buildCorrelationPlots = cms.untracked.bool(False),
    correlationPlotsLimit = cms.untracked.uint32(50),
	correlationPlotsFilter = cms.untracked.string("default=0,1")
)
