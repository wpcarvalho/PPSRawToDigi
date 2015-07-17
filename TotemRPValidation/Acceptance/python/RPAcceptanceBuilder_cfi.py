import FWCore.ParameterSet.Config as cms

RPAcceptanceBuilder = cms.EDAnalyzer("RPAcceptanceBuilder",
	verbosity = cms.untracked.uint32(0),
	hepMCLabel = cms.untracked.string("generator"),
	outputFileName = cms.string("RPAcceptanceBuilder.root"),
	RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit")
)
