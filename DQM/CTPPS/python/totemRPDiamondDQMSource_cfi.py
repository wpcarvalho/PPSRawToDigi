import FWCore.ParameterSet.Config as cms

totemRPDiamondDQMSource = cms.EDAnalyzer("TotemRPDiamondDQMSource",
    tagStatus = cms.InputTag("diamondRPRawToDigi", "RP"),
    tagDigi = cms.InputTag("diamondRPRawToDigi", "RP"),
    tagLocalTrack = cms.InputTag("totemRPLocalTrackFitter"),
  
    verbosity = cms.untracked.uint32(0),
)
