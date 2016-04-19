import FWCore.ParameterSet.Config as cms

TotemDAQTriggerDQMSource = cms.EDAnalyzer("TotemDAQTriggerDQMSource",
  tagFEDInfo = cms.InputTag("TotemRPRawToDigi", "RP"),
  tagTriggerCounters = cms.InputTag("TotemTriggerRawToDigi")
)
