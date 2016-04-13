import FWCore.ParameterSet.Config as cms

TotemTriggerRawToDigi = cms.EDProducer("TotemTriggerRawToDigi",
  rawDataTag = cms.InputTag(""),

  # TODO: remove 577
  # IMPORTANT: leave empty to load the default configuration from
  #    DataFormats/FEDRawData/interface/FEDNumbering.h
  fedId = cms.uint32(577)
)
