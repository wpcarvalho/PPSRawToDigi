import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLETOTEM")

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring("file:MinbiasTest_RECO_TOTEM.root"),
        noEventSort = cms.untracked.bool(True)
   )

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)


process.load("TotemAnalysis.TotemNtuplizer.TotemNtuplizer_cfi")
process.TotemNtuplizer.outputFileName = 'file:MinbiasTest_NTUPLE_TOTEM.root'
process.TotemNtuplizer.verbosity = 2
process.TotemNtuplizer.ProductLabelSimu = 'RPCCRawBits'
process.TotemNtuplizer.ModulLabelSimu = 'RPDataCCProducer'
process.TotemNtuplizer.T1TrackLabel = 't1tracks2'

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_6500GeV_90_cfi")

# Module chain
process.p1 = cms.Path(process.TotemNtuplizer)

#print process.dumpPython()
