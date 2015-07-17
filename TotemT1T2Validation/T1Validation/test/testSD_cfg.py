import FWCore.ParameterSet.Config as cms

process = cms.Process("trackT1T2SD")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('rfio:simu-SD.root')
)

# logging to txt files
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# random number generators seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

process.load("SimTotem.T1Digitizer.T1DigisVFAT_cfi")

process.load("RecoTotemT1T2.T1MakeCluster.T1MakeCluster_cfi")

process.T2MCl = cms.EDFilter("T2MakeCluster",
    Verbosity = cms.int32(0)
)

process.load("RecoTotemT1T2.T1RecHit.T1RecHit_cfi")

process.T2Hits = cms.EDFilter("T2RecHit")

process.load("RecoTotemT1T2.T1RoadProducer.T1RoadProducer_cfi")

process.load("RecoTotemT1T2.T2RoadProducer.T2RoadCollDr15Dz200Dphi17_cfi")

process.load("RecoTotemT1T2.T1TrackProducer.T1TrackProducer_cfi")

process.load("SimTotem.T2Digitizer.T2Digis_cfi")

process.load("RecoTotemT1T2.T2TrackProducer2.T2TrackColl2_cfi")

process.load("RecoTotemT1T2.T1TrackProducer.T1TrackAnalyzer_cfi")

process.load("TotemT1T2Validation.T1Validation.T1Validation_cfi")

process.o1 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('valid-SD_.root')
)

process.p1 = cms.Path(process.T1Digis*process.T2Digis*process.t1cluster*process.T2MCl*process.t1rechit*process.T2Hits*process.T2RoadColl*process.T2TrackColl2*process.t1roads*process.t1tracks*process.t1valid)

process.outpath = cms.EndPath(process.o1)
