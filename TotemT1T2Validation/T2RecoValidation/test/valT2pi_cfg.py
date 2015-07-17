import FWCore.ParameterSet.Config as cms

process = cms.Process("valT1T2pi")

########################### COMMON PART ##########################################

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/totem/offline/berretti/Reco50000Pi5070E10rs51v1c3.root')
)

# logging to txt files 
process.load("Configuration.TotemCommon.LoggerMax_cfi")

# random number generators seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

########################### T1 VALIDATION ##########################################

process.T1Val = cms.EDAnalyzer("T1Validation",
    SIM = cms.double(1.0),
    DIGI = cms.double(1.0),
    RECO = cms.double(1.0),
    TrackLabel = cms.string('t1tracks'),
    OutputFile = cms.string('valT1PlotsPi.root')
)

########################### T2 VALIDATION ##########################################

process.T2G4An = cms.EDAnalyzer("T2G4Analyzer",
    t2SimProdLabel = cms.string('g4SimHits'),
    t2SimInstLabel = cms.string('TotemHitsT2Gem'),
    OutputFile = cms.untracked.string('valT2PlotsG4Pi.root')
)

process.T2DigiAn = cms.EDAnalyzer("T2DigiAnalyzer",
    t2DigiProdLabel = cms.string('T2Digis'),
    t2DigiInstLabel = cms.string('T2StripDigi'),
    t2DigiInstLabelP = cms.string('T2PadDigi'),
    OutputFile = cms.untracked.string('valT2PlotsDigiPi.root')
)

process.T2RecoAn = cms.EDAnalyzer("T2RecoAnalyzer",
    FitdzEnergy = cms.double(10.0),         # choose 10. or 50. GeV
    Chicut = cms.double(3.0),
    DZScale = cms.double(1.0),
    TrkEtamin = cms.double(5.1),            # 4.9 for 10 GeV
    TrkEtaMAX = cms.double(6.9),            # 7.0 for 10 GeV

    PhiChiProbCut = cms.double(0.02),
    RChiProbCut = cms.double(0.02),

    CluLabel = cms.string('T2MCl'),
    HitLabel = cms.string('T2Hits'),
    RoadLabel = cms.string('T2RoadColl'),
    TrackLabel = cms.string('T2TrackColl2'),
    OutputFile = cms.untracked.string('valT2PlotsRecoPirs51.root')
)

process.p1 = cms.Path(process.T2RecoAn)

