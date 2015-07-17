import FWCore.ParameterSet.Config as cms

process = cms.Process("valT1T2pi")

########################### COMMON PART ##########################################

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/exp/totem/scratch/berretti/vers1cand3/CMSSW_1_7_7/src/RecoTotemT1T2/T2TrackProducer/test/gunT1T2pi.root')
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

process.load("TotemT1T2Validation.T2RecoValidation.T2RecoValidation_cfi")
process.T2RecoAn.OutputFile = 'T2RecoValv1c3.root'

process.p1 = cms.Path(process.T2RecoAn)
