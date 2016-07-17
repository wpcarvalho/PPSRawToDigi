import FWCore.ParameterSet.Config as cms

process = cms.Process("TotemStandaloneRawDataTest")

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cout'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# raw data source
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/j/jkaspar/public/run273062_ls0001-2_stream.root')
#)
process.source = cms.Source("NewEventStreamFileReader",
    fileNames = cms.untracked.vstring('/store/t0streamer/Minidaq/A/000/276/395/run276395_ls0001_streamA_StorageManager.dat')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# raw-to-digi conversion
process.load("EventFilter.TotemRawToDigi.totemRawToDigi_cff")

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

# ntuplizer
process.load("TotemAnalysis.TotemNtuplizer.TotemNtuplizer_cfi")
process.totemNtuplizer.outputFileName = "ntuple.root"

process.p = cms.Path(
    process.totemTriggerRawToDigi *
    process.totemRPRawToDigi *
    process.recoCTPPS *
    process.totemNtuplizer
)

# output configuration
from RecoCTPPS.Configuration.RecoCTPPS_EventContent_cff import RecoCTPPSRECO
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:reco.root"),
    outputCommands = RecoCTPPSRECO.outputCommands
)

process.outpath = cms.EndPath(process.output)
