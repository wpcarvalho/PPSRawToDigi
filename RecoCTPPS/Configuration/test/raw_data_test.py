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
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/data/Run2016B/JetHT/RAW/v2/000/273/725/00000/9817A8A7-4F1E-E611-9248-02163E011EB4.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# raw-to-digi conversion
process.load("EventFilter.TotemRawToDigi.totemRawToDigi_cff")

# RP reconstruction chain with standard settings
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")

process.p = cms.Path(
    process.totemTriggerRawToDigi *
    process.totemRPRawToDigi *
    process.recoCTPPS
)

# output configuration
from RecoCTPPS.Configuration.RecoCTPPS_EventContent_cff import RecoCTPPSAOD
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:./AOD.root"),
    outputCommands = RecoCTPPSAOD.outputCommands
)

process.outpath = cms.EndPath(process.output)
