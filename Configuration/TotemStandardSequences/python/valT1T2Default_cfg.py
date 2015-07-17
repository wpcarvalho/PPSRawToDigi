import FWCore.ParameterSet.Config as cms

process = cms.Process("valT1T2Default")

# logging to txt files
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# Configure the input module (read the events from a file)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:prodT1T2pythiaDDbeta11energy1.1TeV.root'),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    skipBadFiles = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)     # -1 means to take all the events from the source as input
)

# random number generators seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# Geometry of T1T2
process.load("Configuration.TotemCommon.geometryT1T2_cfi")

#process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_1180GeV_11_cfi")

######################## T1 Validation modules 

from Configuration.TotemStandardSequences.T1_Validation_cfi import *

process.T1Val = t1valid
process.T1Val.OutputFile = cms.string("file:valT1T2pythiaDDbeta11energy1.1TeV_T1Val.root")


######################## T2 Validation modules 

from Configuration.TotemStandardSequences.T2_Validation_cfi import *

process.T2RecoAn = T2RecoAn
process.T2RecoAn.HepMCProductLabel = cms.string("generator")
process.T2RecoAn.OutputFile = "file:valT1T2pythiaDDbeta11energy1.1TeV_T2RecoAn.root" 

process.load("TotemAnalysis.TotemNtuplizer.TotemNtuplizer_cfi")
process.TotemNtuplizer.outputFileName = 'file:valRPT1T2Default_ntuple.root'
process.TotemNtuplizer.verbosity = 0
process.TotemNtuplizer.ProductLabelSimu = 'RPCCRawBits'
process.TotemNtuplizer.ModulLabelSimu = 'RPDataCCProducer'
process.TotemNtuplizer.T1TrackLabel = 't1tracks2'

#process.p1 = cms.Path(process.T1Val*process.T2RecoAn)
process.p1 = cms.Path(process.T2RecoAn*process.TotemNtuplizer)
