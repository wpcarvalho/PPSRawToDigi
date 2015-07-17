import FWCore.ParameterSet.Config as cms

process = cms.Process("valT1T2Default")

# logging to txt files
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# Configure the input module (read the events from a file)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:valT1T2.root'),
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

process.load("TotemT1T2Validation.T1Validation.T1Validation_cfi")

process.t1valid.OutputFile = cms.string("file:valT1T2pythiaDDbeta11energy1.1TeV_T1Val.root")


######################## T2 Validation modules 

#process.load("TotemT1T2Validation.T2DigiValidation.T2DigiValidation_cfi")

#process.load("TotemT1T2Validation.T2G4Validation.T2G4Validation_cfi")

#process.load("TotemT1T2Validation.T2RecoValidation.T2RecoValidation_cfi")

#process.load("TotemT1T2Validation.T2RecoValidation.T2RecoValidation2_cfi")

#process.T2DigiAn = T2DigiValidation
#process.T2DigiAn.OutputFile = cms.untracked.string("file:valT1T2pythiaDDbeta11energy1.1TeV_T2DigiAn.root")

#process.T2G4An = T2G4Validation
#process.T2G4An.OutputFile = cms.untracked.string("file:valT1T2pythiaDDbeta11energy1.1TeV_T2G4An.root")

#process.T2RecoAn = T2RecoAn
#process.T2RecoAn.HepMCProductLabel = cms.string("generator")
#process.T2RecoAn.OutputFile = "file:valT1T2pythiaDDbeta11energy1.1TeV_T2RecoAn.root" 

#process.T2RecoAn2 = T2RecoAn2
#process.T2RecoAn2.HepMCProductLabel = cms.string("generator")
#process.T2RecoAn2.OutputFile = "file:valT1T2pythiaDDbeta11energy1.1TeV_T2RecoAn2.root"

#process.load("TotemAnalysis/TotemNtuplizer/TotemNtuplizer_cfi")
#process.TotemNtuplizer.outputFileName = 'file:valRPT1T2Default_ntuple.root'
#process.TotemNtuplizer.verbosity = 0
#process.TotemNtuplizer.ProductLabelSimu = 'RPCCRawBits'
#process.TotemNtuplizer.ModulLabelSimu = 'RPDataCCProducer'
#process.TotemNtuplizer.T1TrackLabel = 't1tracks2'

#########################################################################

process.setName_("valT1T2Default")

process.source.fileNames = cms.untracked.vstring('file:valT1T2.root')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)     # -1 means to take all the events from the source as input
)

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_1180GeV_11_cfi")

process.t1valid.OutputFile = cms.string("file:valT1T2phojetMBbeta11energy1.1TeV_T1Val.root")

#process.T2DigiAn.OutputFile = cms.untracked.string("file:valT1T2phojetMBbeta11energy1.1TeV_T2DigiAn.root")

#process.T2G4An.OutputFile = cms.untracked.string("file:valT1T2phojetMBbeta11energy1.1TeV_T2G4An.root")

#process.T2RecoAn.HepMCProductLabel = cms.string("generator")
#process.T2RecoAn.OutputFile = "file:valT1T2phojetMBbeta11energy1.1TeV_T2RecoAn.root" 

#process.T2RecoAn2.HepMCProductLabel = cms.string("generator")
#process.T2RecoAn2.OutputFile = "file:valT1T2phojetMBbeta11energy1.1TeV_T2RecoAn2.root"

############################################################################

process.p1 = cms.Path(process.t1valid
#			*process.T2G4An
#			*process.T2DigiAn
#			*process.T2RecoAn
#			*process.T2RecoAn2
#			*process.TotemNtuplizer
			)

# T2RecoValidation T2DEVNtuplizer TotemNtuplizer

#process.outpath = cms.EndPath(process.o1)

