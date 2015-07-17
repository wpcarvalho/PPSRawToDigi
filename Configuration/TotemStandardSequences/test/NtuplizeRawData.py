import FWCore.ParameterSet.Config as cms

process = cms.Process("nTuplize")

process.maxEvents = cms.untracked.PSet(
input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    noEventSort = cms.untracked.bool(True),
    skipBadFiles = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)
process.source.fileNames.append('file:./rawDataReco.root')
# process.source.fileNames.append('rfio:/castor/cern.ch/user/p/psikora/parallel_output/totcsi_8787.001.root')
# process.source.fileNames.append('rfio:/castor/cern.ch/user/p/psikora/parallel_output/totcsi_8787.002.root')
# process.source.fileNames.append('rfio:/castor/cern.ch/user/p/psikora/parallel_output/totcsi_8787.003.root')
# process.source.fileNames.append('rfio:/castor/cern.ch/user/p/psikora/parallel_output/totcsi_8787.004.root')
# process.source.fileNames.append('rfio:/castor/cern.ch/user/p/psikora/parallel_output/totcsi_8787.005.root')
# process.source.fileNames.append('rfio:/castor/cern.ch/user/p/psikora/parallel_output/totcsi_8787.006.root')


process.load("Configuration.TotemCommon.LoggerMin_cfi")
process.load("Configuration/TotemOpticsConfiguration/OpticsConfig_4000GeV_90_cfi")

process.load("TotemAnalysis/TotemNtuplizer/TotemNtuplizer_cfi")
process.TotemNtuplizer.verbosity = 10
process.TotemNtuplizer.outputFileName = 'file:./NtuplizeRawData.ntuple.root'
process.TotemNtuplizer.ProductLabelSimu = 'rpCCOutput'
process.TotemNtuplizer.ModulLabelSimu = 'Raw2DigiProducer'

process.TotemNtuplizer.RPStripDigiSetLabel = cms.InputTag('Raw2DigiProducer')
process.TotemNtuplizer.RawEventLabel = cms.InputTag('Raw2DigiProducer')
process.TotemNtuplizer.RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit")
process.TotemNtuplizer.RPMulFittedTrackCollectionLabel = cms.InputTag("RPMulTrackNonParallelCandCollFit")
process.TotemNtuplizer.RPReconstructedProtonCollectionLabel = cms.InputTag("RP220Reconst")
process.TotemNtuplizer.RPReconstructedProtonPairCollectionLabel = cms.InputTag("RP2202ArmReconst")

process.p = cms.Path(
     process.TotemNtuplizer
)
