import FWCore.ParameterSet.Config as cms

process = cms.Process("DiamondDigiToRecoTest")
process.maxEvents = cms.untracked.PSet(
input = cms.untracked.int32(10)
)
# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cout'),
    cout = cms.untracked.PSet( threshold = cms.untracked.string('ERROR') )
#    destinations = cms.untracked.vstring('cerr'),
#    cerr = cms.untracked.PSet( threshold = cms.untracked.string('DEBUG') )
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:./DiamondDigi.root')
)

# digi-to-reco conversion
process.load('RecoCTPPS.DiamondLocal.diamondRecHitProducer_cfi')
# process.diamondRecHitProducer.tagDigi = cms.InputTag("diamondRPRawToDigi")

process.p = cms.Path(
    process.diamondRecHitProducer
)


# output configuration
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:./DiamondReco.root"),
     outputCommands = cms.untracked.vstring(
    'keep DiamondRecHitedmDetSetVector_diamondRecHitProducer_*_*'
 )
)

process.outpath = cms.EndPath(process.output)
