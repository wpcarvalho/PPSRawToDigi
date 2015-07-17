import FWCore.ParameterSet.Config as cms

process = cms.Process("TestBeamDataAnalysis")

process.load("Configuration.TotemCommon.LoggerMax_cfi")

process.load("L1TriggerTotem.CoincidenceChip.testbeam_data_pot5_cff")

# file with data simulation 
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:input2.root'),
#    skipEvents = cms.untracked.uint32(100)
#)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# takes raw data and calculates trigger bits
# which may serve as input to CC
process.TriggerBits = cms.EDProducer("RPTriggerBitsProducer",
    verbose = cms.bool(False)
)

# Geometry - we assume beta*=90m for analysis of good tracks
process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90/RP_Dist_Beam_Cent.xml')

process.load("RecoTotemRP.RPClusterizer.RPClusterizationConf_cfi")
# process.RPClustProd.Verbosity = 1

process.load("RecoTotemRP.RPRecoHitProducer.RPRecoHitProdConf_cfi")
# process.RPHecoHitProd.Verbosity = 1

process.load("RecoTotemRP.RPSingleCandidateTrackFinder.RPSingleTrackCandFindConf_cfi")
# process.RPSinglTrackCandFind.Verbosity = 1

process.load("RecoTotemRP.RPTrackCandidateCollectionFitter.RPSingleTrackCandCollFitted_cfi")
# process.RPSingleTrackCandCollFit.Verbosity = 1

# produce simulated response from CC, RPCC module
process.load("L1TriggerTotem.CoincidenceChip.RPCoincidenceProducer_cfi")
#process.RPCC.coincidenceChipConfig.useControlRegisters = False/True
#process.RPCC.coincidenceChipConfig.controlRegister1    = 19
#process.RPCC.coincidenceChipConfig.controlRegister2    = 242
#process.RPCC.coincidenceChipConfig.controlRegister3    = 200

# compares output from real CC with simulated CC
process.TestBeamDataAnalyzer = cms.EDAnalyzer("RPCoincidenceAnalyzer",
    verbose          = cms.uint32(3),
    outputFile       = cms.string('hists1.root'),
    modulLabelRaw    = cms.string('RawToCC'),
    productLabelRaw  = cms.string(''),
    modulLabelSimu   = cms.string('RPCC'),
    productLabelSimu = cms.string('')
)

process.TestBeamDataAnalyzer.productLabelSimu = process.RPCC.productLabelSimu
process.TestBeamDataAnalyzer.productLabelRaw  = process.RawToCC.productLabelRaw


process.o1 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:output2.root')
)

process.p1 = cms.Path(process.RawToDigi*process.RawToCC*process.TriggerBits*process.RPClustProd*process.RPHecoHitProd*process.RPSinglTrackCandFind*process.RPSingleTrackCandCollFit*process.RPCC*process.TestBeamDataAnalyzer)
#process.p1 = cms.Path(process.RawToDigi*process.RawToCC*process.TriggerBits*process.RPClustProd*process.RPHecoHitProd*process.RPSinglTrackCandFind*process.RPSingleTrackCandCollFit*process.RPCC)

process.outpath = cms.EndPath(process.o1)
