import FWCore.ParameterSet.Config as cms

totemNtuplizer = cms.EDAnalyzer("TotemNtuplizer",
    verbosity = cms.untracked.uint32(0),

    outputFileName = cms.untracked.string('totem_ntuple.root'),

    # available modules: raw, trigger, rp
    modules = cms.vstring("raw", "trigger", "rp"),

    # raw metadata part
    FEDInfosLabel = cms.InputTag("totemRPRawToDigi"),

    # trigger-counter part
    TriggerCountersLabel = cms.InputTag("totemTriggerRawToDigi"),

    # RP part
    includeDigi = cms.bool(False),
    includePatterns = cms.bool(False),

    ModulLabelSimu = cms.string('Raw2DigiProducer'),

    TotemRPDigiSetLabel = cms.InputTag("RPSiDetDigitizer"),
    ProductLabelSimu = cms.string('rpCCOutput'),
    RPDigClusterLabel = cms.InputTag("totemRPClusterProducer"),
    RPFittedTrackCollectionLabel = cms.InputTag("totemRPLocalTrackFitter"),
    RPMulFittedTrackCollectionLabel = cms.InputTag("RPMulTrackCandCollFit"),
    RPUVPatternLabel = cms.InputTag("totemRPUVPatternFinder"),

    RPReconstructedProtonCollectionLabel = cms.InputTag("ElasticReconstruction"),
    RPReconstructedProtonPairCollectionLabel = cms.InputTag("ElasticReconstruction"),

#    primaryProtons = cms.bool(True),
#    primaryJets = cms.bool(True),
#    primaryJetsInstance = cms.string("CentralMCJetReco"),
#    primaryJetsLabel = cms.string("kt6algorithm")
)
