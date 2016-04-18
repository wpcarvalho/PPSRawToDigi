import FWCore.ParameterSet.Config as cms

TotemNtuplizer = cms.EDAnalyzer("TotemNtuplizer",
    outputFileName = cms.untracked.string('totem_ntuple.root'),
    verbosity = cms.untracked.uint32(0),

    includeDigi = cms.bool(False),
    includePatterns = cms.bool(False),

    RawEventLabel = cms.InputTag("RPSiDetDigitizer"),

    ModulLabelSimu = cms.string('Raw2DigiProducer'),
    TotemRPDigiSetLabel = cms.InputTag("RPSiDetDigitizer"),

    ProductLabelSimu = cms.string('rpCCOutput'),

    RPDigClusterLabel = cms.InputTag("RPClustProd"),
    RPFittedTrackCollectionLabel = cms.InputTag("TotemRPLocalTrackFitter"),
    RPMulFittedTrackCollectionLabel = cms.InputTag("RPMulTrackCandCollFit"),

    RPReconstructedProtonCollectionLabel = cms.InputTag("ElasticReconstruction"),
    RPReconstructedProtonPairCollectionLabel = cms.InputTag("ElasticReconstruction"),

    RoadLabel = cms.string('NewRoadFinder')

#    primaryProtons = cms.bool(True),
#    primaryJets = cms.bool(True),
#    primaryJetsInstance = cms.string("CentralMCJetReco"),
#    primaryJetsLabel = cms.string("kt6algorithm")
)
