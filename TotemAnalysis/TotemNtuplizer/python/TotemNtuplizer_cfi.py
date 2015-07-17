import FWCore.ParameterSet.Config as cms

TotemNtuplizer = cms.EDAnalyzer("TotemNtuplizer",
    outputFileName = cms.untracked.string('totem_ntuple.root'),
    verbosity = cms.untracked.uint32(0),

    includeDigi = cms.bool(False),
    includePatterns = cms.bool(False),

    ModulLabelSimu = cms.string('Raw2DigiProducer'),

    ProductLabelSimu = cms.string('rpCCOutput'),
    RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit"),
    RPMulFittedTrackCollectionLabel = cms.InputTag("RPMulTrackCandCollFit"),
    RPStripDigiSetLabel = cms.InputTag("RPSiDetDigitizer"),
    RawEventLabel = cms.InputTag("RPSiDetDigitizer"),
    RPDigClusterLabel = cms.InputTag("RPClustProd"),
    RPReconstructedProtonCollectionLabel = cms.InputTag("ElasticReconstruction"),
    RPReconstructedProtonPairCollectionLabel = cms.InputTag("ElasticReconstruction"),

    T2PadClusterCollectionLabel = cms.InputTag("T2MCl", "T2PadClusters"),
    T2StripClusterCollectionLabel = cms.InputTag("T2MCl", "T2StripClusters"),
    HitLabel = cms.InputTag("T2Hits", "T2Hits"),
    T2PadDigiCollectionLabel = cms.InputTag("T2Digis", "T2PadDigi"),
    T2StripDigiCollectionLabel = cms.InputTag("T2Digis", "T2StripDigi"),
    TrackLabel = cms.InputTag("T2TrackColl3", "T2TrackColl"),

    RoadLabel = cms.string('NewRoadFinder'),
    T1TrackLabel = cms.string('t1tracks2'),
    T1DigiWireCollectionLabel = cms.InputTag("T1Digis", "T1DigiWire"),
    T1DigiVfatCollectionLabel= cms.InputTag("T1Digis", "T1DigiVfat"),
    T1RecHit2DCollectionLabel = cms.InputTag("t1rechit")
    

#    primaryProtons = cms.bool(True),
#    primaryJets = cms.bool(True),
#    primaryJetsInstance = cms.string("CentralMCJetReco"),
#    primaryJetsLabel = cms.string("kt6algorithm")
)
