import FWCore.ParameterSet.Config as cms

T2DEVNtuplizer = cms.EDAnalyzer("T2DEVNtuplizer",
    outputFileName = cms.untracked.string('totem_ntuple.root'),
    verbosity = cms.untracked.uint32(0),

    includeDigi = cms.bool(False),
    includePatterns = cms.bool(False),

    ProductLabelSimu = cms.string('rpCCOutput'),
    ModulLabelSimu = cms.string('Raw2DigiProducer'),

    CluLabel = cms.string('T2MCl'),
    HitLabel = cms.string('T2Hits'),

    RoadLabel = cms.string('NewRoadFinder'),
    TrackLabel = cms.string('T2TrackColl3'),
    T1TrackLabel = cms.string('t1tracks2'),

    EasyGeneratorPrimaryDef= cms.bool(False),
    LoadPrimaryGeantInformation =cms.bool(False),
    HepMCProductLabel = cms.string('generator'),

 #   vtxpos= cms.double('vtxpos'),    
#    primaryProtons = cms.bool(True),
#    primaryJets = cms.bool(True),
#    primaryJetsInstance = cms.string("CentralMCJetReco"),
#    primaryJetsLabel = cms.string("kt6algorithm"),

     SimTrackContainerLabel = cms.InputTag("g4SimHits"),
     SimVertexContainerLabel = cms.InputTag("g4SimHits"),
     T2PadDigiCollectionLabel = cms.InputTag("T2Digis", "T2PadDigi"),
     T2StripDigiCollectionLabel = cms.InputTag("T2Digis", "T2StripDigi")
)
