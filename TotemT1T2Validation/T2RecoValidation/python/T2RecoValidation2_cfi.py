import FWCore.ParameterSet.Config as cms

T2RecoAn2 = cms.EDAnalyzer("T2RecoAnalyzer",
    FitdzEnergy = cms.double(50.0),         # choose 10. or 50. GeV
    Chicut = cms.double(3.0),
    DZScale = cms.double(1.0),
    TrkEtamin = cms.double(5.1),            # 4.8 for 10 GeV
    TrkEtaMAX = cms.double(6.9),            # 7.0 for 10 GeV

    singleparticle = cms.bool(False),
    PhiChiProbCut = cms.double(0.01),
    RChiProbCut = cms.double(0.01),

    SimTrackContainerLabel = cms.InputTag("g4SimHits"),
    HepMCProductLabel = cms.string('generator'),
    CluLabel = cms.string('T2MCl'),
    HitLabel = cms.string('T2Hits'),
    RoadLabel = cms.string('T2RoadPadFinder'),
    TrackLabel = cms.string('T2TrackColl3'),

    OutputFile = cms.untracked.string('file:valPythia90T2PlotsRecoCC2.root')
)


