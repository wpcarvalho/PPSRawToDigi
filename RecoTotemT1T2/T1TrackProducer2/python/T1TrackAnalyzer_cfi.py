import FWCore.ParameterSet.Config as cms

t1tracksanal = cms.EDAnalyzer("T1TrackAnalyzer",
    Verbosity = cms.int32(0),
    SeeTrack = cms.double(1.0),
    SeeHits = cms.double(0.0),
    ChiOverNCut = cms.double(20.0),
    ZRange = cms.double(1000.0),
    SimVertexContainerLabel = cms.InputTag("g4SimHits"),
    SimTrackContainerLabel = cms.InpuTag("g4Simhits")
)


