import FWCore.ParameterSet.Config as cms

t1tracksanalTB = cms.EDAnalyzer("T1TrackAnalyzerTB",
    Verbosity = cms.int32(0),
    SeeTrack = cms.double(1.0),
    SeeHits = cms.double(1.0),
    ChiOverNCut = cms.double(20.0),
    ZRange = cms.double(1000.0),
    Zmin = cms.double(-10000.0),
    Zmax = cms.double(10000.0),
    realBeamPosX = cms.double(0.0),
    realBeamPosY = cms.double(0.0),
    realBeamAngle = cms.double(0.0)

)


