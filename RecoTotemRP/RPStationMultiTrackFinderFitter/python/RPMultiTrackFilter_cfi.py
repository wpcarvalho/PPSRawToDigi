import FWCore.ParameterSet.Config as cms

RPMultiTrackFilter = cms.EDFilter('RPMultiTrackFilter',
        generatorLabel = cms.InputTag(''),
        oneRPRecoLabel = cms.InputTag(''),
        stationRecoLabel = cms.InputTag(''),
        trackCount = cms.PSet(
            active = cms.bool(False),
            min = cms.uint32(1),
            max = cms.uint32(2)
        ),
        reconstructableTrackCount = cms.PSet(
            active = cms.bool(False),
            min = cms.uint32(1),
            max = cms.uint32(2)
        ),
        recoTrackCount = cms.PSet(
            active = cms.bool(False),
            min = cms.uint32(1),
            max = cms.uint32(2)
        ),
        minDistBetweenPatterns = cms.PSet(
            active = cms.bool(False),
            min = cms.double(0),
            max = cms.double(1e100)
        ),
        minDistBetweenTracks = cms.PSet(
            active = cms.bool(False),
            min = cms.double(0),
            max = cms.double(1e100)
        ),
)
