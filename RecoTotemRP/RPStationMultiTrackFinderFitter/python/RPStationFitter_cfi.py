import FWCore.ParameterSet.Config as cms

RPStationFitter = cms.EDProducer('RPStationFitter',
        verbosity = cms.untracked.uint32(0),
        RPTrackCandCollProducer = cms.string(''),

        z_h    = cms.double(212698.0),  # mm

        outside_edge_th = cms.double(0.1), # mm

        enableFitting = cms.bool(True)
)
