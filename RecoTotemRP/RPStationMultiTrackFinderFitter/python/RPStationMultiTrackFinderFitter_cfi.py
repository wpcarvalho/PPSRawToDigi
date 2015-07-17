import FWCore.ParameterSet.Config as cms

RPStationMultiTrackFinderFitter = cms.EDProducer('RPStationMultiTrackFinderFitter',
        verbosity = cms.untracked.uint32(0),
        RPTrackCandCollProducer = cms.string(''),

        ax_cut = cms.double(800E-6), # rad
        ax_mean = cms.double(0), # rad
        ay_cut = cms.double(800E-6), # rad
        ay_mean = cms.double(0), # rad

        z_h    = cms.double(212698.0),  # mm

        edge_create_th = cms.double(7.),
        uv_unc = cms.double(20e-3),

        filter_dist_cut = cms.double(4.),
        track_pattern_assoc_cut = cms.double(4.),
        outside_edge_th = cms.double(0.1), # mm

        cluster_merge_si_cut = cms.double(5.),

        minRPsPerStation = cms.uint32(3),

        enableFitting = cms.bool(True)
)
