import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("FlatThetaXYGun",
    verbosity = cms.untracked.uint32(1),
    
    th_x_min = cms.double(-300E-6),
    th_x_max = cms.double(+300E-6),
    th_y_min = cms.double(+20E-6),
    th_y_max = cms.double(+110E-6),
    
    energy = cms.double(4000),      # GeV
    
    particle_id = cms.int32(2212),  # proton by default

    generate_right_arm = cms.bool(True),
    generate_left_arm = cms.bool(True),
)
