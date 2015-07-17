import FWCore.ParameterSet.Config as cms

RPFastFullSimulation = cms.EDProducer("RPFastFullSimulation",
    verbosity = cms.untracked.uint32(0),
    hepMCLabel = cms.string('generator'),
    roundToPitch = cms.bool(True),
    pitch = cms.double(0.066),            # in mm
    dscrWidth = cms.double(10E-3),        # in mm
    thetaLimit = cms.double(0.001),       # in rad

    insensitiveMargin = cms.double(34E-3), #mm, save value as RPActiveEdgePosition in SimTotem/RPDigiProducer/python/RPSiDetConf_cfi.py

    minUVSensorsPerRP = cms.uint32(3),


#	beam_misalignment = cms.VPSet(
#	    cms.PSet(
#	        station = cms.uint32(2),
#	        z0 = cms.double(-214000),     # mm
#	        shift_x = cms.double(0),
#	        shift_y = cms.double(-100),   # um
#	        tilt_x = cms.double(0),
#	        tilt_y = cms.double(1e-2) # mrad
#	    )
#	)

    beam_misalignment = cms.VPSet()
)
