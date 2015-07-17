import FWCore.ParameterSet.Config as cms

RPFastStationSimulation = cms.EDProducer("RPFastStationSimulation",
    verbosity = cms.untracked.uint32(0),

    makeHepMC = cms.bool(False),
    makeHits = cms.bool(True),
    
    particlesPerEvent = cms.uint32(1),

    RPs = cms.vuint32(120, 121, 122, 123, 124, 125),
    z0 = cms.double(214500),

    position_distribution = cms.PSet(
      type = cms.string("box"),
      x_mean = cms.double(0.0),       #in mm
      x_width = cms.double(80.0),
      x_min = cms.double(0.0),
      x_max = cms.double(0.0),

      y_mean = cms.double(0.0),
      y_width = cms.double(80.0),
      y_min = cms.double(0.0),
      y_max = cms.double(0.0)
    ),

    angular_distribution = cms.PSet(
      type = cms.string("gauss"),
      x_mean = cms.double(0.0),       #in rad
      x_width = cms.double(10E-6),
      x_min = cms.double(0E-6),
      x_max = cms.double(0E-6),

      y_mean = cms.double(0.0),
      y_width = cms.double(10E-6),
      y_min = cms.double(0E-6),
      y_max = cms.double(0E-6)
    ),

    pitch = cms.double(66E-3),  # mm
    roundToPitch = cms.bool(True),

    dscrWidth = cms.double(10E-3),  # mm
    dscReduceUncertainty = cms.bool(True),

    insensitiveMargin = cms.double(34E-3), #mm, save value as RPActiveEdgePosition in SimTotem/RPDigiProducer/python/RPSiDetConf_cfi.py

    stationId = cms.uint32(0),
    minUVSensorsPerRP = cms.uint32(3),
    minRPsPerStation = cms.uint32(3)
)
