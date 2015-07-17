import FWCore.ParameterSet.Config as cms

RPFastFullSimulation200 = cms.EDProducer("RPFastFullSimulation200",
    verbosity = cms.untracked.uint32(5),
    hepMCLabel = cms.string('generator'),
    roundToPitch = cms.bool(True),
    pitch = cms.double(0.066),            # in mm
    dscrWidth = cms.double(10E-3),        # in mm
    thetaLimit = cms.double(0.001),       # in rad
    minUVSensorsPerRP = cms.uint32(3),
    minRPsPerStation = cms.uint32(2),
    RPsAllowed = cms.vuint32(),

    savePlots = cms.bool(True),
    plotFileName = cms.string("plots.root"),
    plot_z = cms.double(210E3),           # in mm

    dumpTopRPProtons = cms.bool(False),
    dumpBottomRPProtons = cms.bool(False),
    dumpHorizontalRPProtons = cms.bool(False),
)
