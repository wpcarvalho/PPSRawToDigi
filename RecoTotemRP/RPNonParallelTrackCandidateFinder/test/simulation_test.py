import FWCore.ParameterSet.Config as cms

process = cms.Process("RPNonParallelTrackCandidateFinderSimulationTest")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# random seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# set number of events
process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(300)
)

# beam/optics parameters
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_1535_cfi")

# station simulation
process.load("TotemAlignment.RPFastSimulation.RPFastStationSimulation_cfi")
process.RPFastStationSimulation.verbosity = 10
process.RPFastStationSimulation.particlesPerEvent = 1
process.RPFastStationSimulation.roundToPitch = True

process.RPFastStationSimulation.RPs = cms.vuint32(120, 124)
process.RPFastStationSimulation.z0 = 213000 # mm

process.RPFastStationSimulation.position_distribution.type = "box"
process.RPFastStationSimulation.position_distribution.x_mean = 0
process.RPFastStationSimulation.position_distribution.x_width = 20 # mm
process.RPFastStationSimulation.position_distribution.y_mean = 20
process.RPFastStationSimulation.position_distribution.y_width = 20 # mm

process.RPFastStationSimulation.angular_distribution.type = "gauss"
process.RPFastStationSimulation.angular_distribution.x_mean = 0
process.RPFastStationSimulation.angular_distribution.x_width = 200E-6
process.RPFastStationSimulation.angular_distribution.y_mean = 0
process.RPFastStationSimulation.angular_distribution.y_width = 200E-6

process.RPFastStationSimulation.stationId = 0
process.RPFastStationSimulation.minUVSensorsPerRP = 3
process.RPFastStationSimulation.minRPsPerStation = 3

# base geometry
process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Garage/RP_Dist_Beam_Cent.xml')

# alignment corrections
process.load("TotemAlignment.RPDataFormats.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.MisalignedFiles = cms.vstring("./test/simulation_alignment.xml")
process.TotemRPIncludeAlignments.RealFiles = cms.vstring("./test/simulation_alignment.xml")

# geometry printer
process.GeometryInfo = cms.EDAnalyzer("GeometryInfoModule")

# 1-RP non-parallel pattern recognition
process.load("RecoTotemRP.RPNonParallelTrackCandidateFinder.RPNonParallelTrackCandidateFinder_cfi")
process.NonParallelTrackFinder.DetSetVectorRPRecoHitLabel = cms.InputTag("RPFastStationSimulation")
process.NonParallelTrackFinder.verbosity = 10
process.NonParallelTrackFinder.maxHitsPerPlaneToSearch = 5
process.NonParallelTrackFinder.minPlanesPerProjectionToSearch = 3
process.NonParallelTrackFinder.minPlanesPerProjectionToFit = 3
process.NonParallelTrackFinder.threshold = 2.99

process.p = cms.Path(
   process.GeometryInfo *
    process.RPFastStationSimulation *
    process.NonParallelTrackFinder
)
