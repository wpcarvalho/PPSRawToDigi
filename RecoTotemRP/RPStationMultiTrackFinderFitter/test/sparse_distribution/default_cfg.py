import FWCore.ParameterSet.Config as cms

import os

working_dir = os.path.dirname(os.path.realpath(__file__))

process = cms.Process("RPMultiReco")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# random seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# set number of events
process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# station simulation
process.load("TotemAlignment.RPFastSimulation.RPFastStationSimulation_cfi")
process.RPFastStationSimulation.verbosity = 0
process.RPFastStationSimulation.particlesPerEvent = 2
process.RPFastStationSimulation.roundToPitch = True

process.RPFastStationSimulation.RPs = cms.vuint32(100, 104, 120)
process.RPFastStationSimulation.z0 = 212698 # mm

# more or less uniform distribution of tracks
process.RPFastStationSimulation.position_distribution.type = "box"
process.RPFastStationSimulation.position_distribution.x_mean = 0
process.RPFastStationSimulation.position_distribution.x_width = 17 # mm
process.RPFastStationSimulation.position_distribution.y_mean = 14.5
process.RPFastStationSimulation.position_distribution.y_width = 17 # mm

process.RPFastStationSimulation.angular_distribution.type = "gauss"
process.RPFastStationSimulation.angular_distribution.x_mean = 0
process.RPFastStationSimulation.angular_distribution.x_width = 200E-6
process.RPFastStationSimulation.angular_distribution.y_mean = 0
process.RPFastStationSimulation.angular_distribution.y_width = 200E-6

process.RPFastStationSimulation.stationId = 0
process.RPFastStationSimulation.minUVSensorsPerRP = 3
process.RPFastStationSimulation.minRPsPerStation = 3

# beam/optics parameters
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_1535_cfi")

# base geometry
process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Garage/RP_Dist_Beam_Cent.xml')

# alignment corrections
process.load("TotemAlignment.RPDataFormats.TotemRPIncludeAlignments_cfi")

# RP (not sensor) edge at the beam
alignment_files = cms.vstring(working_dir+"/../alignment_station_200m.xml")

process.TotemRPIncludeAlignments.MisalignedFiles = alignment_files
process.TotemRPIncludeAlignments.RealFiles = alignment_files

# geometry printer
process.GeometryInfo = cms.EDAnalyzer("GeometryInfoModule")

# 1-RP parallel pattern recognition
process.load("RecoTotemRP.RPSingleCandidateTrackFinder.RPSingleTrackCandFindConf_cfi")
process.RPSinglTrackCandFind.Verbosity = 0
process.RPSinglTrackCandFind.RoadSize = 0.05 # mm
process.RPSinglTrackCandFind.MinHitsPerCoord = 3
process.RPSinglTrackCandFind.MaxHitsPerDetector = 15
process.RPSinglTrackCandFind.ReduceWeightsWithMultiplicity = False
process.RPSinglTrackCandFind.RPRecoHitLabel = cms.InputTag("RPFastStationSimulation")

# 1-RP non-parallel pattern recognition
process.load("RecoTotemRP.RPNonParallelTrackCandidateFinder.RPNonParallelTrackCandidateFinder_cfi")
process.NonParallelTrackFinder.verbosity = 0
process.NonParallelTrackFinder.maxHitsPerPlaneToSearch = 5
process.NonParallelTrackFinder.minPlanesPerProjectionToSearch = 3
process.NonParallelTrackFinder.minPlanesPerProjectionToFit = 3
process.NonParallelTrackFinder.clusterSize_a = 0.002
process.NonParallelTrackFinder.clusterSize_b = 0.06
process.NonParallelTrackFinder.threshold = 2.99
process.NonParallelTrackFinder.allowAmbiguousCombination = True
process.NonParallelTrackFinder.DetSetVectorRPRecoHitLabel = cms.InputTag("RPFastStationSimulation")

# multi-track reconstruction
process.load("RecoTotemRP.RPStationMultiTrackFinderFitter.RPStationMultiTrackFinderFitter_cfi")
process.RPStationMultiTrackFinderFitter.verbosity = 0
#process.RPStationMultiTrackFinderFitter.RPTrackCandCollProducer = 'RPSinglTrackCandFind'
process.RPStationMultiTrackFinderFitter.RPTrackCandCollProducer = 'NonParallelTrackFinder'

process.RPStationMultiTrackFinderFitter.ax_cut = 800e-6
process.RPStationMultiTrackFinderFitter.ay_cut = 800e-6
process.RPStationMultiTrackFinderFitter.z_h = 212698.0
process.RPStationMultiTrackFinderFitter.cluster_merge_si_cut = 5.

# analysis
process.load("RecoTotemRP.RPStationMultiTrackFinderFitter.RPStationMultiTrackFinderFitterAnalyzer_cfi")
process.RPStationMultiTrackFinderFitterAnalyzer.verbosity = 0
process.RPStationMultiTrackFinderFitterAnalyzer.simTracksProducer = cms.string("RPFastStationSimulation")
process.RPStationMultiTrackFinderFitterAnalyzer.recoTracksProducer = cms.string("RPStationMultiTrackFinderFitter")
process.RPStationMultiTrackFinderFitterAnalyzer.oneRPPatternsProducer = cms.string("NonParallelTrackFinder")

process.RPStationMultiTrackFinderFitterAnalyzer.outputFileName = cms.string("analyzer.out")
process.RPStationMultiTrackFinderFitterAnalyzer.rootFileName = cms.string("analyzer.root")

# event filter
process.load('TotemRawData.Filters.EventNumberTimeFilter_cfi')
process.EventNumberTimeFilter.eventNumber.active = True
event_number = 4
process.EventNumberTimeFilter.eventNumber.min = event_number
process.EventNumberTimeFilter.eventNumber.max = event_number

#process.dump = cms.EDAnalyzer('EventContentAnalyzer')

process.p = cms.Path(
    process.GeometryInfo *
    process.RPFastStationSimulation *
    process.EventNumberTimeFilter *
    #process.RPSinglTrackCandFind *
    process.NonParallelTrackFinder *
    process.RPStationMultiTrackFinderFitter *
    process.RPStationMultiTrackFinderFitterAnalyzer
)
