import FWCore.ParameterSet.Config as cms

import os

working_dir = os.path.dirname(os.path.realpath(__file__))

process = cms.Process("RPMultiReco4000GeV90m")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# random seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")
process.RandomNumberGeneratorService.generator.initialSeed = 543

# set number of events
process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# multi track generator
process.generator = cms.EDProducer("CachedMultiTrackSource",
    verbosity = cms.untracked.uint32(0),

    tracksPerEvent = cms.uint32(5),

    input = cms.VPSet(
        cms.PSet(
          weight = cms.double(25),
          fileName = cms.string(working_dir+"/sample_ES"),
        ),
        cms.PSet(
          weight = cms.double(7),
          fileName = cms.string(working_dir+"/sample_SD"),
        )
    )
)

# load common settings (optics, geometry)
process.load(".".join(os.path.relpath(working_dir).split("/")+["common_cfi"]))

# fast simulation
process.load("TotemAlignment.RPFastSimulation.RPFastFullSimulation200_cfi")
process.RPFastFullSimulation200.plotFileName = "simulation.root"

# geometry printer
process.GeometryInfo = cms.EDAnalyzer("GeometryInfoModule")

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
process.NonParallelTrackFinder.DetSetVectorRPRecoHitLabel = cms.InputTag("RPFastFullSimulation200")

# MultiTrack reconstruction
process.load("RecoTotemRP.RPStationMultiTrackFinderFitter.RPStationMultiTrackFinderFitter_cfi")
process.RPStationMultiTrackFinderFitter.verbosity = 0
process.RPStationMultiTrackFinderFitter.RPTrackCandCollProducer = 'NonParallelTrackFinder'

#process.RPStationMultiTrackFinderFitter.ax_cut = 49.53*4*1e-6
#process.RPStationMultiTrackFinderFitter.ax_mean = -8.164e-6
#process.RPStationMultiTrackFinderFitter.ay_cut = 87.78*4*1e-6
#process.RPStationMultiTrackFinderFitter.ay_mean = 261.5e-6

# analysis
process.load("RecoTotemRP.RPStationMultiTrackFinderFitter.RPStationMultiTrackFinderFitterAnalyzer_cfi")
process.RPStationMultiTrackFinderFitterAnalyzer.verbosity = 0
process.RPStationMultiTrackFinderFitterAnalyzer.simTracksProducer = cms.string("RPFastFullSimulation200")
process.RPStationMultiTrackFinderFitterAnalyzer.recoTracksProducer = cms.string("RPStationMultiTrackFinderFitter")
process.RPStationMultiTrackFinderFitterAnalyzer.oneRPPatternsProducer = cms.string("NonParallelTrackFinder")

process.p = cms.Path(
    process.GeometryInfo *
    process.generator *
    process.RPFastFullSimulation200 *
    process.NonParallelTrackFinder *
    process.RPStationMultiTrackFinderFitter *
    process.RPStationMultiTrackFinderFitterAnalyzer
)
