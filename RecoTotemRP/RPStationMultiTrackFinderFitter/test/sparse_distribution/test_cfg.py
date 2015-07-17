import FWCore.ParameterSet.Config as cms

from default_cfg import process

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.RPFastStationSimulation.verbosity = 0
process.NonParallelTrackFinder.verbosity = 10
process.RPStationMultiTrackFinderFitter.verbosity = 10
process.RPStationMultiTrackFinderFitterAnalyzer.verbosity = 6
