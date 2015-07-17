import FWCore.ParameterSet.Config as cms
import sys

if len(sys.argv) < 5:
    print "Usage: cmsRun .. <module> <tracks per event> <ev_offset> [<RPId to use> ...]"
    exit(1)

module_name = sys.argv[2]
tracks_per_ev = int(sys.argv[3])
event  = int(sys.argv[4])
RPs = [int(RPId) for RPId in sys.argv[5:]]

module = __import__(module_name, globals(), locals(), ['process'])
process = module.process

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(event)
)

# multi track generator
try:
    process.generator.verbosity = 10
    process.generator.tracksPerEvent = tracks_per_ev
    generatorOK = True
except AttributeError:
    generatorOK = False

# full simulation
try:
    process.RPFastFullSimulation200.verbosity = 10
    process.RPFastFullSimulation200.savePlots = False
    process.RPFastFullSimulation200.RPsAllowed = RPs
    sim200OK = True
except AttributeError:
    sim200OK = False

# station simulation
try:
    process.RPFastStationSimulation.particlesPerEvent = tracks_per_ev
    process.RPFastStationSimulation.RPs = RPs
    simOK = True
except AttributeError:
    simOK = False

#filter events
process.load("TotemRawData.Filters.EventNumberTimeFilter_cfi")
process.EventNumberTimeFilter.RawEventLabel = cms.InputTag("source") # still the source must be modified to ignore the check
process.EventNumberTimeFilter.eventNumber.active = True
process.EventNumberTimeFilter.eventNumber.min = event
process.EventNumberTimeFilter.eventNumber.max = event

# 1-RP non-parallel pattern recognition
process.NonParallelTrackFinder.verbosity = 10

# Hough reconstruction
process.RPStationMultiTrackFinderFitter.verbosity = 10

# Analyzer
process.RPStationMultiTrackFinderFitterAnalyzer.verbosity = 10
process.RPStationMultiTrackFinderFitterAnalyzer.rootFileName = ""
process.RPStationMultiTrackFinderFitterAnalyzer.outputFileName = ""

process.p = cms.Path(process.GeometryInfo)

if generatorOK:
    process.p *= process.generator

process.p *= process.EventNumberTimeFilter

if sim200OK:
    process.p *= process.RPFastFullSimulation200

if simOK:
    process.p *= process.RPFastStationSimulation

process.p *= process.NonParallelTrackFinder
process.p *= process.RPStationMultiTrackFinderFitter
process.p *= process.RPStationMultiTrackFinderFitterAnalyzer
