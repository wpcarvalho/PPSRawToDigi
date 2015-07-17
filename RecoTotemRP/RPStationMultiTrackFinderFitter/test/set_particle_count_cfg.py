import FWCore.ParameterSet.Config as cms
import sys

if len(sys.argv) < 5:
    print "Usage: cmsRun .. <module> <tracks per event> <num of events> [<RPId to use> ...]"
    exit(1)

module_name = sys.argv[2]
tracks_per_ev = int(sys.argv[3])
ev_count = int(sys.argv[4])
RPs = [int(RPId) for RPId in sys.argv[5:]]

module = __import__(module_name, globals(), locals(), ['process'])
process = module.process

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(ev_count)
)

# multi track generator
try:
    process.generator.tracksPerEvent = tracks_per_ev
    generatorOK = True
except AttributeError:
    generatorOK = False

# full simulation
try:
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

#process.load("RecoTotemRP.RPStationMultiTrackFinderFitter.RPStationFitter_cfi")
#process.RPStationFitter.RPTrackCandCollProducer = "NonParallelTrackFinder"

# filtering
process.load("RecoTotemRP.RPStationMultiTrackFinderFitter.RPMultiTrackFilter_cfi")
if sim200OK:
    process.RPMultiTrackFilter.generatorLabel = cms.InputTag('RPFastFullSimulation200')
elif simOK:
    process.RPMultiTrackFilter.generatorLabel = cms.InputTag('RPFastStationSimulation')
process.RPMultiTrackFilter.oneRPRecoLabel = cms.InputTag('NonParallelTrackFinder')
#process.RPMultiTrackFilter.stationRecoLabel = cms.InputTag('RPStationFitter')
process.RPMultiTrackFilter.stationRecoLabel = cms.InputTag('RPStationMultiTrackFinderFitter')
process.RPMultiTrackFilter.reconstructableTrackCount = cms.PSet(
    active = cms.bool(True),
    min = cms.uint32(tracks_per_ev),
    max = cms.uint32(tracks_per_ev+1)
)
#process.RPMultiTrackFilter.minDistBetweenTracks = cms.PSet(
    #active = cms.bool(True),
    #min = cms.double(0),
    #max = cms.double(1000)
#)

# analysis
process.RPStationMultiTrackFinderFitterAnalyzer.verbosity = 10
process.RPStationMultiTrackFinderFitterAnalyzer.rootFileName = "hists_"+str(tracks_per_ev)+".root"
process.RPStationMultiTrackFinderFitterAnalyzer.outputFileName = "stats_"+str(tracks_per_ev)
#process.RPStationMultiTrackFinderFitterAnalyzer.recoTracksProducer = cms.string("RPStationFitter")

process.p = cms.Path(process.GeometryInfo)

if generatorOK:
    process.p *= process.generator
if sim200OK:
    process.p *= process.RPFastFullSimulation200
if simOK:
    process.p *= process.RPFastStationSimulation

process.p *= process.NonParallelTrackFinder
process.p *= process.RPStationMultiTrackFinderFitter
#process.p *= process.RPStationFitter
process.p *= process.RPMultiTrackFilter
process.p *= process.RPStationMultiTrackFinderFitterAnalyzer
