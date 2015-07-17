import FWCore.ParameterSet.Config as cms

emptySource = $emptySource

process = cms.Process("TestBeamAlignment")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# select test beam data for a RP
if emptySource:
    process.source = cms.Source("EmptySource")
    process.maxEvents = cms.untracked.PSet(
      input = cms.untracked.int32(0)
    )
else:
    import sys
    sys.path.append('/data2/totem/RPtest/python')
    process.load('$testbeamData')
    process.source.verbosity = 10
    process.RawToDigi.verbosity = 10

# clusterization
process.load("RecoTotemRP.RPClusterSigmaService.ClusterSigmaServiceConf_cfi")
process.load("RecoTotemRP.RPClusterizer.RPClusterizationConf_cfi")

# reco hit production
process.load("RecoTotemRP.RPRecoHitProducer.RPRecoHitProdConf_cfi")

process.TotemRPIncludeAlignments = cms.ESProducer("TotemRPIncludeAlignments"
    $misalignedString
    $realString
)

# TB geometry
process.load("Configuration.TotemCommon.geometryH8_cfi")
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

#process.GI = cms.EDAnalyzer("GeometryInfoModule")

# track search
process.load("RecoTotemRP.RPNonParallelTrackCandidateFinder.RPNonParallelTrackCandidateFinder_cfi")
process.NonParallelTrackFinder.minPlanesPerProjectionToSearch = 3
process.NonParallelTrackFinder.minPlanesPerProjectionToFit = 4

# track fitting
process.load("RecoTotemRP.RPTrackCandidateCollectionFitter.RPSingleTrackCandCollFitted_cfi")

# alignment
process.load("TotemAlignment.RPTrackBased.RPStraightTrackAligner_cfi")
process.RPStraightTrackAligner.verbosity = 5
process.RPStraightTrackAligner.maxEvents = 0

process.RPStraightTrackAligner.resolveShR = $resolveShR
process.RPStraightTrackAligner.resolveRotZ = $resolveRotZ
process.RPStraightTrackAligner.resolveShZ = False
process.RPStraightTrackAligner.resolveRPShZ = False
process.RPStraightTrackAligner.RPIds = [120]
process.RPStraightTrackAligner.z0 = 0

process.RPStraightTrackAligner.algorithms = cms.vstring($algorithms)

process.RPStraightTrackAligner.constraintsType = "$constraintsType"
process.RPStraightTrackAligner.fixedDetectorsConstraints.ShR.ids = cms.vuint32(1200, 1201, 1208, 1209)
process.RPStraightTrackAligner.fixedDetectorsConstraints.RotZ.ids = cms.vuint32(1200, 1201)

process.RPStraightTrackAligner.useExtendedRotZConstraint = $useExtendedConstraints
process.RPStraightTrackAligner.useExtendedShZConstraints = $useExtendedConstraints
process.RPStraightTrackAligner.useExtendedRPShZConstraint = $useExtendedConstraints

process.RPStraightTrackAligner.singularLimit = 1E-6
process.RPStraightTrackAligner.chiSqPerNdfCut = 100
process.RPStraightTrackAligner.minimumHitsPerProjectionPerRP = 3
process.RPStraightTrackAligner.useExternalFitter = False

process.RPStraightTrackAligner.selectionStatisticsFile = cms.string('statistics.root')
process.RPStraightTrackAligner.fileNamePrefix = cms.string('results_')
process.RPStraightTrackAligner.cumulativeFileNamePrefix = cms.string('cumulative_')
process.RPStraightTrackAligner.JanAlignmentAlgorithm.residuaHistogramsFile = 'residuals.root'

if emptySource:
  process.p = cms.Path(process.RPClustProd*process.RPHecoHitProd*process.NonParallelTrackFinder*process.RPStraightTrackAligner)
else:
  process.p = cms.Path(process.RawToDigi*process.RPClustProd*process.RPHecoHitProd*process.NonParallelTrackFinder*process.RPStraightTrackAligner)
