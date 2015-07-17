import FWCore.ParameterSet.Config as cms

NonParallelTrackFinder = cms.EDProducer("RPNonParallelTrackCandidateFinder",

    verbosity = cms.untracked.uint32(0),

    # (full) cluster size in slope-intercept space
    clusterSize_a = cms.double(0.02), # rad
    clusterSize_b = cms.double(0.3),  # mm

    maxHitsPerPlaneToSearch = cms.uint32(7),
    minPlanesPerProjectionToSearch = cms.uint32(3),
    minPlanesPerProjectionToFit = cms.uint32(3),

    threshold = cms.double(2.99),

    allowAmbiguousCombination = cms.bool(False),

    # maximal angle (in any projection) to mark the candidate as fittable -> controls track parallelity with beam
	# huge value -> no constraint
    max_a_toFit = cms.double(10.0),

    DetSetVectorRPRecoHitLabel = cms.InputTag("RPHecoHitProd")
)
