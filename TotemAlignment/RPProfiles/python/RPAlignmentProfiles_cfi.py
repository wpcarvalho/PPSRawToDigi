import FWCore.ParameterSet.Config as cms
from math import sqrt

pitch = 66E-3/sqrt(2) # mm
n = int(25/pitch)
span = pitch * (n + 0.5)

vov_binsize = pitch
vov_n = n

RPAlignmentProfiles = cms.EDAnalyzer("RPAlignmentProfiles",
    verbosity = cms.untracked.uint32(0),

    # the numbers (= last digits) of RPs that will be processed
    rpNums = cms.vuint32(0, 1, 2, 3, 4, 5),

    # will accept LocalTrackFitter results with two units only
    requireTwoUnits = cms.bool(True),

    # the name of the track-candidate producer module
    #   use `RPSinglTrackCandFind' for parallel finder
    #   use `NonParallelTrackFinder' for non-parallel finder
    RPTrackCandidateCollectionLabel = cms.InputTag("RPSinglTrackCandFind"),

    # name of the track fit producer
    RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit"),

    outputFile = cms.string('profiles.root'),
    resultFile = cms.string('profiles.xml'),

    hit_profiles = cms.PSet(
        x_min = cms.double(-span),       # in mm
        x_max = cms.double(+span),
        y_min = cms.double(-span),       # in mm
        y_max = cms.double(+span),
        bins_1d = cms.uint32(2*n + 1),
        bins_2d = cms.uint32(50),

        x_slices_number = cms.uint32(30),
        x_slices_start = cms.double(0),  # in mm
        x_slice_width = cms.double(0.5), # in mm
        
        y_slices_number = cms.uint32(60),
        y_slices_start = cms.double(-15),  # in mm
        y_slice_width = cms.double(0.5), # in mm
    ),

    vertOfVert = cms.PSet(
        y_min = cms.double(-vov_n*vov_binsize),       # in mm
        y_max = cms.double(+vov_n*vov_binsize),
        bins_1d = cms.uint32(2*vov_n),
        
        x_slices_number = cms.uint32(40),
        x_slices_start = cms.double(-10),  # in mm
        x_slice_width = cms.double(0.5) # in mm
    ),

    angular_profiles = cms.PSet(
        x_min = cms.double(-20),       # in mrad
        x_max = cms.double(20),
        y_min = cms.double(-20),       # in mrad
        y_max = cms.double(20),
        bins_1d = cms.uint32(4000)
    ),

    thresholdToMax_main = cms.double(0.6),
    thresholdToMax_secondary = cms.double(0.8),

    maxResidualToSigma = cms.double(3),
    minimumHitsPerProjectionPerRP = cms.uint32(4)
)
