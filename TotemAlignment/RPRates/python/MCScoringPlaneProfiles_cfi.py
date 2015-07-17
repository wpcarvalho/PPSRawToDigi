import FWCore.ParameterSet.Config as cms

MCScoringPlaneProfiles = cms.EDAnalyzer("MCScoringPlaneProfiles",

    verbosity = cms.untracked.uint32(0),

    # maximum theta for which a proton is considered as forward
    forwardThLimit = cms.double(0.001),
    
    # maximum xi for which a proton is considered as forward
    forwardXiLimit = cms.double(0.3),

    # position of scoring plane
    scoringPlaneAtRP = cms.uint32(124),
    
    # offsets range
    offsets_N = cms.uint32(10),
    offsets_from = cms.double(0),
    offsets_to = cms.double(9),

    # file where the results will be written
    outputFile = cms.string('')
)


