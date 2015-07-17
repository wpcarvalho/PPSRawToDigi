import FWCore.ParameterSet.Config as cms

ScoringPlaneProfiles = cms.EDAnalyzer("ScoringPlaneProfiles",
    verbosity = cms.untracked.uint32(0),
    
    # offsets range
    offsets_N = cms.uint32(10),
    offsets_from = cms.double(0),
    offsets_to = cms.double(9),
    
    # nothing saved when empty
    outputFile = cms.string('')	
)


