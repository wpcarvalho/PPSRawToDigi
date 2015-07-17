import FWCore.ParameterSet.Config as cms

TransportValidation = cms.EDAnalyzer("TransportValidation",
    verbosity = cms.untracked.uint32(1),

    # maximum theta for which a proton is considered as forward
    forwardThLimit = cms.double(0.001),

    # file where the results will be written
    outputFile = cms.string(''),
    
    RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit")
)
