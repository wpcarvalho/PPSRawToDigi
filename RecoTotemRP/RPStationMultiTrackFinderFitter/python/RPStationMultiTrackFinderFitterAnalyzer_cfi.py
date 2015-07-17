import FWCore.ParameterSet.Config as cms

RPStationMultiTrackFinderFitterAnalyzer = cms.EDAnalyzer('RPStationMultiTrackFinderFitterAnalyzer',
        verbosity = cms.untracked.uint32(0),

        simTracksProducer  = cms.string(''),
        recoTracksProducer = cms.string(''),
        oneRPPatternsProducer = cms.string(''),

        outputFileName = cms.string('stats'),
        rootFileName   = cms.string('results.root'),
        histName       = cms.string('missing_fake'),
        missing_max    = cms.uint32(7),
        fake_max       = cms.uint32(7),

        si_de_ax   = cms.double(3.), # urad
        si_de_ay   = cms.double(3.), # urad
        si_de_bx   = cms.double(0.02), # mm
        si_de_by   = cms.double(0.02), # mm
        acceptance = cms.double(9.),

        outside_edge_th = cms.double(0.1) # mm, DUMMY
)

