import FWCore.ParameterSet.Config as cms

LinearTransportTest = cms.EDAnalyzer("LinearTransportTest",
    verbosity = cms.untracked.uint32(0),

    samples = cms.uint32(10),
    t_min = cms.double(1.8),
    t_max = cms.double(5.0),
    vertex_n_sigma = cms.double(1.0),
    xi_min = cms.double(-1e-05),
    xi_max = cms.double(1e-05),

    RPs = cms.vuint32(100, 104, 120, 124),

    outputFile = cms.string('linearApproximationTest.root')
)
