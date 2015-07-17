import FWCore.ParameterSet.Config as cms

ElasticRecoVal = cms.EDAnalyzer("ElasticRecoVal",
    verbosity = cms.untracked.uint32(2),

    # HepMC product labels
    generatorLabel = cms.string('generator'),
    originalLabel = cms.string('original'),

    # output file names (set empty for no output)
    validationOutputFile = cms.string('elastic_reconstruction_validation.root'),
    acceptanceOutputFile = cms.string('elastic_acceptances.root'),

    # default limits (plot range) values
    lim_vtx_dx = cms.untracked.double(1E-3),
    lim_vtx_dy = cms.untracked.double(1E-3),
    lim_th_dx = cms.untracked.double(10E-6),
    lim_th_dy = cms.untracked.double(10E-6),
    lim_t_min = cms.untracked.double(0),
    lim_t_max = cms.untracked.double(1),
    
    RPDetTriggerSetLabel = cms.InputTag("RPSiDetDigitizer"),
    RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit"),
    RPRecoElasticEventLabel = cms.InputTag("ElasticReconstruction")
)
