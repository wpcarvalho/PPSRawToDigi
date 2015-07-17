import FWCore.ParameterSet.Config as cms

from GeneratorInterface.Pythia6Interface.pythiaDefault_cff import *
generator = cms.EDFilter("PomwigGeneratorFilter",
    doMPInteraction = cms.bool(False),
    useJimmy = cms.bool(False),
    herwigHepMCVerbosity = cms.untracked.bool(False),
    filterEfficiency = cms.untracked.double(1.0),
    herwigVerbosity = cms.untracked.int32(1),
    comEnergy = cms.double(8000.0),
    printCards = cms.untracked.bool(True),
    crossSection = cms.untracked.double(-1.0),
    maxEventsToPrint = cms.untracked.int32(2),
    survivalProbability = cms.double(0.05),
    HerwigParameters = cms.PSet(
        DPEInclusiveJets = cms.vstring('NSTRU      = 14         ! H1 Pomeron Fit B', 
            'Q2WWMN     = 1E-6       ! Minimum |t|', 
            'Q2WWMX     = 4.0        ! Maximum |t|', 
            'YWWMIN     = 1E-6       ! Minimum xi', 
            'YWWMAX     = 0.2        ! Maximum xi', 
            'IPROC      = 11500      ! Process PomPom -> jets', 
            'PTMIN      = 30         ! 2->2 PT min', 
            'MODPDF(1)  = -1         ! Set MODPDF', 
            'MODPDF(2)  = -1         ! Set MODPDF'),
        parameterSets = cms.vstring('DPEInclusiveJets')
    ),
    h1fit = cms.int32(2),
    doPDGConvert = cms.bool(False),
    diffTopology = cms.int32(0)
)
