import FWCore.ParameterSet.Config as cms
generator = cms.EDFilter("Pythia8GeneratorFilter",
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(14000.0),
    PythiaParameters = cms.PSet(
        pythiaMinBias = cms.vstring('SoftQCD:singleDiffractive = on',
        	'PhaseSpace:pTHatMin = 15.'
        ),
        parameterSets = cms.vstring('pythiaMinBias')	
    )
)




