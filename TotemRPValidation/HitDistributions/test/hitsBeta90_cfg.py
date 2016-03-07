import FWCore.ParameterSet.Config as cms

process = cms.Process("RPinelasticBeta90Energy7TeV")

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Configure the output module (save the result in a file -- RPinelastic90.root)
process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *', 
        'drop *_*_TrackerHits*_*', 
        'drop *_*_Muon*_*', 
        'drop *_*_Ecal*_*', 
        'drop *_*_Hcal*_*', 
        'drop *_*_Calo*_*', 
        'drop *_*_Castor*_*', 
        'drop *_*_FP420SI_*', 
        'drop *_*_ZDCHITS_*'),
    fileName = cms.untracked.string('file:RPinelastic90.root')
)

process.load("Configuration.TotemCommon.LoggerMin_cfi")

# Use random number generator service
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# Use particle table
process.load("SimGeneral.HepPDTESSource.pdt_cfi")

# Geometry - beta* specific
process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90/RP_Dist_Beam_Cent.xml')

# declare optics parameters
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_90_cfi")

# Magnetic Field
# by default we have 3.8T
process.load("Configuration.StandardSequences.MagneticField_cff")

process.source = cms.Source("EmptySource")

# Monte Carlo gun 
process.load("IOMC.FlatProtonLogKsiLogTGun.Beta90_cfi")
# process.FlatProtonLogKsiLogTGun.Verbosity = 1

# Smearing
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")
process.SmearingGenerator.originalLabel = 'original'
process.SmearingGenerator.modifyLabel = 'generator'
process.SmearingGenerator.verbosity = 0

# Oscar - G4 simulation & proton transport
process.load("Configuration.TotemCommon.g4SimHits_cfi")
process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup
# The input source for G4 module is connected to "process.generator".
process.g4SimHits.Generator.HepMCProductLabel = 'generator'    # energy+vertex smearing
process.g4SimHits.G4EventManagerVerbosity = 0
process.g4SimHits.G4StackManagerVerbosity = 0
process.g4SimHits.G4TrackingManagerVerbosity = 0
process.g4SimHits.MagneticField.Verbosity = False
process.g4SimHits.Physics.Verbosity = 0
process.g4SimHits.Physics.BeamProtTransportSetup.Verbosity = False
process.g4SimHits.Generator.Verbosity = 0
process.g4SimHits.SteppingAction.Verbosity = 0
process.g4SimHits.Totem_RP_SD.Verbosity = 0
process.g4SimHits.TotemSD.Verbosity = 0

# No pile up for the mixing module
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

# RP Strip digitization
process.load("SimTotem.RPDigiProducer.RPSiDetConf_cfi")
# process.RPSiDetDigitizer.RPVerbosity = 1

process.load("RecoTotemRP.RPClusterizer.RPClusterizationConf_cfi")
# process.RPClustProd.Verbosity = 1

process.load("RecoTotemRP.RPRecoHitProducer.RPRecoHitProdConf_cfi")
# process.RPHecoHitProd.Verbosity = 1

process.load("RecoTotemRP.RPSingleCandidateTrackFinder.RPSingleTrackCandFindConf_cfi")
# process.RPSinglTrackCandFind.Verbosity = 1

process.load("RecoTotemRP.RPTrackCandidateCollectionFitter.RPSingleTrackCandCollFitted_cfi")
# process.RPSingleTrackCandCollFit.Verbosity = 1

process.load("RecoTotemRP.RPInelasticReconstruction.RPRec220_cfi")
process.RP220Reconst.BeamProtTransportSetup = process.BeamProtTransportSetup
# process.RP220Reconst.Verbosity = 1

process.load("TotemRPValidation.HitDistributions.HitDistributions_cfi")
process.RPHitDists.outputFile = cms.string('file:hits.root')

process.p1 = cms.Path(process.generator*process.SmearingGenerator*process.g4SimHits*process.mix*process.RPSiDetDigitizer*process.RPClustProd*process.RPHecoHitProd*process.RPSinglTrackCandFind*process.RPSingleTrackCandCollFit*process.RP220Reconst*process.RPHitDists)

process.outpath = cms.EndPath(process.o1)
