import FWCore.ParameterSet.Config as cms

process = cms.Process("RPBeta90Energy7TeV")

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
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
    fileName = cms.untracked.string('file:input2.root')
)

# Configure if you want to detail or simple log information.
# LoggerMax -- detail log info output including:
#   - errors.log
#   - warnings.log
#   - infos.log
#   - debugs.log
# LoggerMin -- simple log info output to the standard output (e.g. screen)
#process.load("Configuration.TotemCommon.LoggerMax_cfi")
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

# Monte Carlo gun 
process.load("IOMC.FlatProtonLogKsiLogTGun.Beta90_cfi")
# process.FlatProtonLogKsiLogTGun.Verbosity = 1

# Smearing
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")
process.SmearingGenerator.originalLabel = 'original'
process.SmearingGenerator.modifyLabel = 'source'
process.SmearingGenerator.verbosity = 0

# Oscar - G4 simulation & proton transport
process.load("Configuration.TotemCommon.g4SimHits_cfi")
process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup
process.g4SimHits.Generator.HepMCProductLabel = 'source'    # energy+vertex smearing
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

process.p1 = cms.Path(process.SmearingGenerator*process.g4SimHits*process.mix*process.RPSiDetDigitizer*process.RPClustProd*process.RPHecoHitProd*process.RPSinglTrackCandFind*process.RPSingleTrackCandCollFit)

process.outpath = cms.EndPath(process.o1)
