import FWCore.ParameterSet.Config as cms

process = cms.Process("RPinelasticBeta2.5Energy3.5TeV")

################################################
#production

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Configure the output module (save the result in a file -- RPinelasticBeta2.5Energy3.5TeV.root)
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
    fileName = cms.untracked.string('file:RPinelasticBeta2.5Energy3.5TeV.root')
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
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_3500_Beta_2.5_220/RP_Dist_Beam_Cent.xml')

# declare optics parameters
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_3500GeV_2p5_cfi")

# Magnetic Field
# by default we have 3.8T
process.load("Configuration.StandardSequences.MagneticField_cff")

# Monte Carlo gun 
process.load("IOMC.FlatProtonLogKsiLogTGun.Beta2p5Energy3500GeV_cfi")
process.FlatProtonLogKsiLogTGun.Verbosity = 1
process.FlatProtonLogKsiLogTGun.MinKsi = -0.015
process.FlatProtonLogKsiLogTGun.MaxKsi = -0.25

# Smearing
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")
process.SmearingGenerator.verbosity = 0
process.SmearingGenerator.originalLabel = 'original'
process.SmearingGenerator.modifyLabel = 'source'

# Oscar - G4 simulation & proton transport
process.load("Configuration.TotemCommon.g4SimHits_cfi")
process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup
#bel 
# The input source for G4 module is connected to "process.source".
#
process.g4SimHits.Generator.HepMCProductLabel = 'source'    # energy+vertex smearing
#
# verbosity control
#
#process.g4SimHits.G4EventManagerVerbosity = 1
#process.g4SimHits.G4StackManagerVerbosity = 1
#process.g4SimHits.G4TrackingManagerVerbosity = 1
#process.g4SimHits.MagneticField.Verbosity = True
#process.g4SimHits.Physics.Verbosity = 1
#process.g4SimHits.Physics.BeamProtTransportSetup.Verbosity = True
#process.g4SimHits.Generator.Verbosity = 1
#process.g4SimHits.SteppingAction.Verbosity = 1
#process.g4SimHits.Totem_RP_SD.Verbosity = 1
#process.g4SimHits.TotemSD.Verbosity = 1

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
process.RPSingleTrackCandCollFit.Verbosity = 1

process.load("RecoTotemRP.RPInelasticReconstruction.Rec_3500GeV_beta_2p5_220_cfi")
process.RP220Reconst.BeamProtTransportSetup = process.BeamProtTransportSetup
# process.RP220Reconst.Verbosity = 1




################################################
#validation

process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_3500_Beta_11_220/RP_Dist_Beam_Cent.xml')

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_3500GeV_2p5_cfi")

# logging to txt files 
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# random number generators seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

########################### RP VALIDATION ##########################################

# module RPHitDists
process.load("TotemRPValidation.HitDistributions.HitDistributions_cfi")
process.RPHitDists.outputFile = cms.string('file:valRPinelasticBeta2.5Energy3.5TeVhits.root')

# module RPInelProtRecVal
process.load("TotemRPValidation.InelasticReconstructionValidation.InelasticReconstVal_2p5_3500GeV_cfi")
process.RPInelProtRecVal.HistogramFileName = 'file:valRPinelasticBeta2.5Energy3.5TeVplots.root'
process.RPInelProtRecVal.Verbosity = 1

# module RPRecTracksVal TODO
process.load("TotemRPValidation.RPReconstructedTracksValidation.ReconstTracksVal_beta_2p5_3500GeV_cfi")
process.RPRecTracksVal.HistogramFileName = 'file:valRPinelasticBeta2.5Energy3.5TeVtracks.root'

process.p1 = cms.Path(process.SmearingGenerator*process.g4SimHits*process.mix*process.RPSiDetDigitizer*process.RPClustProd*process.RPHecoHitProd*process.RPSinglTrackCandFind*process.RPSingleTrackCandCollFit*process.RP220Reconst*process.RPInelProtRecVal*process.RPHitDists*process.RPRecTracksVal)
process.outpath = cms.EndPath(process.o1)

