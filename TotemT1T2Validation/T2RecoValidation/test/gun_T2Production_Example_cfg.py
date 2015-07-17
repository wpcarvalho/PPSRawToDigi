import FWCore.ParameterSet.Config as cms

process = cms.Process("gunT1T2mu")

# Specify the maximum event to simulate
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Configure the output module (save the result in a file -- gunT1T2mu.root)
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
    fileName = cms.untracked.string('file:gunT1T2mu.root')
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

# Customized random number generator service used by
#   - Gun module (FlatRandomEGunProducer)
#   - G4 Simulation module (SimG4Core.Application.g4SimHits_cfi)
#   - Smearing module (Configuration.StandardSequences.VtxSmearedGauss_cff)
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    moduleSeeds = cms.PSet(
        generator = cms.untracked.uint32(456789),
        g4SimHits = cms.untracked.uint32(12340),
        VtxSmeared = cms.untracked.uint32(12340567)
    ),
    sourceSeed = cms.untracked.uint32(123456789)
)

# Load particle table
process.load("SimGeneral.HepPDTESSource.pdt_cfi")

# Smearing module
process.load("Configuration.StandardSequences.VtxSmearedGauss_cff")

########################### SIMU ###############################################

# Geometry
process.load("Configuration.TotemCommon.geometryT1T2_cfi")

### initialize magnetic field ###
# by default we have 3.8T
process.load("Configuration.StandardSequences.MagneticField_cff")

# NOTE:
#   In CMSSW_3_1_1, it is obligatory to add the "EmptySource" as the source before using
# the "FlatRandomEGunProducer" as the event producer. For CMSSW convention, this event
# producer will be labeled "generator.generator".
#   Ref: CMSSW/Configuration/EcalTB/python/reco_application_tbsim_DetSim-Digi_cfg.py, Line 57, 59
# -----------------------------------------------------------------------------------
#   In CMSSW_1_7_7. "FlatRandomEGunProducer" is directly used as the source and also
# labeled "generator.source".
process.source = cms.Source("EmptySource",
    firstRun = cms.untracked.uint32(1),
    firstEvent = cms.untracked.uint32(1)
)

process.generator = cms.EDProducer("FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(13),            # 13 (particle id) = mu^-
        MinEta = cms.double(5.0),
        MaxEta = cms.double(7.0),
        MinPhi = cms.double(-3.14159265358979323846),  # in radians
        MaxPhi = cms.double(3.14159265358979323846),
        MinE = cms.double(1.0),
        MaxE = cms.double(50.01)
    ),
    AddAntiParticle = cms.bool(False),
    Verbosity = cms.untracked.int32(0)  # set to 1 (or greater) for printouts
)

# Geant4-based CMS Detector simulation
process.load("Configuration.TotemCommon.g4SimHits_cfi")
process.g4SimHits.Watchers = cms.VPSet()
process.g4SimHits.Physics.type = cms.string('SimG4Core/Physics/QGSP_BERT_EML')
process.g4SimHits.Physics.DefaultCutValue = 100.
#
# NOTE:
#   In CMSSW_3_1_1, the Gun ("FlatRandomEGunProducer") is now created as a EDProducer and labeled
# as "process.generator". So the data source for Geant 4 simulation module (g4SimHits.Generator.HepMCProductLabel) 
# should now be connected to "process.generator".
# -----------------------------------------------------------------------------------
#   In CMSSW_1_7_7. "g4SimHits.Generator.HepMCProductLabel" is conntected to "process.source", because
# "FlatRandomEGunProducer" is directly used as the source and also labeled "generator.source".
#
process.g4SimHits.Generator.HepMCProductLabel = 'generator'
# process.g4SimHits.Generator.ApplyPtCuts = False (TypeError: ApplyPtCuts does not already exist)
process.g4SimHits.Generator.ApplyPCuts = False
process.g4SimHits.Generator.ApplyEtaCuts = False
process.g4SimHits.Generator.LeaveScatteredProtons = cms.bool(False)
process.g4SimHits.UseMagneticField = True


# Mixing module
process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.mix.playback = True
process.mix.useCurrentProcessOnly = True


########################### DIGI + RECO T2 ##########################################

process.load("SimTotem.T2Digitizer.T2Digis_cfi")

process.load("RecoTotemT1T2.T2MakeCluster.T2MakeCluster_cfi")

process.load("RecoTotemT1T2.T2RecHit.T2RecHit_cfi")

process.load("RecoTotemT1T2.T2RoadProducer.T2RoadCollDr15Dz200Dphi17_cfi")

process.load("RecoTotemT1T2.T2TrackProducer2.T2TrackColl2_cfi")

# NOTE:
# Remember to add process.generator (based on "pythia") to the path.
process.p1 = cms.Path(process.generator*process.VtxSmeared*process.g4SimHits*process.mix*process.T2Digis*process.T2MCl*process.T2Hits*process.T2RoadColl*process.T2TrackColl2)


process.outpath = cms.EndPath(process.o1)
