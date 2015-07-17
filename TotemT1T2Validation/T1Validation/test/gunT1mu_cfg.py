import FWCore.ParameterSet.Config as cms

process = cms.Process("gunT1T2mu")

# Specify the maximum event to simulate
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
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

process.load("Configuration.TotemCommon.LoggerMin_cfi")

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    moduleSeeds = cms.PSet(
        generator = cms.untracked.uint32(456789),
        g4SimHits = cms.untracked.uint32(12340),
        SmearingGenerator = cms.untracked.uint32(12340567)
    ),
    sourceSeed = cms.untracked.uint32(123456789)
)

# Load particle table
process.load("SimGeneral.HepPDTESSource.pdt_cfi")

# Smearing module
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")

# optics
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_3500GeV_2p5_cfi")

########################### SIMU ###############################################

# Geometry
process.load("Configuration.TotemCommon.geometryT1T2_cfi")

### initialize magnetic field ###
# by default we have 3.8T
process.load("Configuration.StandardSequences.MagneticField_cff")

process.source = cms.Source("EmptySource",
    firstRun = cms.untracked.uint32(1),
    firstEvent = cms.untracked.uint32(1)
)

process.generator = cms.EDProducer("FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(13),            # 13 (particle id) = mu^-
        MinEta = cms.double(3.99),
        MaxEta = cms.double(4.3),
        MinPhi = cms.double(2.65),  #  -3.14159265358979323846),  # in radians
        MaxPhi = cms.double(3.65),  #   3.14159265358979323846),
        MinE = cms.double(50.180),
        MaxE = cms.double(50.1801)
    ),
    AddAntiParticle = cms.bool(False),
    Verbosity = cms.untracked.int32(0)  # set to 1 (or greater) for printouts
)

# Geant4-based CMS Detector simulation
process.load("Configuration.TotemCommon.g4SimHits_cfi")
process.g4SimHits.Watchers = cms.VPSet()
process.g4SimHits.Physics.type = cms.string('SimG4Core/Physics/QGSP_BERT_EML')
process.g4SimHits.Physics.DefaultCutValue = 100.
process.g4SimHits.Generator.HepMCProductLabel = 'generator'
process.g4SimHits.Generator.ApplyPCuts = False
process.g4SimHits.Generator.ApplyEtaCuts = False
process.g4SimHits.Generator.LeaveScatteredProtons = cms.bool(False)
process.g4SimHits.UseMagneticField = True


# Mixing module
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

########################### DIGI + RECO T1 ##########################################

process.load("SimTotem.T1Digitizer.T1DigisVFAT_cfi")

process.load("RecoTotemT1T2.T1MakeCluster.T1MakeCluster_cfi")

#process.load("RecoTotemT1T2.T1RecHit.T1RecHit_cfi")

process.t1rechit = cms.EDProducer("T1RecHitProducer",
    Verbosity = cms.int32(0),
    ChiSquare = cms.double(2.0)
)

process.load("RecoTotemT1T2.T1RoadProducer.T1RoadProducer_cfi")

process.load("RecoTotemT1T2.T1TrackProducer.T1TrackProducer_cfi")

process.t1valid = cms.EDAnalyzer("T1Validation",
    Verbosity = cms.int32(0),
    SIM = cms.double(1.0),
    DIGI = cms.double(1.0),
    RECO = cms.double(1.0),
    TrackLabel = cms.string('t1tracks'),
    OutputFile = cms.string('./T1Valid_gunmu.root')
)

process.T1RecoAnal = cms.EDAnalyzer("T1RecHitAnalyzer",
                                    OutputFile = cms.string('./T1RecoAnal.root') )

process.p1 = cms.Path(process.generator*process.SmearingGenerator*process.g4SimHits*process.mix*process.T1Digis*process.t1cluster*process.t1rechit*process.t1roads*process.t1tracks*process.t1valid*process.T1RecoAnal)

process.outpath = cms.EndPath(process.o1)
