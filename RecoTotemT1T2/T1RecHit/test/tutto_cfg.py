import FWCore.ParameterSet.Config as cms

process = cms.Process("tutto")

process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('cout', 
        'warnings', 
        'errors', 
        'infos', 
        'debugs'),
    categories = cms.vstring('ForwardSim'),
    debugModules = cms.vstring('*'),
    infos = cms.PSet(
        threshold = cms.string('INFO')
    ),
    warnings = cms.PSet(
        threshold = cms.string('WARNING')
    ),
    errors = cms.PSet(
        threshold = cms.string('ERROR')
    ),
    debugs = cms.PSet(
        threshold = cms.string('DEBUG'),
        INFO = cms.PSet(
            limit = cms.int32(0)
        ),
        DEBUG = cms.PSet(
            limit = cms.int32(0)
        ),
        ForwardSim = cms.PSet(
            limit = cms.int32(1000000)
        ),
        CaloSim = cms.PSet(
            limit = cms.int32(0)
        )
    )
)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    moduleSeeds = cms.PSet(
        g4SimHits = cms.untracked.uint32(12340),
        VtxSmeared = cms.untracked.uint32(12340567)
    ),
    sourceSeed = cms.untracked.uint32(135790974)
)

process.load("SimGeneral.HepPDTESSource.pdt_cfi")

process.source = cms.Source("FlatRandomEGunSource",
    maxEvents = cms.untracked.int32(5000),
    PGunParameters = cms.untracked.PSet(
        PartID = cms.untracked.vint32(13),
        MinEta = cms.untracked.double(3.3),
        MaxEta = cms.untracked.double(4.5),
        MinPhi = cms.untracked.double(0.0),
        MaxPhi = cms.untracked.double(6.28),
        MinE = cms.untracked.double(99.99),
        MaxE = cms.untracked.double(100.01)
    ),
    Verbosity = cms.untracked.int32(0)
)

# Vtx Smearing
process.load("IOMC.EventVertexGenerators.VtxSmearedGauss_cfi")

# Geometry
process.load("SimG4CMS.Forward.test.testGeometryXML_cfi")

#Magnetic Field
### Full field map, static configuration for each field value
# process.load("Configuration.StandardSequences.MagneticField_20T_cff")
# process.load("Configuration.StandardSequences.MagneticField_30T_cff")
# process.load("Configuration.StandardSequences.MagneticField_35T_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
# process.load("Configuration.StandardSequences.MagneticField_40T_cff") 

process.o1 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('event.root')
)

process.g4SimHits = cms.EDProducer("OscarProducer",
    NonBeamEvent = cms.bool(False),
    G4EventManagerVerbosity = cms.untracked.int32(0),
    G4StackManagerVerbosity = cms.untracked.int32(0),
    G4TrackingManagerVerbosity = cms.untracked.int32(0),
    UseMagneticField = cms.bool(True),
    OverrideUserStackingAction = cms.bool(True),
    StoreRndmSeeds = cms.bool(False),
    RestoreRndmSeeds = cms.bool(False),
    PhysicsTablesDirectory = cms.string('PhysicsTables'),
    StorePhysicsTables = cms.bool(False),
    RestorePhysicsTables = cms.bool(False),

    Watchers = cms.VPSet(),

    MagneticField = cms.PSet(
        delta = cms.double(1.0)
    ),

    Physics = cms.PSet(
        # NOTE : if you want EM Physics only,
        #        please select "SimG4Core/Physics/DummyPhysics" for type
        #        and turn ON DummyEMPhysics
        #
        type = cms.string('SimG4Core/Physics/QGSP'),
        DummyEMPhysics = cms.bool(False),
        CutsPerRegion = cms.bool(True),
        DefaultCutValue = cms.double(10.0), # cuts in cm, i.e. 10m
        G4BremsstrahlungThreshold = cms.double(0.5),    # # cut in GeV
        Verbosity = cms.untracked.int32(0), # 1 will print cuts as they get set fdrom DD
                                            # 2 will do as 1 + will dump Geant4 table of cuts
        GFlashEmin = cms.double(1.0),
        GFlashEmax = cms.double(1000000.0),
        GFlashEToKill = cms.double(0.1)
    ),

    Generator = cms.PSet(
        HepMCProductLabel = cms.string('VtxSmeared'),   # alternative would be : 
                                                        # InputTag MCProductLabel = source
                                                        # to get the "original" on (without vtx smearing);
                                                        # NO other MCProduct's/labels in edm::Event
        ApplyPtCuts = cms.bool(False),
        ApplyEtaCuts = cms.bool(False),
        ApplyPhiCuts = cms.bool(False),
        MinPhiCut = cms.double(-3.14159265359),     # in radians
        MaxPhiCut = cms.double(3.14159265359),      # according to CMS conventions
        MinEtaCut = cms.double(-5.5),
        MaxEtaCut = cms.double(5.5),
        MinPtCut = cms.double(4e-05),               # the pt-cut is in GeV (CMS conventions)
        MaxPtCut = cms.double(99999.0),             # the ptmax=99.TeV in this case
        DecLenCut = cms.untracked.double(10.0),     # the minimum decay length in cm (!) for mother tracking
        Verbosity = cms.untracked.int32(0)
    ),

    RunAction = cms.PSet(
        StopFile = cms.string('StopRun')
    ),

    EventAction = cms.PSet(
        CollapsePrimaryVertices = cms.bool(False),
        StopFile = cms.string('StopRun'),
        debug = cms.untracked.bool(False)
    ),

    StackingAction = cms.PSet(
        SavePrimaryDecayProductsAndConversions = cms.untracked.bool(False)
    ),

    TrackingAction = cms.PSet(
        DetailedTiming = cms.untracked.bool(False)
    ),

    SteppingAction = cms.PSet(
        KillBeamPipe = cms.bool(True),
        CriticalEnergyForVacuum = cms.double(2.0),
        CriticalDensity = cms.double(1e-15),
        Verbosity = cms.untracked.int32(0)
    ),
    
    TrackerSD = cms.PSet(
        ZeroEnergyLoss = cms.bool(False),
        NeverAccumulate = cms.bool(False),
        PrintHits = cms.bool(False),
        ElectronicSigmaInNanoSeconds = cms.double(12.06),
        EnergyThresholdForPersistencyInGeV = cms.double(0.5),
        EnergyThresholdForHistoryInGeV = cms.double(0.05)
    ),
    
    MuonSD = cms.PSet(
        EnergyThresholdForPersistency = cms.double(1.0),
        AllMuonsPersistent = cms.bool(False),
        PrintHits = cms.bool(False)
    ),
    
    CaloSD = cms.PSet(
        EminTrack = cms.double(1.0),
        CheckHits = cms.untracked.int32(25),
        UseMap = cms.untracked.bool(True),
        Verbosity = cms.untracked.int32(0),
        DetailedTiming = cms.untracked.bool(False)
    ),
    
    
    ECalSD = cms.PSet(
        UseBirkLaw = cms.bool(False),
        BirkC1 = cms.double(0.013),
        BirkC2 = cms.double(9.6e-06)
    ),

    HCalSD = cms.PSet(
        UseBirkLaw = cms.bool(False),
        BirkC1 = cms.double(0.013),
        BirkC2 = cms.double(9.6e-06),
        UseShowerLibrary = cms.bool(True),
        TestNumberingScheme = cms.bool(False),
        UseHF = cms.untracked.bool(True),
        forTBH2 = cms.untracked.bool(False),
    ),

    CaloTrkProcessing = cms.PSet(
        TestBeam = cms.bool(False),
        EminTrack = cms.double(0.01)
    ),

    HFShower = cms.PSet(
        ProbMax = cms.double(0.7268),
        CFibre = cms.double(0.5)
    ),

    HFShowerLibrary = cms.PSet(
        FileName = cms.FileInPath('SimG4CMS/Calo/data/hfshowerlibrary_lhep.root'),
        TreeEMID = cms.string('h3'),
        TreeHadID = cms.string('h8')
    ),

    HFCherenkov = cms.PSet(
        RefIndex = cms.double(1.459),
        Lambda1 = cms.double(280.0),
        Lambda2 = cms.double(700.0),
        Aperture = cms.double(0.33),
        Gain = cms.double(0.33),
        ApertureTrapped = cms.double(0.22),
        CheckSurvive = cms.bool(False),
    ),

    CastorSD = cms.PSet(
        Verbosity = cms.untracked.int32(0)
    ),

    TotemSD = cms.PSet(
        Verbosity = cms.untracked.int32(0)
    ),
    
    ZdcSD = cms.PSet(
        Verbosity = cms.int32(0),
        FiberDirection = cms.double(0.0)
    ),
    
    HcalTB06BeamSD = cms.PSet(
        UseBirkLaw = cms.bool(False),
        BirkC1 = cms.double(0.013),
        BirkC2 = cms.double(9.6e-06)
    ) 
)

process.load("SimTotem.T1Digitizer.T1DigisVFAT_cfi")

process.t1cluster = cms.EDFilter("T1MakeCluster",
    Electronics = cms.string('VFAT'),
    WEIGHT1 = cms.double(500.0),
    WEIGHT2 = cms.double(500.0)
)

process.t1rechit = cms.EDProducer("T1RecHitProducer",
    ChiSquare = cms.double(2.0)
)

process.t1recanal = cms.EDAnalyzer("T1GhostAnalyzer",
    outFile = cms.string('histoGh.root'),
    Event2Plot = cms.int32(3)
)

process.p1 = cms.Path(process.VtxSmeared*process.g4SimHits*process.T1Digis*process.t1cluster*process.t1rechit*process.t1recanal)

process.outpath = cms.EndPath(process.o1)

