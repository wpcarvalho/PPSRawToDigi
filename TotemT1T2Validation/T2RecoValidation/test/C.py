import FWCore.ParameterSet.Config as cms



process = cms.Process("RecoTrk2011April") #remind to change name 







process.maxEvents = cms.untracked.PSet(

    input = cms.untracked.int32(20)  #-1

)













process.o1 = cms.OutputModule("PoolOutputModule",                 

       fileName = cms.untracked.string('file:run_0_.root'))





# Use particle table

process.load("SimGeneral.HepPDTESSource.pdt_cfi")







process.source = cms.Source("EmptySource")



process.generator = cms.EDFilter("Pythia8GeneratorFilter",

     maxEventsToPrint = cms.untracked.int32(1),

     pythiaPylistVerbosity = cms.untracked.int32(1),

     filterEfficiency = cms.untracked.double(1.0),

     pythiaHepMCVerbosity = cms.untracked.bool(False),

     comEnergy = cms.double(8000.),

     #useUserHook = cms.bool(True),

     PythiaParameters = cms.PSet(

        processParameters = cms.vstring(

            'Main:timesAllowErrors    = 10000',

            'ParticleDecays:limitTau0 = on',

            'ParticleDecays:tauMax = 10',

            'SoftQCD:minBias = on',

            'SoftQCD:singleDiffractive = on',

            'SoftQCD:doubleDiffractive = on',

            'Tune:pp 2',

            'Tune:ee 3'),

        parameterSets = cms.vstring('processParameters')

    )

)

 









process.load("SimTotem.T2Digitizer.T2DigisSTDGeoEffi_cfi")





#process.load("SimTotem.T2Digitizer.T2DigisOnlyPrim_cfi")

process.T2Digis.saveDigiVFAT=cms.bool(True)

process.T2Digis.Misalignment.simulatemisalign=cms.untracked.bool(True)#UseAlign

#Decomment for alignment multi-scenario validation

process.T2Digis.Misalignment.inputFileNameMisal="SimTotem/T2Digitizer/data/2012IntGlobAlignV1.dat"#DigiMisalFilename

  















process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",

    moduleSeeds = cms.PSet(

        generator = cms.untracked.uint32(47844),

        g4SimHits = cms.untracked.uint32(15391),

        VtxSmeared = cms.untracked.uint32(22271),

        T2Digis = cms.untracked.uint32(94682),

        T2MCl = cms.untracked.uint32(78164),

    ),

    sourceSeed = cms.untracked.uint32(16785),

   

)





#from SimGeneral.MixingModule.mixNoPU_cfi  import *

process.load("SimGeneral.MixingModule.mixNoPU_cfi")



#process.load("Geometry.CMSCommonData.cmsExtendedGeometryGFlashXMLMy_cfi.py")

process.load("Configuration.StandardSequences.MagneticField_cff")



# G4 simulation & proton transport Configuration.TotemCommon.cmsIdealGeometryGFlashXML_cfi 

#process.load("SimG4Core.Application.g4SimHitsMy_cfi")





process.XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",

    geomXMLFiles = cms.vstring('Geometry/CMSCommonData/data/materials.xml',

        'Geometry/CMSCommonData/data/rotations.xml',

        'Geometry/CMSCommonData/data/extend/cmsextent.xml',

        'Geometry/CMSCommonData/data/cms.xml',

        'Geometry/CMSCommonData/data/cmsMother.xml',

        'Geometry/CMSCommonData/data/caloBase.xml',

        'Geometry/CMSCommonData/data/cmsCalo.xml',

        'Geometry/CMSCommonData/data/muonBase.xml',

        'Geometry/TrackerCommonData/data/trackermaterial.xml',

        'Geometry/CMSCommonData/data/mgnt.xml',



        'Geometry/CMSCommonData/data/cmsMuon.xml',



                    

                               

  #      'Geometry/CMSCommonData/data/beampipe.xml',

        'Geometry/ForwardCommonData/data/beampipe.xml',                       

        'Geometry/CMSCommonData/data/cmsBeam.xml',

                               

        'Geometry/CMSCommonData/data/muonMB.xml',

        'Geometry/CMSCommonData/data/muonMagnet.xml',





        'Geometry/CMSCommonData/data/cavern.xml',

        'Geometry/TrackerCommonData/data/pixfwdMaterials.xml',

        'Geometry/TrackerCommonData/data/pixfwdCommon.xml',

        'Geometry/TrackerCommonData/data/pixfwdPlaq.xml',

        'Geometry/TrackerCommonData/data/pixfwdPlaq1x2.xml',

        'Geometry/TrackerCommonData/data/pixfwdPlaq1x5.xml',

        'Geometry/TrackerCommonData/data/pixfwdPlaq2x3.xml',

        'Geometry/TrackerCommonData/data/pixfwdPlaq2x4.xml',

        'Geometry/TrackerCommonData/data/pixfwdPlaq2x5.xml',

        'Geometry/TrackerCommonData/data/pixfwdPanelBase.xml',

        'Geometry/TrackerCommonData/data/pixfwdPanel.xml',

        'Geometry/TrackerCommonData/data/pixfwdBlade.xml',

        'Geometry/TrackerCommonData/data/pixfwdNipple.xml',

        'Geometry/TrackerCommonData/data/pixfwdDisk.xml',

        'Geometry/TrackerCommonData/data/pixfwdCylinder.xml',

        'Geometry/TrackerCommonData/data/pixfwd.xml',

        'Geometry/TrackerCommonData/data/pixbarmaterial.xml',

        'Geometry/TrackerCommonData/data/pixbarladder.xml',

        'Geometry/TrackerCommonData/data/pixbarladderfull.xml',

        'Geometry/TrackerCommonData/data/pixbarladderhalf.xml',

        'Geometry/TrackerCommonData/data/pixbarlayer.xml',

        'Geometry/TrackerCommonData/data/pixbarlayer0.xml',

        'Geometry/TrackerCommonData/data/pixbarlayer1.xml',

        'Geometry/TrackerCommonData/data/pixbarlayer2.xml',

        'Geometry/TrackerCommonData/data/pixbar.xml',

        'Geometry/TrackerCommonData/data/trackermaterial.xml',

        'Geometry/TrackerCommonData/data/tracker.xml',

        'Geometry/TrackerCommonData/data/trackerpixbar.xml',

        'Geometry/TrackerCommonData/data/trackerpixfwd.xml',

 #       'Geometry/TrackerCommonData/data/trackertibtidservices.xml',

 #       'Geometry/TrackerCommonData/data/trackertib.xml',

 #       'Geometry/TrackerCommonData/data/trackertid.xml',

 #       'Geometry/TrackerCommonData/data/trackertob.xml',

 #       'Geometry/TrackerCommonData/data/trackertec.xml',

 #        'Geometry/TrackerCommonData/data/trackerbulkhead.xml',

        'Geometry/TrackerCommonData/data/trackerother.xml',

        'Geometry/HcalSimData/data/hf.xml', 

        'Geometry/HcalSimData/data/hffibre.xml',                       



                               

#        'Geometry/CMSCommonData/data/cavern.xml',

        'Geometry/ForwardCommonData/data/forward.xml',

        'Geometry/ForwardCommonData/data/forwardshield.xml',

                      

#        'Geometry/ForwardCommonData/data/bundle/forwardshield.xml',

#        'Geometry/ForwardCommonData/data/brmrotations.xml',

#        'Geometry/ForwardCommonData/data/brm.xml',

         'Geometry/ForwardCommonData/data/totemMaterials.xml',

         'Geometry/ForwardCommonData/data/totemRotations.xml',

         'Geometry/ForwardCommonData/data/totemt1.xml',

         'Geometry/ForwardCommonData/data/totemt2.xml',

         'Geometry/ForwardCommonData/data/ionpump.xml',

         'Geometry/ForwardCommonData/data/castor.xml',

         'Geometry/ForwardSimData/data/totemsensGem.xml',







         'Geometry/HcalCommonData/data/hcalrotations.xml', 

         'Geometry/HcalCommonData/data/hcalalgo.xml', 

         'Geometry/HcalCommonData/data/hcalbarrelalgo.xml', 

         'Geometry/HcalCommonData/data/hcalendcapalgo.xml', 

         'Geometry/HcalCommonData/data/hcalouteralgo.xml', 

         'Geometry/HcalCommonData/data/hcalforwardalgo.xml', 

         'Geometry/HcalCommonData/data/hcalforwardfibre.xml', 

         'Geometry/HcalCommonData/data/hcalforwardmaterial.xml',

         'Geometry/MuonCommonData/data/mbCommon.xml', 

         'Geometry/MuonCommonData/data/mb1.xml', 

         'Geometry/MuonCommonData/data/mb2.xml', 

         'Geometry/MuonCommonData/data/mb3.xml', 

         'Geometry/MuonCommonData/data/mb4.xml', 

         'Geometry/MuonCommonData/data/muonYoke.xml', 

         'Geometry/MuonCommonData/data/mf.xml',                      





                               

 #        'Geometry/HcalSimData/data/CaloUtil.xml', 

  #       'Geometry/HcalSimData/data/hf.xml', 

#         'Geometry/HcalSimData/data/hffibre.xml', 

                               

#        'Geometry/ForwardCommonData/data/lumimaterials.xml',

#        'Geometry/ForwardCommonData/data/zdcrotations.xml',

#        'Geometry/ForwardCommonData/data/lumirotations.xml',

#        'Geometry/ForwardCommonData/data/zdc.xml',

#        'Geometry/ForwardCommonData/data/zdclumi.xml',      

#        'Geometry/HcalCommonData/data/hcalsenspmf.xml',

#        'Geometry/HcalSimData/data/hf.xml',

#       'Geometry/HcalSimData/data/hfpmt.xml',

#        'Geometry/HcalSimData/data/hffibrebundle.xml',

#        'Geometry/HcalSimData/data/CaloUtil.xml',

#        'Geometry/ForwardCommonData/data/brmsens.xml',

#        'SimG4Core/GFlash/data/gflashCaloProdCuts.xml',

 #       'Geometry/EcalSimData/data/ESProdCuts.xml'

 #       'Geometry/TrackerSimData/data/trackerProdCuts.xml',

 #       'Geometry/TrackerSimData/data/trackerProdCutsBEAM.xml',

 #       'Geometry/MuonSimData/data/muonProdCuts.xml',

 #       'Geometry/ForwardSimData/data/zdcProdCuts.xml',

 #       'Geometry/ForwardSimData/data/ForwardShieldProdCuts.xml'

        ),

    rootNodeName = cms.string('cms:OCMS')

)











from SimG4Core.Application.hectorParameter_cfi import *



common_heavy_suppression = cms.PSet(

    NeutronThreshold = cms.double(30.0),

    ProtonThreshold = cms.double(30.0),

    IonThreshold = cms.double(30.0)

)



#common_maximum_time = cms.PSet(

#    MaxTrackTime  = cms.double(500.0),

#    MaxTimeNames  = cms.vstring('ZDCRegion','QuadRegion','InterimRegion'),

#    MaxTrackTimes = cms.vdouble(2000.0,0.,0.)

#)

common_maximum_time = cms.PSet(

    MaxTrackTime  = cms.double(1000.0),

    MaxTimeNames  = cms.vstring('ZDCRegion','CastorRegion','QuadRegion','InterimRegion'),

    MaxTrackTimes = cms.vdouble(2000.0,2000.,2000.,2000.)

)

common_UsePMT = cms.PSet(

    UseR7600UPMT  = cms.bool(False)

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

    CheckOverlap = cms.untracked.bool(False),

    G4Commands = cms.vstring(),

    Watchers = cms.VPSet(),

    theLHCTlinkTag = cms.InputTag("LHCTransport"),

    MagneticField = cms.PSet(

        UseLocalMagFieldManager = cms.bool(False),

        Verbosity = cms.untracked.bool(False),

        ConfGlobalMFM = cms.PSet(

            Volume = cms.string('OCMS'),

            OCMS = cms.PSet(

                Stepper = cms.string('G4ClassicalRK4'),

                Type = cms.string('CMSIMField'),

                G4ClassicalRK4 = cms.PSet(

                    MaximumEpsilonStep = cms.untracked.double(0.01), ## in mm



                    DeltaOneStep = cms.double(0.001), ## in mm



                    MaximumLoopCounts = cms.untracked.double(1000.0),

                    DeltaChord = cms.double(0.001), ## in mm



                    MinStep = cms.double(0.1), ## in mm



                    DeltaIntersectionAndOneStep = cms.untracked.double(-1.0),

                    DeltaIntersection = cms.double(0.0001), ## in mm



                    MinimumEpsilonStep = cms.untracked.double(1e-05) ## in mm



                )

            )

        ),

        delta = cms.double(1.0)

    ),

    Physics = cms.PSet(

        # NOTE : if you want EM Physics only,

        #        please select "SimG4Core/Physics/DummyPhysics" for type

        #        and turn ON DummyEMPhysics

        #

        type = cms.string('SimG4Core/Physics/QGSP_FTFP_BERT_EML'),

        DummyEMPhysics = cms.bool(False),

        CutsPerRegion = cms.bool(False),

        DefaultCutValue = cms.double(0.001), ## cuts in cm 1.0 orig

        G4BremsstrahlungThreshold = cms.double(0.5), ## cut in GeV

        Verbosity = cms.untracked.int32(0),

        # 1 will print cuts as they get set from DD

        # 2 will do as 1 + will dump Geant4 table of cuts

        MonopoleCharge       = cms.untracked.int32(1),

        MonopoleDeltaRay     = cms.untracked.bool(True),

        MonopoleMultiScatter = cms.untracked.bool(False),

        MonopoleTransport    = cms.untracked.bool(True),

        Region      = cms.string(' '),

	TrackingCut = cms.bool(True),

        SRType      = cms.bool(True),

        EMPhysics   = cms.untracked.bool(True),

        HadPhysics  = cms.untracked.bool(True),

        FlagBERT    = cms.untracked.bool(False),

        FlagCHIPS   = cms.untracked.bool(False),

        FlagFTF     = cms.untracked.bool(False),

        FlagGlauber = cms.untracked.bool(False),

        FlagHP      = cms.untracked.bool(False),

        GflashEcal  = cms.bool(False),

        bField      = cms.double(3.8),

        energyScaleEB = cms.double(1.032),

        energyScaleEE = cms.double(1.024),

        GflashHcal  = cms.bool(False)

        #GFlash = cms.PSet(

        #    GflashHistogram = cms.bool(True),

        #    GflashEMShowerModel = cms.bool(True),

        #    GflashHadronPhysics = cms.string('QGSP_BERT_EMV'),

        #    GflashHadronShowerModel = cms.bool(False)

        #)

    ),

    Generator = cms.PSet(

        HectorEtaCut,

        # string HepMCProductLabel = "VtxSmeared"

        HepMCProductLabel = cms.string('generator'),

        ApplyPCuts = cms.bool(False),

        MinPCut = cms.double(0.04), ## the pt-cut is in GeV (CMS conventions)

        MaxPCut = cms.double(99999.0), ## the ptmax=99.TeV in this case

        ApplyEtaCuts = cms.bool(True),

        MinEtaCut = cms.double(-5.5),

        MaxEtaCut = cms.double(5.5),

        ApplyPhiCuts = cms.bool(False),

        MinPhiCut = cms.double(-3.14159265359), ## in radians

        MaxPhiCut = cms.double(3.14159265359), ## according to CMS conventions

        RDecLenCut = cms.double(2.9), ## the minimum decay length in cm (!) for mother tracking

        Verbosity = cms.untracked.int32(0)

    ),

    RunAction = cms.PSet(

        StopFile = cms.string('StopRun')

    ),

    EventAction = cms.PSet(

        debug = cms.untracked.bool(False),

        StopFile = cms.string('StopRun'),

        CollapsePrimaryVertices = cms.bool(False)

    ),

    StackingAction = cms.PSet(

        common_heavy_suppression,

        common_maximum_time,

        KillDeltaRay  = cms.bool(False),

        TrackNeutrino = cms.bool(False),

        KillHeavy     = cms.bool(False),

        SaveFirstLevelSecondary = cms.untracked.bool(True),

        SavePrimaryDecayProductsAndConversionsInTracker = cms.untracked.bool(True),

        SavePrimaryDecayProductsAndConversionsInCalo = cms.untracked.bool(False),

        SavePrimaryDecayProductsAndConversionsInMuon = cms.untracked.bool(False)

    ),

    TrackingAction = cms.PSet(

        DetailedTiming = cms.untracked.bool(False)

    ),

    SteppingAction = cms.PSet(

        common_maximum_time,

        KillBeamPipe            = cms.bool(True),

        CriticalEnergyForVacuum = cms.double(2.0),

        CriticalDensity         = cms.double(1e-15),

        EkinNames               = cms.vstring(),

        EkinThresholds          = cms.vdouble(),

        EkinParticles           = cms.vstring(),

        Verbosity = cms.untracked.int32(0)

    ),

    

  

   

    HFShower = cms.PSet(

        common_UsePMT,

        ProbMax         = cms.double(1.0),

        CFibre          = cms.double(0.5),

        PEPerGeV        = cms.double(0.25),#(0.31)

        TrackEM         = cms.bool(False),

        UseShowerLibrary= cms.bool(False),

        UseHFGflash     = cms.bool(True),

        EminLibrary     = cms.double(0.0),

        RefIndex        = cms.double(1.459),

        Lambda1         = cms.double(280.0),

        Lambda2         = cms.double(700.0),

        Aperture        = cms.double(0.33),

        ApertureTrapped = cms.double(0.22),

        Gain            = cms.double(0.33),

        OnlyLong        = cms.bool(True),

        LambdaMean      = cms.double(350.0),

        CheckSurvive    = cms.bool(False),

        ApplyFiducialCut= cms.bool(True),

        ParametrizeLast = cms.untracked.bool(False)

    ),

    HFShowerLibrary = cms.PSet(

        FileName        = cms.FileInPath('SimG4CMS/Calo/data/hfshowerlibrary_lhep_140_edm.root'),

        BackProbability = cms.double(0.2),

        TreeEMID        = cms.string('emParticles'),

        TreeHadID       = cms.string('hadParticles'),

        Verbosity       = cms.untracked.bool(False),

        ApplyFiducialCut= cms.bool(True),

        BranchPost      = cms.untracked.string('_R.obj'),

        BranchEvt       = cms.untracked.string('HFShowerLibraryEventInfos_hfshowerlib_HFShowerLibraryEventInfo'),

        BranchPre       = cms.untracked.string('HFShowerPhotons_hfshowerlib_')

    ),

    HFShowerPMT = cms.PSet(

        common_UsePMT,

        PEPerGeVPMT     = cms.double(1.0),

        RefIndex        = cms.double(1.52),

        Lambda1         = cms.double(280.0),

        Lambda2         = cms.double(700.0),

        Aperture        = cms.double(0.99),

        ApertureTrapped = cms.double(0.22),

        Gain            = cms.double(0.33),

        CheckSurvive    = cms.bool(False)

    ),

    HFShowerStraightBundle = cms.PSet(

        common_UsePMT,

        FactorBundle    = cms.double(1.0),

        RefIndex        = cms.double(1.459),

        Lambda1         = cms.double(280.0),

        Lambda2         = cms.double(700.0),

        Aperture        = cms.double(0.33),

        ApertureTrapped = cms.double(0.22),

        Gain            = cms.double(0.33),

        CheckSurvive    = cms.bool(False)

    ),

    HFShowerConicalBundle = cms.PSet(

        common_UsePMT,

        FactorBundle    = cms.double(1.0),

        RefIndex        = cms.double(1.459),

        Lambda1         = cms.double(280.0),

        Lambda2         = cms.double(700.0),

        Aperture        = cms.double(0.33),

        ApertureTrapped = cms.double(0.22),

        Gain            = cms.double(0.33),

        CheckSurvive    = cms.bool(False)

    ),

    HFGflash = cms.PSet(

        BField          = cms.untracked.double(3.8),

        WatcherOn       = cms.untracked.bool(True),

        FillHisto       = cms.untracked.bool(True)

    ),

   

   # CastorShowerLibrary =  cms.PSet(

   #     FileName  = cms.FileInPath('SimG4CMS/Forward/data/CastorShowerLibrary_CMSSW500_Standard.root'),

   #     BranchEvt = cms.untracked.string('hadShowerLibInfo.'),

    #    BranchEM  = cms.untracked.string('emParticles.'),

   #     BranchHAD = cms.untracked.string('hadParticles.'),

   #     Verbosity = cms.untracked.bool(False)

  #  ),

    TotemSD = cms.PSet(

        Verbosity = cms.untracked.int32(0)

    ),

    BscSD = cms.PSet(

         Verbosity = cms.untracked.int32(0)

    ),

    CaloTrkProcessing = cms.PSet(

         TestBeam   = cms.bool(False),

         EminTrack  = cms.double(0.01),

         PutHistory = cms.bool(False)

    ),

   CaloSD = cms.PSet(

        common_heavy_suppression,

        SuppressHeavy = cms.bool(False),

        EminTrack = cms.double(1.0),

        TmaxHit   = cms.double(1000.0),

        HCNames   = cms.vstring('EcalHitsEB','EcalHitsEE','EcalHitsES','HcalHits','ZDCHITS'),

        EminHits  = cms.vdouble(0.015,0.010,0.0,0.0,0.0),

        EminHitsDepth = cms.vdouble(0.0,0.0,0.0,0.0,0.0),

        TmaxHits  = cms.vdouble(500.0,500.0,500.0,500.0,2000.0),

        UseResponseTables = cms.vint32(0,0,0,0,0),

        BeamPosition      = cms.double(0.0),

        CorrectTOFBeam    = cms.bool(False),

        DetailedTiming    = cms.untracked.bool(False),

        UseMap            = cms.untracked.bool(False),

        Verbosity         = cms.untracked.int32(0),

        CheckHits         = cms.untracked.int32(25)

    ),                                

    ECalSD = cms.PSet(

        UseBirkLaw      = cms.bool(True),

        BirkL3Parametrization = cms.bool(True),

        BirkSlope       = cms.double(0.253694),

        BirkCut         = cms.double(0.1),

        BirkC1          = cms.double(0.03333),

        BirkC3          = cms.double(1.0),

        BirkC2          = cms.double(0.0),

        SlopeLightYield = cms.double(0.02),

        StoreSecondary  = cms.bool(False),

        TimeSliceUnit   = cms.int32(1),

        IgnoreTrackID   = cms.bool(False),

        XtalMat         = cms.untracked.string('E_PbWO4'),

        TestBeam        = cms.untracked.bool(False),

        NullNumbering   = cms.untracked.bool(False),

        StoreRadLength  = cms.untracked.bool(False)

    ),

    MuonSD = cms.PSet(

        EnergyThresholdForPersistency = cms.double(1.0),

        PrintHits = cms.bool(False),

        AllMuonsPersistent = cms.bool(True)

    ),

     TrackerSD = cms.PSet(

        ZeroEnergyLoss = cms.bool(False),

        PrintHits = cms.bool(False),

        ElectronicSigmaInNanoSeconds = cms.double(12.06),

        NeverAccumulate = cms.bool(False),

        EnergyThresholdForPersistencyInGeV = cms.double(0.2),

        EnergyThresholdForHistoryInGeV = cms.double(0.05)

    ),                               

   HCalSD = cms.PSet(

        UseBirkLaw          = cms.bool(True),

        BirkC3              = cms.double(1.75),

        BirkC2              = cms.double(0.142),

        BirkC1              = cms.double(0.0052),

        UseShowerLibrary    = cms.bool(False),

        UseParametrize      = cms.bool(True),

        UsePMTHits          = cms.bool(True),

        UseFibreBundleHits  = cms.bool(True),

        TestNumberingScheme = cms.bool(False),

        EminHitHB           = cms.double(0.0),

        EminHitHE           = cms.double(0.0),

        EminHitHO           = cms.double(0.0),

        EminHitHF           = cms.double(0.0),

        BetaThreshold       = cms.double(0.7),

        TimeSliceUnit       = cms.int32(1),

        IgnoreTrackID       = cms.bool(False),

        UseHF               = cms.untracked.bool(True),

        ForTBH2             = cms.untracked.bool(False),

        UseLayerWt          = cms.untracked.bool(False),

        WtFile              = cms.untracked.string('None')

    )                                

)

















process.load("RecoTotemT1T2.T2MakeCluster.T2MakeCluster_cfi")







process.T2MCl.TakeCleanEventOnly=cms.bool(False) #IMPORTANT

#process.T2MCl.maskvect = cms.vuint32(21,23,25,27,29)

process.T2MCl.SimuClusterEfficiency=cms.bool(True)#SimulateClusterEffi

process.T2MCl.EffiGeoRootFileName= cms.string("RecoTotemT1T2/T2MakeCluster/data/Effi-8334.root")



process.load("RecoTotemT1T2.T2RecHit.T2RecHit_cfi")







#from RecoTotemT1T2.T2RecHit.T2RecHit_cfi import *

process.T2Hits.Cl1MaxPad = cms.uint32(25) #Tune better

process.T2Hits.Cl1MaxStrip = cms.uint32(25)

process.T2Hits.IncludeClass0Hits = True

#T2Hits.inputFileNameMisal=cms.untracked.string('/afs/cern.ch/exp/totem/scratch/berretti/tmp/testSplitMerge/CMSSW_3_1_1/src/$

#process.T2Hits.inputFileNameMisal=cms.untracked.string('SimTotem/T2Digitizer/data/run_5525-5535_IntAlignHX50000Evt_XovY0.3HIP_ANDGLOB_BBConf.dat')

process.T2Hits.inputFileNameMisal=cms.untracked.string('SimTotem/T2Digitizer/data/2012IntGlobAlignV1.dat')#DigiFNameForHit



process.T2Hits.useTXTfile=cms.bool(True) #UseAlign

process.T2Hits.InsertAlignmentbyCFG=cms.bool(True) # UseAlign True for data

process.T2Hits.verbosity=cms.untracked.bool(False)

process.T2Hits.CorrectWithResolution=cms.bool(True) #False:Old Strategy True:New Strategy



#from RecoTotemT1T2.T2RoadPadFinder.NewLabelT2RoadPadFinder_cfi import *

process.load("RecoTotemT1T2.T2RoadPadFinder.NewLabelT2RoadPadFinder_cfi")#T2RoadPadFinder_cfi



process.T2RoadPadFinder.HitLabel=cms.string("T2Hits")#T2HitsSTD for effi reco |  #T2Hits for std reco

process.T2RoadPadFinder.CluLabel=cms.string("T2MCl")#T2MClSTD for effi reco |   T2MCl  for std reco



process.T2RoadPadFinder.verbosity = 0

process.T2RoadPadFinder.TwoPointsTubesAngularCollinearity=0.09#0.07 default

process.T2RoadPadFinder.MinCluSize_considered_asBlobs = cms.int32(5)

process.T2RoadPadFinder.MinimumNumCl1Hit= 3

process.T2RoadPadFinder.chi2XYProb_Thr= 0.01

process.T2RoadPadFinder.Nmin_padsFinal= 4

process.T2RoadPadFinder.T2RoadCollProdName="NewRoadFinderRELOAD"







process.T2RoadPadFinder.AllowsPadReAssociation=False

process.T2RoadPadFinder.AllowsConcurrentBranches=False

process.T2RoadPadFinder.useStraightPadTowers= cms.bool(True)#False for internal alignment studies

process.T2RoadPadFinder.ResolveOverlapDoubleCount = cms.bool(False) #Default is True, False for shadow alignment and dndeta An

process.T2RoadPadFinder.OverlapDoubleCountDR = cms.double(2.0) #Depend on your alignment Resol

process.T2RoadPadFinder.OverlapDoubleCountDPhi =cms.double(3.5)

process.T2RoadPadFinder.OverlapDoubleCountDTheta =  cms.double(0.01)



#process.T2RoadPadFinder.VplaneToExclude = cms.vint32(9,19,29,39)





process.T2RoadPadFinder.QuartedSelected = cms.vint32(0,1,2,3)

process.T2RoadPadFinder.BiggestTubeAngleConsidered =cms.double(0.3)

process.T2RoadPadFinder.NumSigma= cms.double(6.)#Important for ALignment

#TolleranceDX

process.T2RoadPadFinder.NumPadCluOccupancyAlert= cms.double(50.)

process.T2RoadPadFinder.InefficiencyMaxJump= cms.int32(3)#2 is default

process.T2RoadPadFinder.Nmin_padsFinal= 4



process.load("RecoTotemT1T2.T2TrackProducer3.T2TrackColl3_cfi")

#process.load("TotemT1T2Validation.T2RecoValidation.T2RecoValidation2_cfi")



#from RecoTotemT1T2.T2TrackProducer3.T2TrackColl3_cfi import *

process.T2TrackColl3.StripFitting=cms.bool(False)

process.T2TrackColl3.RoadModuleLabel="T2RoadPadFinder"

process.T2TrackColl3.RoadInstanceLabel="NewRoadFinderRELOAD"

process.T2TrackColl3.verbosity=False

process.T2TrackColl3.RemoveOutliers=True #False for Internal  ALignment studies

process.T2TrackColl3.GhostSuppression=True

process.T2TrackColl3.PickUpDisplacedHit=True

process.T2TrackColl3.PickUpRadius=2.0















#process.g4SimHits.Generator.HepMCProductLabel = 'generator'

#process.g4SimHits.Physics.type = cms.string('SimG4Core/Physics/QGSP_BERT_EML')

#process.g4SimHits.Watchers = cms.VPSet()





#process.g4SimHits.Physics.DefaultCutValue = 0.001

#process.g4SimHits.Generator.LeaveScatteredProtons = cms.bool(False)



process.load("IOMC.EventVertexGenerators.VtxSmearedGauss_cfi")

process.VtxSmeared.MeanX = 0.0

process.VtxSmeared.MeanY = 0.0

process.VtxSmeared.MeanZ = 0.0



process.VtxSmeared.SigmaZ = 5.

                     

#process.p1 = cms.Path(process.generator*process.VtxSmeared*process.g4SimHits*process.mix)

process.p1 = cms.Path(process.generator*process.VtxSmeared*process.g4SimHits*process.mix*process.T2Digis*process.T2MCl*process.T2Hits*process.T2RoadPadFinder*process.T2TrackColl3)#*process.T2Digis





process.outpath = cms.EndPath(process.o1)













