import FWCore.ParameterSet.Config as cms


T2Digis = cms.EDProducer("T2DigiProducer",
   previousModule = cms.string('mix'),
   instanceLabel = cms.string('g4SimHitsTotemHitsT2Gem'),

   saveDigiVFAT=cms.bool(False),                      
   SimTrackContainerLabel = cms.InputTag("g4SimHits"),
   SimVertexContainerLabel = cms.InputTag("g4SimHits"),
   TakeOnlyPrim=cms.bool(False),
   TakeOnlySec=cms.bool(False),
   Misalignment = cms.PSet(                    
      simulatemisalign= cms.untracked.bool(False),
      generaterandom= cms.untracked.bool(False),
      sigmaDispl= cms.untracked.double(0.20),
      sigmaPhiDispl= cms.untracked.double(0.02), # In radians
      sigmaGlobalShiftXY= cms.untracked.double(2.0),
      sigmaGlobalThetaXY= cms.untracked.double(0.003),
      verbosity= cms.untracked.bool(False),
      inputFileNameMisal = cms.untracked.string('SimTotem/T2Digitizer/test/testMisal.txt')
      ),


                         

    stripHit = cms.PSet(
       # diffCoeff = cms.double(0.26),           # diffusion coefficient of driftzone
        NUMBER_OF_SIM_STEPS = cms.int32(30),    # number of uniformly distributed electron-hole pairs
        eIonE = cms.double(28.0),               # Energy to create electron/ionpair
        #gain = cms.double(8000.0),              # Gain of GEM detector
        z_max = cms.double(0.9),                # dimensions of driftzone in cm
        z_min = cms.double(0.6),
        StripWidth = cms.vdouble(0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036,0.036),
	diffCoeff = cms.vdouble(0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34),     
	gain =  cms.vdouble(30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000,30000)
        
    ),

    padHit = cms.PSet(
        #diffCoeff = cms.double(0.26),           # diffusion coefficient of driftzone
        NUMBER_OF_SIM_STEPS = cms.int32(30),    # number of uniformly distributed electron-hole pairs
        eIonE = cms.double(28.0),               # Energy to create electron/ionpair
        #gain = cms.double(8000.0),              # Gain of GEM detector
        z_max = cms.double(0.9),                # dimensions of driftzone in cm
        z_min = cms.double(0.6),
        diffCoeff = cms.vdouble(0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34,0.34),     
	gain =  cms.vdouble(15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000,15000)
        
    ),

                         

    stripElec = cms.PSet(
        bins = cms.vint32(30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30),
	capaNoiseFactorStrip = cms.vdouble(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1),
	sigmaExtraNoise = cms.vint32 (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	simpleThreshold = cms.vint32 (400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400),
        UseDBInfo=cms.bool(False),
        UseCFGInfo=cms.bool(True),
        UseVFATs=cms.bool(True),
        SetVfatEfficiency= cms.bool(True),
        inputFileNameEffi= cms.string('SimTotem/T2Digitizer/data/VfatEffiH0123_4380.dat'),
        inputFileNameCorrupt= cms.string('SimTotem/T2Digitizer/data/Zero_Corruption.dat'),
        inputFileNameDeadSect= cms.string('SimTotem/T2Digitizer/data/3706_DeadSector.dat'),
        inputFileNameDeadChannels= cms.string('SimTotem/T2Digitizer/data/Zero_DeadChannel.dat')                    
    ),

   
    
    padElec = cms.PSet(

    bins = cms.vint32(30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30),
	capaNoiseFactorPad = cms.vdouble(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1),
	sigmaExtraNoise = cms.vint32 (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	simpleThreshold = cms.vint32 (400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400),
        UseDBInfo=cms.bool(False), # osolete, this flag is no longer supported, because of removed TotemDatabaseService module
        UseCFGInfo=cms.bool(True),
        UseVFATs=cms.bool(True),
        SetVfatEfficiency= cms.bool(True),
        inputFileNameEffi= cms.string('SimTotem/T2Digitizer/data/VfatEffiH0123_4380.dat'),
        inputFileNameCorrupt= cms.string('SimTotem/T2Digitizer/data/Zero_Corruption.dat'),
        inputFileNameDeadSect= cms.string('SimTotem/T2Digitizer/data/3706_DeadSector.dat'),
        inputFileNameDeadChannels= cms.string('SimTotem/T2Digitizer/data/Zero_DeadChannel.dat')                 
    )
)


