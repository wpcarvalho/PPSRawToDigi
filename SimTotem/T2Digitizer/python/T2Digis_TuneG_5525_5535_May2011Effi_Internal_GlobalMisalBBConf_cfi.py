import FWCore.ParameterSet.Config as cms


T2Digis = cms.EDProducer("T2DigiProducer",

   saveDigiVFAT=cms.bool(False),                      
   previousModule = cms.string('mix'),
   instanceLabel = cms.string('g4SimHitsTotemHitsT2Gem'),
   SimTrackContainerLabel = cms.InputTag("g4SimHits"),
   SimVertexContainerLabel = cms.InputTag("g4SimHits"),
   TakeOnlyPrim=cms.bool(False),
   TakeOnlySec=cms.bool(False),

   Misalignment = cms.PSet(                    
      simulatemisalign= cms.untracked.bool(True),   #reput True !!
      generaterandom= cms.untracked.bool(False),
      sigmaDispl= cms.untracked.double(0.20),
      sigmaPhiDispl= cms.untracked.double(0.02), # In radians
      sigmaGlobalShiftXY= cms.untracked.double(2.0),
      sigmaGlobalThetaXY= cms.untracked.double(0.003),
      verbosity= cms.untracked.bool(False),
      inputFileNameMisal = cms.untracked.string('SimTotem/T2Digitizer/data/run_5525-5535_IntAlignHX50000Evt_XovY0.3HIP_ANDGLOB_BBConf.dat')
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
        #3706_2by2_H01_23old
        inputFileNameEffi= cms.string('SimTotem/T2Digitizer/data/2011_RawVFATEffi_5525_5535.dat'),#,
        inputFileNameCorrupt= cms.string('SimTotem/T2Digitizer/data/run_5535_CorruptProb.txt'),#('SimTotem/T2Digitizer/data/run_5535_CorruptProb.txt')
        #SimTotem/T2Digitizer/data/Zero_DeadSector.dat
        inputFileNameDeadSect= cms.string('SimTotem/T2Digitizer/data/Zero_DeadSector.dat'),
        #SimTotem/T2Digitizer/data/DeadCH_3706_Apr2011_20000Evts_PL09TrkMult_4-6H1.dat
        inputFileNameDeadChannels= cms.string('SimTotem/T2Digitizer/data/DeadCH_5525-5535.dat'),#
        #SimTotem/T2Digitizer/data/NoisyCH_3706_Apr2011_20000Evts_PL09TrkMult_4-6H1.dat
        inputFileNameNoisyChannels= cms.string(''),
        #SimTotem/T2Digitizer/data/3706_Apr2011_20000Evts_PL09TrkMult_1-1H1Correl.root
        inputFileNameNoiseCorrelation= cms.string('')
        
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
        inputFileNameEffi= cms.string('SimTotem/T2Digitizer/data/2011_RawVFATEffi_5525_5535.dat'),#('SimTotem/T2Digitizer/data/2011_RawVFATEffi_5525_5535.dat'),
        inputFileNameCorrupt= cms.string('SimTotem/T2Digitizer/data/run_5535_CorruptProb.txt'),#SimTotem/T2Digitizer/data/run_5535_CorruptProb.txt
        inputFileNameDeadSect= cms.string('SimTotem/T2Digitizer/data/Zero_DeadSector.dat'),
        #SimTotem/T2Digitizer/data/DeadCH_3706_Apr2011_20000Evts_PL09TrkMult_4-6H1.dat
        inputFileNameDeadChannels= cms.string('SimTotem/T2Digitizer/data/DeadCH_5525-5535.dat'),#,
        #SimTotem/T2Digitizer/data/NoisyCH_3706_Apr2011_20000Evts_PL09TrkMult_4-6H1.dat
        inputFileNameNoisyChannels= cms.string(''),
        inputFileNameNoiseCorrelation= cms.string('')
    )
)



