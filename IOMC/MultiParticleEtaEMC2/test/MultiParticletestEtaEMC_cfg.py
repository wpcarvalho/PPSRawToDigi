import FWCore.ParameterSet.Config as cms

process = cms.Process("T1T2gun")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.FrameworkJobReport.FwkJob.limit = 0

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    sourceSeed = cms.untracked.uint32(201456789),
    moduleSeeds = cms.PSet(
        VtxSmeared = cms.untracked.uint32(20340567),
        g4SimHits = cms.untracked.uint32(20340)
    )
)

process.load("SimGeneral.HepPDTESSource.pdt_cfi")


# Geometry
process.load("Configuration.TotemCommon.geometryT1T2_cfi")

################## STEP 1
process.source = cms.Source("EmptySource")

################## STEP 2 - process.generator
process.generator = cms.EDProducer("MultiParticleGunSourceEtaEMC2",
PGunParameters =cms.untracked.PSet(
     PartID = cms.untracked.vint32(-211,-211), 
	 MinPhis = cms.untracked.vdouble(-3.14159265),
	 MaxPhis = cms.untracked.vdouble(3.14159265),
	 InputMCFileName = cms.untracked.string("$CMSSW_BASE/src/IOMC/MultiParticleEtaEMC2/data/T2Pion.root"),
	 InputMCHistoName =  cms.untracked.string("PmenoEtaEnergy2")
      ),      
      Verbosity = cms.untracked.int32(0) # for printouts, set it to 1 (or greater) 
      #untracked bool  AddAntiParticle = false   # if you turn it ON, for each particle
                                                # an anti-particle will be generated,
                                                # with 3-mom opposite to the particle's   
   )

# Use random number generator service
process.load("Configuration.TotemCommon.RandomNumbers_cfi")


process.o1 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('./MultitestEtaEMC2.root')
)

process.p1 = cms.Path(process.generator)
process.outpath = cms.EndPath(process.o1)

   
