# Auto generated configuration file
# using: 
# Revision: 1.372.2.1 
# Source: /local/reps/CMSSW.admin/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/EightTeV/MinBias_Tune4C_8TeV_pythia8_cff.py --step GEN,SIM --beamspot Realistic8TeVCollision --conditions START52_V7::All --pileup NoPileUp --datamix NODATAMIXER --eventcontent RAWSIM --datatier GEN-SIM --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
#process.load("SimGeneral.MixingModule.mixLowLumPU_cfi")
#process.load("SimGeneral.MixingModule.mixProdStep2_cfi"). 
#process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(510)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('PYTHIA8-MinBias Tune 4C at 8TeV'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/CMSSW/Configuration/GenProduction/python/EightTeV/MinBias_Tune4C_8TeV_pythia8_cff.py,v $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('/tmp/berretti/testMu.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START44_V13::All'

#process.generator = cms.EDFilter("CosmicRayGeneratorFilter",
#    beammomentum = cms.double(4000),
#    targetmomentum = cms.double(-4000),
#    beamid = cms.int32(1),
#    targetid = cms.int32(1),
#    model = cms.int32(7),
#    paramFileName = cms.string("/afs/cern.ch/user/b/berretti/public/cms.param")
#)
process.generator = cms.EDProducer("FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(13),            # 13 (particle id) = mu^-
        MinEta = cms.double(-6.4),
        MaxEta = cms.double(-5.35),
        MinPhi = cms.double(-3.14159265358979323846),  # in radians
        MaxPhi = cms.double(3.14159265358979323846),
        MinE = cms.double(49.0),
        MaxE = cms.double(50.01)
    ),
    AddAntiParticle = cms.bool(False),
    Verbosity = cms.untracked.int32(1)  # set to 1 (or greater) for printouts
)

#process.generator = cms.EDFilter("ReggeGribovPartonMCGeneratorFilter",
#                    beammomentum = cms.double(4000),
#                    targetmomentum = cms.double(-4000),
#                    beamid = cms.int32(1),
#                    targetid = cms.int32(1),
#                    model = cms.int32(0),#7 qsjetII-4, 0 epos, 6 syb, 2 QGSJET01 ReggeGribovPartonMC.param				 
#                    paramFileName = cms.string("Configuration/Generator/data/ReggeGribovPartonMC_8TeV.param")
#                   )
#ReggeGribovPartonMC_8TeV.param works for epos retune

#process.generator = cms.EDProducer("FlatRandomEGunProducer",
#    PGunParameters = cms.PSet(
#        PartID = cms.vint32(-211),            # 13 (particle id) = mu^- ||||||  -211=pi-
#        MinEta = cms.double(3.5),
#        MaxEta = cms.double(4.7),
#        MinPhi = cms.double(-3.14159265358979323846),  # in radians
#        MaxPhi = cms.double(3.14159265358979323846),
#        MinE = cms.double(1.0),
#        MaxE = cms.double(50.01)
#    ),
#    AddAntiParticle = cms.bool(False),
#    Verbosity = cms.untracked.int32(0)  # set to 1 (or greater) for printouts
#)


#process.generator = cms.EDFilter("Pythia8GeneratorFilter",
#    pythiaPylistVerbosity = cms.untracked.int32(1),
#    filterEfficiency = cms.untracked.double(1.0),
#    pythiaHepMCVerbosity = cms.untracked.bool(False),
#    comEnergy = cms.double(8000.0),
#    crossSection = cms.untracked.double(72850000000),
#    maxEventsToPrint = cms.untracked.int32(0),
#    PythiaParameters = cms.PSet(
#        processParameters = cms.vstring('Main:timesAllowErrors    = 10000',
#            'ParticleDecays:limitTau0 = on',
#            'ParticleDecays:tauMax = 10',
#            'SoftQCD:minBias = on',
#            'SoftQCD:singleDiffractive = on',
#            'SoftQCD:doubleDiffractive = on',
#            'Tune:pp 5',
#            'Tune:ee 3'),
#        parameterSets = cms.vstring('processParameters')
#    )
#)



process.g4SimHits.Physics.DefaultCutValue = 0.001
process.g4SimHits.Generator.LeaveScatteredProtons = cms.bool(False)



#process.g4SimHits.HCalSD.UseShowerLibrary = cms.bool(True)
#process.g4SimHits.HCalSD.UseParametrize = cms.bool(False)
#process.g4SimHits.HCalSD.UsePMTHits = cms.bool(False)
#process.g4SimHits.HCalSD.UseFibreBundleHits = cms.bool(True)
#process.g4SimHits.HFShower.UseShowerLibrary = cms.bool(True)
#process.g4SimHits.HFShower.UseHFGflash = cms.bool(False)
#process.g4SimHits.HFShower.ApplyFiducialCut = cms.bool(False)
#process.g4SimHits.CastorSD.nonCompensationFactor = cms.double(0.85)

process.g4SimHits.HCalSD.UseShowerLibrary = cms.bool(False)
process.g4SimHits.HCalSD.UseParametrize = cms.bool(False)#Should be always false
process.g4SimHits.HCalSD.UsePMTHits = cms.bool(False)
process.g4SimHits.HCalSD.UseFibreBundleHits = cms.bool(True)
process.g4SimHits.HFShower.UseShowerLibrary = cms.bool(False)
process.g4SimHits.HFShower.UseHFGflash = cms.bool(True)
process.g4SimHits.HFShower.ApplyFiducialCut = cms.bool(False)
process.g4SimHits.CastorSD.nonCompensationFactor = cms.double(0.85)






process.g4SimHits.Generator.ApplyPCuts = cms.bool(False)
process.g4SimHits.Generator.ApplyEtaCuts = cms.bool(False)
process.g4SimHits.Generator.MinPCut = cms.double(0.001)
process.g4SimHits.Generator.ApplyPCuts = cms.bool(False)
process.g4SimHits.RDecLenCut = cms.double(0.2)
process.g4SimHits.EtaCutForHector = cms.double(20)

process.ProductionFilterSequence = cms.Sequence(process.generator)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)

process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)







process.load("SimTotem.T2Digitizer.T2DigisSTDGeoEffi_cfi")


#process.load("SimTotem.T2Digitizer.T2DigisOnlyPrim_cfi")
process.T2Digis.saveDigiVFAT=cms.bool(True)
process.T2Digis.Misalignment.simulatemisalign=cms.untracked.bool(True)#UseAlign
#Decomment for alignment multi-scenario validation
process.T2Digis.Misalignment.inputFileNameMisal="SimTotem/T2Digitizer/data/2012IntGlobAlignV2.dat"#DigiMisalFilename
  







process.load("RecoTotemT1T2.T2MakeCluster.T2MakeCluster_cfi")



process.T2MCl.TakeCleanEventOnly=cms.bool(False) #IMPORTANT
#process.T2MCl.maskvect = cms.vuint32(21,23,25,27,29)
process.T2MCl.SimuClusterEfficiency=cms.bool(True)#SimulateClusterEffi
process.T2MCl.EffiGeoRootFileName= cms.string("RecoTotemT1T2/T2MakeCluster/data/Geom_effiOutput_All8372_Effi_V2Pl.root")

process.load("RecoTotemT1T2.T2RecHit.T2RecHit_cfi")



#from RecoTotemT1T2.T2RecHit.T2RecHit_cfi import *
process.T2Hits.Cl1MaxPad = cms.uint32(25) #Tune better
process.T2Hits.Cl1MaxStrip = cms.uint32(25)
process.T2Hits.IncludeClass0Hits = True
#T2Hits.inputFileNameMisal=cms.untracked.string('/afs/cern.ch/exp/totem/scratch/berretti/tmp/testSplitMerge/CMSSW_3_1_1/src/$
#process.T2Hits.inputFileNameMisal=cms.untracked.string('SimTotem/T2Digitizer/data/run_5525-5535_IntAlignHX50000Evt_XovY0.3HIP_ANDGLOB_BBConf.dat')
process.T2Hits.inputFileNameMisal=cms.untracked.string('SimTotem/T2Digitizer/data/2012IntGlobAlignV2.dat')#DigiFNameForHit

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
process.T2RoadPadFinder.NumSigma= cms.double(2.)#Important for ALignment
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
process.T2TrackColl3.verbosity=True
process.T2TrackColl3.RemoveOutliers=True #False for Internal  ALignment studies
process.T2TrackColl3.GhostSuppression=True
process.T2TrackColl3.PickUpDisplacedHit=False
process.T2TrackColl3.PickUpRadius=2.0







#process.g4SimHits.Generator.HepMCProductLabel = 'generator'
#process.g4SimHits.Physics.type = cms.string('SimG4Core/Physics/QGSP_BERT_EML')
#process.g4SimHits.Watchers = cms.VPSet()


process.g4SimHits.Physics.DefaultCutValue = 0.001
process.g4SimHits.Generator.LeaveScatteredProtons = cms.bool(False)
process.g4SimHits.TotemSD.Verbosity=11

process.load("IOMC.EventVertexGenerators.VtxSmearedGauss_cfi")

#process.VtxSmeared.MeanX = 0.0
#process.VtxSmeared.MeanY = 0.0
process.VtxSmeared.MeanZ = 0.0
process.VtxSmeared.MeanX = 0.0
process.VtxSmeared.MeanY = 0.0
#process.VtxSmeared.MeanZ = 1125.
process.VtxSmeared.SigmaX = 0.01
process.VtxSmeared.SigmaY = 0.01
process.VtxSmeared.SigmaZ = 5.



process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    moduleSeeds = cms.PSet(
        generator = cms.untracked.uint32(22223),#Seed1
        g4SimHits = cms.untracked.uint32(22223),
        VtxSmeared = cms.untracked.uint32(22223),
        T2Digis = cms.untracked.uint32(22223),
        T2MCl = cms.untracked.uint32(22223),
    ),
    sourceSeed = cms.untracked.uint32(22223),
   
)



process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/ForwardSimData/data/totemsensGem.xml")

process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("SimGeneral.MixingModule.cfwriter_cfi")
process.load("TotemT1T2Validation.T2BackgroundAnalysis.T2BackgroundAn_cfi")

process.T2BackgroundAn.TrackLabel=cms.string('T2TrackColl3')#T2TrackColl2XYSecondary T2TrackColl3 T2TrackColl2
process.T2BackgroundAn.VtxClusterDistance=cms.double(5.)
process.T2BackgroundAn.selected_event=2
process.T2BackgroundAn.verbosity=2
#Full:(0,1,2,3) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #0, for DN/Deta Analysis             
process.T2BackgroundAn.T2_QuarterUsed = cms.vint32(0,1,2,3)       #0,1,2,3 for Vtx efficiency
#Full:(0,1,2,3) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
process.T2BackgroundAn.fastSimulation=cms.bool(True)
process.T2BackgroundAn.fastAnalysis=cms.bool(True)

process.T2BackgroundAn.PadRoadFinderAnalysis=cms.bool(True)
process.T2BackgroundAn.RoadInstanceLabel="NewRoadFinderRELOAD"
process.T2BackgroundAn.RoadLabel= cms.string('T2RoadPadFinder') #needed only for pad-road analysis T2RoadPadFinder T2RoadColl2
process.T2BackgroundAn.outputFileName=cms.string('Mueffi.root')
process.T2BackgroundAn.numhitRequiredFormMatching = 2
process.T2BackgroundAn.vtxseednumber = 56659

process.T2BackgroundAn.PtCutinPrimaryEfficiency= 0.04 # ex 0.03 Gev/c
process.T2BackgroundAn.EnergyCutinPrimaryEfficiency= 5.0
process.T2BackgroundAn.ZEffiCutImpact= cms.vdouble(6600, 6600, 6300, 6150)#(4400.0,4400.0,4200.0,4100.0)#From analysis, this is the simu 95% area
#Here I feed the cmssw process
process.T2BackgroundAn.MaxPadCluOfFittingCurves = 65
process.T2BackgroundAn.NameOfGenerator = cms.string('Py8')


#process.simulation_step = cms.Path(process.psim*process.mix)
process.simulation_step = cms.Path(process.psim*process.mix*process.T2Digis*process.T2MCl*process.T2Hits*process.T2RoadPadFinder*process.T2TrackColl3*process.T2BackgroundAn)

process.RAWSIMEventContent.outputCommands.append("keep *_T2Hits_*_*")
process.RAWSIMEventContent.outputCommands.append("keep *_T2TrackColl3_*_*")
process.RAWSIMEventContent.outputCommands.append("keep *_T2MCl_*_*")
process.RAWSIMEventContent.outputCommands.append("keep *_T2Digis_*_*")
#process.RECOSIMoutput.outputCommands.append("keep *_CFWriter_*_*")
process.RAWSIMEventContent.outputCommands.append("keep *_mix_*_*")

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 


#print process.dumpConfig()


#*process.T2Digis*process.T2MCl*process.T2Hits*process.T2RoadPadFinder*process.T2TrackColl3
