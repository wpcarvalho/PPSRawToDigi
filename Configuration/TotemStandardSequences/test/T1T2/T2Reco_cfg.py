import FWCore.ParameterSet.Config as cms

process = cms.Process("T2Reco")

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)

# Configure the output module (save the result in a file)
process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *', 'drop *_*_*TrackerHits*_*', 'drop *_*_*Muon*_*', 'drop *_*_*Ecal*_*', 'drop *_*_*Hcal*_*', 'drop *_*_*Calo*_*', 'drop *_*_*Castor*_*', 'drop *_*_*FP420SI_*', 'drop *_*_*ZDCHITS_*', 'drop *_*_*BSCHits_*', 'drop *_*_*ChamberHits_*', 'drop *_*_*FibreHits_*', 'drop *_*_*WedgeHits_*'),
    fileName = cms.untracked.string('file:T2Reco.root')
)

# Configure if you want to detail or simple log information.
# LoggerMax -- detail log info output including: errors.log, warnings.log, infos.log, debugs.log
# LoggerMin -- simple log info output to the standard output (e.g. screen)
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# Use particle table
process.load("SimGeneral.HepPDTESSource.pdt_cfi")


# raw to digi conversion
process.load('TotemCondFormats/DAQInformation/DAQMappingSourceXML_cfi')
#process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_220.xml')
#process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_147.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/t1_all_run1.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/t2_4quarters.xml')


################## STEP 0 - empty source

process.source = cms.Source("EmptySource")

################## STEP 1 - process.generator

# Use random number generator service
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# Pythia source labeled as "process.generator" @ CMSSW_4_2_4
#process.load("Configuration.TotemCommon.PythiaMB_cfi")
#process.generator.comEnergy = cms.double(2360)

import IOMC.Phojet.Phojet_cfi
process.generator = IOMC.Phojet.Phojet_cfi.generator

################## STEP 2 process.SmearingGenerator

# declare optics parameters
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_1180GeV_11_cfi")

# Smearing
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")

################## STEP 3 process.g4SimHits

# Geometry
process.load("Configuration.TotemCommon.geometryT1T2_cfi")

# Magnetic Field, by default we have 3.8T
process.load("Configuration.StandardSequences.MagneticField_cff")

# G4 simulation & proton transport
process.load("Configuration.TotemCommon.g4SimHits_cfi")
process.g4SimHits.Generator.HepMCProductLabel = 'generator'
process.g4SimHits.Physics.type = cms.string('SimG4Core/Physics/QGSP_BERT_EML')
process.g4SimHits.Watchers = cms.VPSet()
process.g4SimHits.Generator.LeaveScatteredProtons = cms.bool(False)

################## STEP 4 process.mix*process.T1Digis*process.t1cluster*process.t1rechit*process.t1roads*process.t1tracks

# # No pile up for the mixing module
# process.load("SimGeneral.MixingModule.mixNoPU_cfi")

# No pile up for the mixing module
process.load("Configuration.TotemCommon.mixNoPU_cfi")

########################### DIGI + RECO T1 ##########################################



################## STEP 5 process.mix*process.T2Digis*process.T2MCl*process.T2Hits*process.T2RoadColl*process.T2TrackColl2

# No pile up for the mixing module

########################### DIGI + RECO T2 ##########################################

#process.load("SimTotem.T2Digitizer.T2Digis_TuneG_5525_5535_May2011Effi_Internal_GlobalMisalBBConf_cfi")
process.load("SimTotem.T2Digitizer.T2Digis_cfi")

#2011-Tune Digitized data.
#use T2Digis for generic purpose simulation

process.load("RecoTotemT1T2.T2MakeCluster.T2MakeCluster_cfi")
process.T2MCl.maskvect = cms.vuint32()  #no masking by default
#Use T2MCl.maskvect = cms.vuint32(21,23,25,27,29) For the 2011 (After June) Run

process.load("RecoTotemT1T2.T2RecHit.T2RecHit_cfi")
process.T2Hits.inputFileNameMisal=cms.untracked.string('SimTotem/T2Digitizer/data/run_5525-5535_IntAlignHX50000Evt_XovY0.3HIP_ANDGLOB_BBConf.dat')
# put useTXTfile and InsertAlignmentbyCFG to False in order to avoid Hit alignment correction.
process.T2Hits.useTXTfile=cms.bool(False)   
process.T2Hits.InsertAlignmentbyCFG=cms.bool(False)  
process.T2Hits.verbosity=cms.untracked.bool(True)
#T2Hits.CorrectWithResolution=cms.bool(True) #False:Old Strategy

process.load("RecoTotemT1T2.T2RoadPadFinder.NewLabelT2RoadPadFinder_cfi")

process.T2RoadPadFinder.T2RoadCollProdName="NewRoadFinderRELOAD"
process.T2RoadPadFinder.useStraightPadTowers= cms.bool(True)#False: Old Strategy
process.T2RoadPadFinder.ResolveOverlapDoubleCount = cms.bool(False) # False for shadow alignment and dndeta An
process.T2RoadPadFinder.BiggestTubeAngleConsidered =cms.double(0.3)


#process.load("RecoTotemT1T2.T2TrackProducer3.T2TrackColl3_cfi")
#T2TrackColl3.RoadModuleLabel="T2RoadPadFinder"
#T2TrackColl3.RoadInstanceLabel="NewRoadFinderRELOAD"
#T2TrackColl3.RemoveOutliers=True
#T2TrackColl3.GhostSuppression=True


################## STEP 6 process.T2CC

#process.load("L1TriggerTotem.CoincidenceChip.T2CoincidenceProducer_cfi")


process.T2Digis.TakeOnlyPrim = cms.bool(False)
process.T2Digis.TakeOnlySec = cms.bool(False)

process.T2MCl.T2PadDigiCollectionLabel = cms.InputTag("T2Digis", "T2PadDigi")
process.T2MCl.T2StripDigiCollectionLabel = cms.InputTag("T2Digis", "T2StripDigi")

#process.T2Hits.LabelproductInstanceName=cms.string("T2Hits")
#process.T2Hits.ModuleLabelInput=cms.string("T2MCl")

process.p1 = cms.Path(process.generator
                      *process.SmearingGenerator
                      *process.g4SimHits
                      *process.mix
                      *process.T2Digis
                      *process.T2MCl
                      *process.T2Hits
                      *process.T2RoadPadFinder
#                      *process.T2TrackColl3
#                      *process.T2CC
                      )

process.outpath = cms.EndPath(process.o1)
