import FWCore.ParameterSet.Config as cms

process = cms.Process("prodT1T2RPDefault")

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2)
)

# Configure the output module (save the result in a file)
process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *', 'drop *_*mix*_*_*',  'drop *_*_*TrackerHits*_*', 'drop *_*_*Muon*_*', 'drop *_*_*Ecal*_*', 'drop *_*_*Hcal*_*', 'drop *_*_*Calo*_*', 'drop *_*_*Castor*_*', 'drop *_*_*FP420SI_*', 'drop *_*_*ZDCHITS_*', 'drop *_*_*BSCHits_*', 'drop *_*_*ChamberHits_*', 'drop *_*_*FibreHits_*', 'drop *_*_*WedgeHits_*'),
    fileName = cms.untracked.string('file:prodT1T2RPDefault.root')
)

# Configure if you want to detail or simple log information.
# LoggerMax -- detail log info output including: errors.log, warnings.log, infos.log, debugs.log
# LoggerMin -- simple log info output to the standard output (e.g. screen)
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# raw to digi conversion
process.load('TotemCondFormats/DAQInformation/DAQMappingSourceXML_cfi')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_220.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_147.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/t1_all_run1.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/t2_4quarters.xml')


################## STEP 0 - empty source

process.source = cms.Source("EmptySource")

################## STEP 1 - process.generator

# Use random number generator service
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# Pythia source labeled as "process.generator" @ CMSSW_3_1_1
process.load("Configuration.TotemCommon.PythiaSD_cfi")
process.generator.comEnergy = cms.double(2360)

################## STEP 2 process.SmearingGenerator

# declare optics parameters
#process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_1180GeV_11_cfi")

# Smearing
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")

################## STEP 3 process.g4SimHits

# Geometry
#process.load("Configuration.TotemCommon.geometryGlobal_cfi")
#process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_1180_Beta_11_220/RP_Dist_Beam_Cent.xml')

# Magnetic Field, by default we have 3.8T
process.load("Configuration.StandardSequences.MagneticField_cff")

# G4 simulation & proton transport
process.load("Configuration.TotemCommon.g4SimHits_cfi")
process.g4SimHits.Generator.HepMCProductLabel = 'generator'
#process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup
process.g4SimHits.Physics.DefaultCutValue = 100.


################## STEP 4 process.mix*process.T1Digis*process.t1cluster*process.t1rechit*process.t1roads*process.t1tracks

process.load("Configuration.TotemStandardSequences.T1_Digi_and_TrackReconstruction_cfi")

################## STEP 5 process.mix*process.T2Digis*process.T2MCl*process.T2Hits*process.T2RoadColl*process.T2TrackColl2

process.load("Configuration.TotemStandardSequences.T2_Digi_and_TrackReconstruction_cfi")

################## STEP 6 process.mix*process.RPSiDetDigitizer*process.RPClustProd*process.RPHecoHitProd*process.RPSinglTrackCandFind*process.RPSingleTrackCandCollFit

process.load("Configuration.TotemStandardSequences.RP_Digi_and_TrackReconstruction_cfi")

################## STEP 7 process.RP220Reconst

#process.load("RecoTotemRP.RPInelasticReconstruction.Rec_1180GeV_beta_11_220_cfi")
#process.RP220Reconst.BeamProtTransportSetup = process.BeamProtTransportSetup

################## STEP 8 process.RPCC*process.T2CC

process.load("L1TriggerTotem.CoincidenceChip.RPCoincidenceProducer_cfi")
process.load("L1TriggerTotem.CoincidenceChip.T2CoincidenceProducer_cfi")


process.p1 = cms.Path(process.generator*process.SmearingGenerator*process.g4SimHits*process.mix*process.T1Digis*process.t1cluster*process.t1rechit*process.t1roads*process.t1tracks2*process.T2Digis*process.T2MCl*process.T2Hits*process.T2RoadPadFinder*process.T2TrackColl3*process.RPSiDetDigitizer*process.RPClustProd*process.RPHecoHitProd*process.RPSinglTrackCandFind*process.RPSingleTrackCandCollFit*process.RPCC*process.T2CC)
process.outpath = cms.EndPath(process.o1)

