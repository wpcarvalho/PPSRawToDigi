import FWCore.ParameterSet.Config as cms

process = cms.Process("T1Check")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.load("Configuration.TotemCommon.LoggerMax_cfi")

# raw data source
process.load("TotemRawData.Readers.RawDataSource_cfi")
#process.source.fileNames.append('/afs/cern.ch/user/f/fferro/WORKSPACE/CMSSW_6_2_0/src/Configuration/TotemStandardSequences/test/run_6917.000.vmea')

process.source.fileNames.append('/castor/cern.ch/totem/LHCRawData/2011/Physics/run_6917.000.vmea')

# mapping files
process.load('TotemCondFormats.DAQInformation.DAQMappingSourceXML_cfi')
#process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_220.xml')
#process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_147.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/t1_all_run1.xml')
#process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/t2_4quarters.xml')

# mask files
process.DAQMappingSourceXML.maskFileNames.append('TotemCondFormats/DAQInformation/data/T1DeadNoisyChannelsListRun1.xml')

# raw to digi 
process.load('TotemRawData.RawToDigi.Raw2DigiProducer_cfi')
process.Raw2DigiProducer.verbosity = 0


process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *', 'drop *_*_*TrackerHits*_*', 'drop *_*_*Muon*_*', 'drop *_*_*Ecal*_*', 'drop *_*_*Hcal*_*', 'drop *_*_*Calo*_*', 'drop *_*_*Castor*_*', 'drop *_*_*FP420SI_*', 'drop *_*_*ZDCHITS_*', 'drop *_*_*BSCHits_*', 'drop *_*_*ChamberHits_*', 'drop *_*_*FibreHits_*', 'drop *_*_*WedgeHits_*'),
    fileName = cms.untracked.string('file:T1_6917.root')
)

# Configure if you want to detail or simple log information.
# LoggerMax -- detail log info output including: errors.log, warnings.log, infos.log, debugs.log
# LoggerMin -- simple log info output to the standard output (e.g. screen)
process.load("Configuration.TotemCommon.LoggerMax_cfi")

# Use particle table
process.load("SimGeneral.HepPDTESSource.pdt_cfi")




# Use random number generator service
process.load("Configuration.TotemCommon.RandomNumbers_cfi")


process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_1180GeV_11_cfi")



# Geometry
process.load("Configuration.TotemCommon.geometryT1T2_cfi")

# Magnetic Field, by default we have 3.8T
process.load("Configuration.StandardSequences.MagneticField_cff")



process.load("RecoTotemT1T2.T1MakeCluster.T1MakeCluster_cfi")

process.load("RecoTotemT1T2.T1RecHit.T1RecHit_cfi")

process.load("RecoTotemT1T2.T1RoadProducer.T1RoadProducer_cfi")

process.load("RecoTotemT1T2.T1TrackProducer2.T1TrackProducer2_cfi")





process.t1cluster.T1DigiVfatCollectionLabel = cms.InputTag("Raw2DigiProducer", "t1DataOutput")
process.t1cluster.ActivateDeadChannels = cms.bool(True)

process.t1rechit.T1ClusterCollectionLabel = cms.InputTag("t1cluster")
process.t1rechit.T1DigiWireCollectionLabel = cms.InputTag("Raw2DigiProducer", "t1DataOutput")
process.t1rechit.T1RecHit2DCollectionLabel = cms.InputTag("t1rechit")

process.t1roads.T1RecHit2DCollectionLabel = cms.InputTag("t1rechit")
process.t1roads.Alignment = cms.bool(True)

process.t1tracks2.T1RoadCollectionLabel = cms.InputTag("t1roads")




process.load("TotemT1T2Validation.T1Validation.T1Validation_cfi")
process.t1valid.OutputFile = cms.string("file:valT1_6917.root")
process.t1valid.SIM = cms.double(0)
process.t1valid.TrackLabel = cms.string('t1tracks2')
process.t1valid.T1DigiWireCollectionLabel = cms.InputTag("Raw2DigiProducer", "t1DataOutput")
process.t1valid.T1DigiVfatCollectionLabel = cms.InputTag("Raw2DigiProducer", "t1DataOutput")

process.p1 = cms.Path(

                       process.Raw2DigiProducer
                      *process.t1cluster
                      *process.t1rechit
                      *process.t1roads
                      *process.t1tracks2
                       *process.t1valid
                      )

process.outpath = cms.EndPath(process.o1)




