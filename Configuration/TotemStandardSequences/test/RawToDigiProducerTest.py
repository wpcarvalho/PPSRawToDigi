import FWCore.ParameterSet.Config as cms

process = cms.Process("T1ChannelsHitRate")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.load("Configuration.TotemCommon.LoggerMin_cfi")

# raw data source
process.load("TotemRawData.Readers.RawDataSource_cfi")
process.source.fileNames.append('/castor/cern.ch/user/k/kmielnik/run_8318.000.vmea')

# mapping files
process.load('TotemCondFormats.DAQInformation.DAQMappingSourceXML_cfi')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_220.xml')
#process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_147.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/t1_all_run1.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/t2_4quarters.xml')

# mask files
process.DAQMappingSourceXML.maskFileNames.append('TotemCondFormats/DAQInformation/data/T1DeadNoisyChannelsListRun1.xml')

# raw to digi 
process.load('TotemRawData.RawToDigi.Raw2DigiProducer_cfi')
process.Raw2DigiProducer.verbosity = 2

process.o1 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('./recoT1_091101_mux.root'),
)

process.p1 = cms.Path(
    process.Raw2DigiProducer
)

process.outpath = cms.EndPath(process.o1)
