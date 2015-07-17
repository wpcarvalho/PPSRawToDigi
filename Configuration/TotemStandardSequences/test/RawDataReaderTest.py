import FWCore.ParameterSet.Config as cms

process = cms.Process("RawDataDumper")
process.load("Configuration.TotemCommon.LoggerMin_cfi")

process.load('TotemRawData.Readers.RawDataSource_cfi')

process.source.fileNames.append('/castor/cern.ch/user/k/kmielnik/run_8318.000.vmea')
#process.source.fileNames.append('/castor/cern.ch/user/k/kmielnik/run_8318.001.vmea')

process.RawDataDumper = cms.EDAnalyzer("RawDataDumper")
process.RawDataDumper.RawEventLabel = cms.InputTag("source")

process.p = cms.Path(
    process.RawDataDumper
)