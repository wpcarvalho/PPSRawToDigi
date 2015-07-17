import FWCore.ParameterSet.Config as cms

process = cms.Process("RealDataMonitorXML")

# minimum of logs
#process.load("Configuration.TotemCommon.LoggerMax_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20000)
)



# input of raw data
process.source = cms.Source("RawDataSource",
         verbosity = cms.untracked.uint32(1),                   
	eventsToCheck = cms.uint32(6),
        skipCorruptedEvents = cms.bool(False),
        performChecks = cms.bool(False),                       
        fileNames = cms.vstring('TotemRawData/RawToDigi/python/run_151.vmea')
)


process.o1 = cms.OutputModule("PoolOutputModule",
     fileName = cms.untracked.string('file:TotemRawData/RawToDigi/python/run_151_BeforeHit20000.root'),
	outputCommands = cms.untracked.vstring('drop TotemRawEvent_*_*_*','keep T2DetIdT2DigiVfatTotemDigiCollection_*_*_*','keep *_*_*T2VfatInformation*_*','keep T1T2Tracks_*_*_*','keep T2Hits_*_*_*','keep T2Roads_*_*_*','keep T2DetIdT2Clustersstdmap_*_*_*','keep T2DetIdT2StripDigiTotemDigiCollection_*_*_*','keep T2DetIdT2PadDigiTotemDigiCollection_*_*_*')
    )


#process.dump = cms.EDAnalyzer("EventContentAnalyzer")

#T2GeoMapIP5_4quarter_cmssw.xml
#T2MapIP5_D3_D2_NewFormatBIS.xml
#Load T2 vfat mapping
process.DAQInformationSourceXML = cms.ESSource("DAQInformationSourceXML",
  xmlFileName = cms.string('TotemRawData/RawToDigi/python/T2GeoMapIP5_4quarter_vmea_Full.xml')
  #  xmlFileName = cms.string('TotemRawData/RawToDigi/python/T2GeoMapIP5_4quarter_vmea_cmssw.xml')
)

#  outputCommands = cms.untracked.vstring('drop TotemRawEvent_*_*_*')
# Fill T2 digi and vfat object
process.RawToDigi = cms.EDProducer("T2XMLDataDigiProducer",
	verbosity = cms.untracked.uint32(10)
)

	
#process.T2Hits.Cl1MaxPad =20
#process.T2Hits.Cl1MaxStrip =20
#process.T2Hits.IncludeClass0Hits = True


process.load("RecoTotemT1T2.T2MakeCluster.T2MakeCluster_cfi")



#process.p = cms.Path(process.RawToDigi)
process.p = cms.Path(process.RawToDigi*process.T2MCl)
#*process.dqm1 *process.T2Hits*process.T2RoadColl*process.T2TrackColl2
process.outpath = cms.EndPath(process.o1)
