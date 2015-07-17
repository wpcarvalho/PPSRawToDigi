import FWCore.ParameterSet.Config as cms

process = cms.Process("RealDataMonitorXML")

# minimum of logs
#process.load("Configuration.TotemCommon.LoggerMax_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)



# input of raw data
process.source = cms.Source("RawDataSource",
         verbosity = cms.untracked.uint32(1),                   
	eventsToCheck = cms.uint32(6),
        skipCorruptedEvents = cms.bool(False),
        performChecks = cms.bool(False),                       
        fileNames = cms.vstring('TotemRawData/RawToDigi/python/calpulseScan_Threshold100_channelsOneverywhere_128events.vme2')
)


process.o1 = cms.OutputModule("PoolOutputModule",
     fileName = cms.untracked.string('file:/home/mirko/SL/WorkingArea/CMSSW_311SVN2/CMSSW_3_1_1/src/TotemRawData/RawToDigi/python/calpulseScan_Threshold100_channelsOneverywhere_128events.root'),
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

process.load("RecoTotemT1T2.T2RecHit.T2RecHit_cfi")

process.load("RecoTotemT1T2.T2RoadProducer.T2RoadCollDr15Dz200Dphi17_cfi")

process.load("RecoTotemT1T2.T2TrackProducer2.T2TrackColl25H_cfi")


process.T2Hits.Cl1MaxPad =20
process.T2Hits.Cl1MaxStrip =20
process.T2Hits.IncludeClass0Hits = True
 #   process.T2TrackColl2.verbosity = cms.bool(False),







process.load("RecoTotemT1T2.T2MakeCluster.T2MakeCluster_cfi")

process.load("RecoTotemT1T2.T2RecHit.T2RecHit_cfi")

#process.load("RecoTotemT1T2.T2RoadProducer.T2RoadCollDr15Dz200Dphi17_cfi")

process.load("RecoTotemT1T2.T2TrackProducer2.T2TrackColl25H_cfi")


	


process.load("RecoTotemT1T2.T2RoadProducer.T2RoadCollDr15Dz200Dphi17_cfi")



#process.load("RecoTotemT1T2.T2TrackProducer2.T2TrackColl2_cfi")


process.T2Hits.Cl1MaxPad = cms.uint32(20)
process.T2Hits.Cl1MaxStrip = cms.uint32(20)
process.T2Hits.IncludeClass0Hits = False
process.T2Hits.inputFileNameMisal=cms.untracked.string('/home/mirko/SL/WorkingArea/CMSSW_311SVN2/CMSSW_3_1_1/src/TotemAlignment/T2TrkBasedInternalAlignment/test/MisalOut/run_152_RecoXY_IntAlCfg1.dat')
process.T2Hits.useTXTfile=cms.bool(False)
process.T2Hits.InsertAlignmentbyCFG=cms.bool(False)
process.T2Hits.verbosity=cms.untracked.bool(False)



process.T2RoadColl.DeltaR=1.5
process.T2RoadColl.DeltaPhi=2.0
process.T2RoadColl.DeltaZ = 300.0

process.T2RoadColl.MinHits = 5  
#process.T2RoadColl.MaxHits = 10

process.T2TrackColl2.forceRZfit=False
process.T2TrackColl2.forceXYfit=True
process.T2TrackColl2.MinHitFinalTrk=5
process.T2TrackColl2.MaxHitsInRoad=15
process.T2TrackColl2.verbosity=False
# DQM modules
#HT 0-1 Arm0 |---| HT 2-3 Arm1
#RecoWithAlign/Reco_InternalAlCorr_MergedHalf0TrackFindingXY_2.py






process.p = cms.Path(process.RawToDigi*process.T2MCl*process.T2Hits*process.T2RoadColl*process.T2TrackColl2) 
#*process.dqm1
process.outpath = cms.EndPath(process.o1)
