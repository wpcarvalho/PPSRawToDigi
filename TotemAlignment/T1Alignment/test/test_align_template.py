import FWCore.ParameterSet.Config as cms

process = cms.Process("RealDataMonitorXML")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(3000)
)



# input of raw data
process.source = cms.Source("RawDataSource",
        verbosity = cms.untracked.uint32(1),                   
        eventsToCheck = cms.uint32(6),
        skipCorruptedEvents = cms.bool(False),
        performChecks = cms.bool(True),setRunNumberFromFileName=cms.bool(False),                       
        fileNames=cms.vstring(
             '/project/gruppo1/totem/Testbeam2010/november2010/vme_outfile_2127.vme2'         
)
)

process.o1 = cms.OutputModule("PoolOutputModule",
     fileName = cms.untracked.string('file:./analtracks_nov_yellow_mu_3500.root'),
	outputCommands = cms.untracked.vstring('drop TotemRawEvent_*_*_*','keep T1T2Tracks_*_*_*',
'keep T1DetIdT1RecHit2DsOwnedRangeMap_*__*',
                                               'keep T1RecHitGlobalss_*__*',
                                               'keep T1DetIdT1ClusterTotemDigiCollection_*__*',
                                               'keep T1DetIdT1DigiVfatTotemDigiCollection_*__*',
                                               'keep T1DetIdT1DigiWireTotemDigiCollection_*__*' ) )


process.DAQInformationSourceXML = cms.ESSource("DAQInformationSourceXML",
  xmlFileName = cms.string('/project/gruppo1/totem/work/sep2010/CMSSW_3_1_1/src/TotemDQM/Modules/test/T1_yellow_nov2010.xml')
)


process.RawToDigi = cms.EDProducer("T1XMLDataDigiProducer",
	verbosity = cms.untracked.uint32(0),
                                    GET_NOISY_MAP = cms.untracked.uint32(0),
                                    NOISY_MAP = cms.untracked.string("NOISY_MAP.txt") 
)



process.load("RecoTotemT1T2.T1MakeCluster.T1MakeCluster_cfi")

process.load("RecoTotemT1T2.T1RecHit.T1RecHit_cfi")

process.load("RecoTotemT1T2.T1RoadProducer.T1RoadProducer_cfi")


process.t1tracks = cms.EDProducer("T1TrackProducer",
    Verbosity = cms.int32(0),
    ChiRidThr = cms.double(40.0)
)

process.t1internalalign = cms.EDAnalyzer("T1InternalAlignment",
 Verbosity = cms.int32(0)
                                          )



process.p = cms.Path(process.RawToDigi*process.t1cluster*process.t1rechit*process.t1roads*process.t1tracks*process.t1internalalign)



process.outpath = cms.EndPath(process.o1)

