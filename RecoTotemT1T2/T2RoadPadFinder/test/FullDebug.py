import FWCore.ParameterSet.Config as cms

process = cms.Process("RecoAndAlign")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMax_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)  #-1
)
 



# input of raw data
process.source = cms.Source("RawDataSource",
        verbosity = cms.untracked.uint32(1),
        eventsToCheck = cms.uint32(6),
        skipCorruptedEvents = cms.bool(True),
        performChecks = cms.bool(True),
      # fileNames = cms.vstring('rfio:/castor/cern.ch/totem/LHCRawData/2010/T2/run_2811.000.vmea')
   #/castor/cern.ch/totem/LHCRawData/2010/T2/ION_RUNS/run_3791.000.vmea
        fileNames = cms.untracked.vstring('/tmp/berretti/run_3706.000.vmea'),#run_3686.023.vmea run_3791.000.vmea run_3706.000
        setRunNumberFromFileName = cms.bool(False),                      
       #fileNames = cms.vstring('rfio:/castor/cern.ch/user/b/berretti/DataAnalysis2010/RawData010/run_647.vmea')
         #fileNames = cms.vstring('/afs/cern.ch/exp/totem/scratch/data/T2/14June2010/run_647.vmea')
     #     fileNames = cms.vstring('/afs/cern.ch/exp/totem/scratch/berretti/tmp/BatchSubmitter/headrun_647.vmea')
      #  fileNames = cms.vstring('file:/afs/cern.ch/exp/totem/scratch/berretti/tmp/BatchSubmitter/run2003.vmea')
        #skipEvents = cms.untracked.uint32(0)                   
#rfcp /castor/cern.ch/user/b/berretti/DataAnalysis2010/RawData010/run_647.vmea /tmp/berretti/
) 


process.o1 = cms.OutputModule("PoolOutputModule",
   #   fileName = cms.untracked.string('rfio:/castor/cern.ch/user/b/berretti/DataAnalysis2010/RecoData010/TrackDevelop/run2003Dev1.root'),
        fileName = cms.untracked.string('file:/tmp/berretti/testRecoOLD.root'),
       	outputCommands = cms.untracked.vstring(	'drop TotemRawEvent_*_*_*',
						'keep T2DetIdT2DigiVfatTotemDigiCollection_*_*_*',
						'keep *_*_*T2VfatInformation*_*',
						'keep T1T2Tracks_*_*_*',
						'keep T2Hits_*_*_*',
						'keep T2Roads_*_*_*',
						'keep T2DetIdT2Clustersstdmap_*_*_*',
						'keep T2DetIdT2StripDigiTotemDigiCollection_*_*_*',
						'keep T2DetIdT2PadDigiTotemDigiCollection_*_*_*')
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
        verbosity = cms.untracked.uint32(0),
        discardHighOccupancyVfatverbosity= cms.untracked.bool(False)#IMPORTANT
)


process.load("RecoTotemT1T2.T2MakeCluster.T2MakeCluster_cfi")

process.T2MCl.TakeCleanEventOnly=cms.bool(False) #IMPORTANT
#process.T2MCl.BlobMinSize=cms.double(30.);

process.load("RecoTotemT1T2.T2RecHit.T2RecHit_cfi")

process.load("RecoTotemT1T2.T2RoadPadFinder.T2RoadPadFinder_cfi")
process.T2RoadColl.verbosity=False
process.T2RoadColl.TwoPointsTubesAngularCollinearity=0.07
process.T2RoadColl.AllowsPadReAssociation=True
process.T2RoadColl.Nmin_padsFinal = cms.int32(4)
process.T2RoadColl.MinimumNumCl1Hit= cms.int32(3)
process.T2RoadColl.chi2XYProb_Thr= cms.double(0.01)


process.load("RecoTotemT1T2.T2TrackProducer3.T2TrackColl3_cfi")
process.T2TrackColl3.StripFitting=cms.bool(False)


process.T2Hits.Cl1MaxPad = cms.uint32(25) #Tune better
process.T2Hits.Cl1MaxStrip = cms.uint32(25)
process.T2Hits.IncludeClass0Hits = True
#process.T2Hits.inputFileNameMisal=cms.untracked.string('/home/mirko/SL/WorkingArea/CMSSW_311SVN2/CMSSW_3_1_1/src/TotemAlignment/T2TrkBasedInternalAlignment/test/MisalOut/run_152_RecoXY_IntAlCfg1.dat')
#/afs/cern.ch/exp/totem/scratch/berretti/AnalysisInternalAlignment/run_648InternalAlignNOREF30000_NOGLOBAL.dat
process.T2Hits.inputFileNameMisal=cms.untracked.string('RecoTotemT1T2/T2RecHit/data/run_648InternalAlignNOREF30000GLOBAL_cf24_6_6_6.dat')
#run_648InternalAlignNOREF30000GLOBAL_cf24_6_6_6.dat
process.T2Hits.useTXTfile=cms.bool(True)
process.T2Hits.InsertAlignmentbyCFG=cms.bool(True)
process.T2Hits.verbosity=cms.untracked.bool(False)


process.p = cms.Path(process.RawToDigi*process.T2MCl*process.T2Hits*process.T2RoadColl*process.T2TrackColl3)

#process.p = cms.Path(process.RawToDigi*process.T2MCl*process.T2Hits*process.T2RoadColl)*process.T2TrackColl3
process.outpath = cms.EndPath(process.o1)
