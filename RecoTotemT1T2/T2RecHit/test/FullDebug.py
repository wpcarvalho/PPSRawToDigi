import FWCore.ParameterSet.Config as cms

process = cms.Process("RecoAndAlign")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMax_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)  #-1
)




# input of raw data
process.source = cms.Source("RawDataSource",
        verbosity = cms.untracked.uint32(1),
        eventsToCheck = cms.uint32(6),
        skipCorruptedEvents = cms.bool(True),
        performChecks = cms.bool(True),      
        fileNames = cms.untracked.vstring('/afs/cern.ch/exp/totem/scratch/data/T2/14June2010/run_647.vmea'),#run_3686.023.vmea
        setRunNumberFromFileName = cms.bool(False),                        
)


process.o1 = cms.OutputModule("PoolOutputModule",
       fileName = cms.untracked.string('file:testT2Reco.root'),
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

process.load("RecoTotemT1T2.T2RoadProducer.T2RoadCollDr15Dz200Dphi17_cfi")

process.load("RecoTotemT1T2.T2TrackProducer2.T2TrackColl25H_cfi")

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


process.T2RoadColl.Useonehalf = cms.bool(False)
process.T2RoadColl.DeltaR=5.0
process.T2RoadColl.DeltaPhi=4.5
process.T2RoadColl.DeltaZ = 400.0                       #default has 300

process.T2RoadColl.MinHits = 3   #non standard (5)      #default has 4
#process.T2RoadColl.MaxHits = 10
process.T2RoadColl.IncludeClass0PadCluster=cms.bool(True)

process.T2TrackColl2.UseRefittingProcedure=cms.bool(True)
process.T2TrackColl2.forceRZfit=False
process.T2TrackColl2.forceXYfit=True
process.T2TrackColl2.MinHitFinalTrk=3 #non standard (5)  #default have 4
process.T2TrackColl2.MaxHitsInRoad=30                   #default have 15
process.T2TrackColl2.verbosity=False
process.T2TrackColl2.MinHitFinalTrkCl1Hit=3


process.load("TotemAnalysis.T2ValidRawData.T2RecoValidation2_cfi")

process.T2ValidRaw.OutputFile=cms.untracked.string('EffiFilename.root')


process.T2ValidRaw.MaxEvents=29999#14000 #1500#15000


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

process.T2ValidRaw.SelectedHalf= (0)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@





process.T2ValidRaw.OnlycorruptionAnalysis=cms.bool(False)
#To be set true only for simulation file
process.T2ValidRaw.simufile=cms.bool(False)

#Used for hit error definition
process.T2ValidRaw.UseJointProb=0 #1=true; 0=false

#Used only for alignment histograms
process.T2ValidRaw.HitNumb4Align=4
process.T2ValidRaw.MeasuredXYResol=0.05
process.T2ValidRaw.SHIFTprescale=1.0
process.T2ValidRaw.MaxStepalignstep=400	
process.T2ValidRaw.Idreferencedet=13
process.T2ValidRaw.AlignmentHitRMax=140.0

#Used only for efficiency calculation (reference track cuts)
process.T2ValidRaw.NumHitGood=4
process.T2ValidRaw.useRZforResol=0
process.T2ValidRaw.MaxPad=7
process.T2ValidRaw.MaxStrip=7
process.T2ValidRaw.MaxDphi=7.0	
process.T2ValidRaw.maxdphihit=5.5#2.5
process.T2ValidRaw.maxdrhit=5.5
process.T2ValidRaw.NoiseDphiMAX=30.0
process.T2ValidRaw.NoiseDrMAX=30.0


process.T2ValidRaw.MaxPadAllowedInQuarter=40 #default tested with Visualizer was 40
	  
#Be careful: the following cut define the efficiency and the value Effmaxdphihit Effmaxdrhit depend on your fit choosing
#dr dphi could be taken as dx-dy
process.T2ValidRaw.Effgoodhitnumber=4 #era 7
process.T2ValidRaw.EffMaxPad=7      #era 9
process.T2ValidRaw.EffMaxStrip= 15     #era 9
process.T2ValidRaw.Effmaxdrhit=2.0   #final cut for efficiency calculation (this is really R)
process.T2ValidRaw.Effmaxdphihit=4.0	#final cut for efficiency calculation (this is really Phi)

process.T2ValidRaw.chiRredCut=100.0  
process.T2ValidRaw.chiPhiredCut=100.0
process.T2ValidRaw.AllowedDRTrackDistance=1.0

process.T2ValidRaw.verbosity=cms.bool(False)
process.T2ValidRaw.produceVfatEffiFile=cms.bool(True)
#/afs/cern.ch/exp/totem/scratch/berretti/tmp/testSplitMerge/CMSSW_3_1_1/src
process.T2ValidRaw.xmlfilenameFull=cms.string('TotemRawData/RawToDigi/python/T2GeoMapIP5_4quarter_vmea_Full.xml')
process.T2ValidRaw.xmlfilenameUsed_NotDead=cms.string('TotemRawData/RawToDigi/python/T2GeoMapIP5_4quarter_vmea_VFATsAlive_noTMC-14apr10.xml')

process.T2ValidRaw.MaxTrkInProcess=1


process.T2ValidRaw.skipSelectedEvents=cms.bool(False)
process.T2ValidRaw.skipEventFileName=cms.string('TotemAnalysis/T2ValidRawData/test/CorruptedEv.txt')

process.T2ValidRaw.LookToRawEvent=cms.bool(False)
process.T2ValidRaw.UseUncorrupetdEventMap=cms.bool(True)

# DQM modules
#HT 0-1 Arm0 |---| HT 2-3 Arm1
#RecoWithAlign/Reco_InternalAlCorr_MergedHalf0TrackFindingXY_2.py

#process.p = cms.Path(process.T2RoadCollH0*process.T2RoadCollH1*process.T2TrackColl2H0*process.T2TrackColl2H1)
process.p = cms.Path(process.RawToDigi*process.T2MCl*process.T2Hits*process.T2RoadColl*process.T2TrackColl2*process.T2ValidRaw)
#process.p = cms.Path(process.RawToDigi*process.T2MCl*process.T2Hits*process.T2RoadColl*process.T2TrackColl2*process.T2TrackColl2VTX)
#process.p = cms.Path(process.RawToDigi*process.T2MCl)
process.outpath = cms.EndPath(process.o1)
