import FWCore.ParameterSet.Config as cms

process = cms.Process("HighEnergyReco")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMax_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(15000)
)

process.source = cms.Source("PoolSource",  
    fileNames = cms.untracked.vstring('file:TotemRawData/RawToDigi/python/run_151_BeforeHit20000.root')
)

process.o1 = cms.OutputModule("PoolOutputModule",
     fileName = cms.untracked.string('file:TotemRawData/RawToDigi/python/outp.root')
)


process.load("RecoTotemT1T2.T2MakeCluster.T2MakeCluster_cfi")

process.load("RecoTotemT1T2.T2RecHit.T2RecHit_cfi")

#process.load("RecoTotemT1T2.T2RoadProducer.T2RoadCollDr15Dz200Dphi17_cfi")

process.load("RecoTotemT1T2.T2TrackProducer2.T2TrackColl25H_cfi")


	
#process.T2Hits.Cl1MaxPad =20
#process.T2Hits.Cl1MaxStrip =20
#process.T2Hits.IncludeClass0Hits = True

#process.T2RoadColl.DeltaR=4.
#process.T2RoadColl.DeltaRSecond=4.
#process.T2RoadColl.DeltaPhi=3.4
        

#process.T2RoadColl.DeltaPhi = 6.0
#process.T2RoadColl.DeltaZ = 450.0

#process.T2RoadColl.MinHits = 2          
#process.T2RoadColl.MaxHits = 20


#process.T2TrackColl2.verbosity= False #cms.bool(False)
#process.T2TrackColl2.forceRZfit=True
#process.T2TrackColl2.forceXYfit=False
#process.T2TrackColl2.MaxHitsInRoad=30





process.load("RecoTotemT1T2.T2RoadProducer.T2RoadCollDr15Dz200Dphi17_cfi")



#process.load("RecoTotemT1T2.T2TrackProducer2.T2TrackColl2_cfi")


process.T2Hits.Cl1MaxPad = cms.uint32(150)
process.T2Hits.Cl1MaxStrip = cms.uint32(150)
process.T2Hits.IncludeClass0Hits = True
process.T2Hits.inputFileNameMisal=cms.untracked.string('/home/mirko/SL/WorkingArea/CMSSW_311SVN2/CMSSW_3_1_1/src/TotemAlignment/T2TrkBasedInternalAlignment/test/MisalOut/Run_647_20000FalseChk_conf00.dat')
process.T2Hits.useTXTfile=cms.bool(True)
process.T2Hits.InsertAlignmentbyCFG=cms.bool(True)
process.T2Hits.verbosity=cms.untracked.bool(True)



process.T2RoadColl.DeltaR=4.0
process.T2RoadColl.DeltaPhi=3.0
process.T2RoadColl.DeltaZ = 300.0

process.T2RoadColl.MinHits = 4 
#process.T2RoadColl.MaxHits = 10 Not Used Anymore.

process.T2TrackColl2.forceRZfit=False
process.T2TrackColl2.forceXYfit=True
process.T2TrackColl2.MinHitFinalTrk=4
process.T2TrackColl2.MaxHitsInRoad=20
process.T2TrackColl2.verbosity=False
# DQM modules
#HT 0-1 Arm0 |---| HT 2-3 Arm1
#RecoWithAlign/Reco_InternalAlCorr_MergedHalf0TrackFindingXY_2.py






process.p = cms.Path(process.T2Hits*process.T2RoadColl*process.T2TrackColl2) 
#*process.dqm1
process.outpath = cms.EndPath(process.o1)
