import FWCore.ParameterSet.Config as cms

process = cms.Process("valT1T2mu")

########################### COMMON PART ##########################################

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(15000) # -1 means to take all the events from the source as input
)



# Configure the input module (read the events from a file -- gunT1T2mu.root)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/totem/offline/Reco/2010/T2/run_1552Reco05Jul010_NoClean_IntGlobAL24666_.root')
)
#('file:/media/usbdisk1/RecoDataApr010/Recorun647XYIntAL_conf00.root')
#('file:/afs/cern.ch/exp/totem/scratch/data/run_647.root')

process.load("Configuration.TotemCommon.LoggerMin_cfi")



# random number generators seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")
process.load("TotemAnalysis.T2ValidRawData.T2RecoValidation2_cfi")

process.T2ValidRaw.OutputFile=cms.untracked.string('outputEff.root')
process.T2ValidRaw.MaxEvents=14000 #1500#15000


#To be set true only for simulation file
process.T2ValidRaw.simufile=cms.bool(False)

#Used for hit error definition
process.T2ValidRaw.UseJointProb=0 #1=true; 0=false

#Used only for alignment histograms
process.T2ValidRaw.HitNumb4Align=6
process.T2ValidRaw.MeasuredXYResol=0.05
process.T2ValidRaw.SHIFTprescale=1.0
process.T2ValidRaw.MaxStepalignstep=400	
process.T2ValidRaw.Idreferencedet=13
process.T2ValidRaw.AlignmentHitRMax=140.0

#Used only for efficiency calculation (reference track cuts)
process.T2ValidRaw.NumHitGood=5	 
process.T2ValidRaw.useRZforResol=0
process.T2ValidRaw.MaxPad=7
process.T2ValidRaw.MaxStrip=7
process.T2ValidRaw.MaxDphi=7.0	
process.T2ValidRaw.maxdphihit=5.5#2.5
process.T2ValidRaw.maxdrhit=5.5
process.T2ValidRaw.NoiseDphiMAX=30.0
process.T2ValidRaw.NoiseDrMAX=30.0
 
	  
#Be careful: the following cut define the efficiency and the value Effmaxdphihit Effmaxdrhit depend on your fit choosing
#dr dphi could be taken as dx-dy
process.T2ValidRaw.Effgoodhitnumber=5 #era 7
process.T2ValidRaw.EffMaxPad=7      #era 9
process.T2ValidRaw.EffMaxStrip= 15     #era 9
process.T2ValidRaw.Effmaxdrhit=2.0   #final cut for efficiency calculation (this is really R)
process.T2ValidRaw.Effmaxdphihit=4.0	#final cut for efficiency calculation (this is really Phi)

process.T2ValidRaw.chiRredCut=100.0  
process.T2ValidRaw.chiPhiredCut=100.0
process.T2ValidRaw.AllowedDRTrackDistance=1.0

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

process.T2ValidRaw.SelectedHalf=0

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

process.p1 = cms.Path(process.T2ValidRaw)

#process.outpath = cms.EndPath(process.o1)
