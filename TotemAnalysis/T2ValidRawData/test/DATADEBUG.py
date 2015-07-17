import FWCore.ParameterSet.Config as cms

process = cms.Process("valT1T2mu")

########################### COMMON PART ##########################################

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100) # -1 means to take all the events from the source as input
)


pool_output_dir ='/castor/cern.ch/totem/offline/Reco/2011/T2/May_PHYS_RUNS_Tracks_NTuple_2KeV_AlignmentConfBB_Geo91mm_GhostSuppression/run_run_5530.006.vmea.root__startEvt0.root'



finalCastordir = pool_output_dir + '/' 

commandls = 'nsls ' + finalCastordir + ' | grep .root '

print (commandls)
print '@@@@@@@@@@@@@@@ \n'
#echo $finalCastordir
import commands
stat, out = commands.getstatusoutput(commandls)
if not stat:
    print out


#Here I split the list of files

import string
lines = string.split(out, '\n')

new_list = []

for oneline in lines:
 print 'rfio:', oneline
 OnecastorFilePath = 'rfio:'  + finalCastordir + '/' + oneline
 new_list += [OnecastorFilePath]

#print '@@@@@@@@@@@@@@@ \n'

for oneline2 in new_list:
 print oneline2


#Here I feed the cmssw process new_list
#'file:prodT1T2_PythiaMirkobeta2energy3.5TeV_VFatEffi_Tune2_VtxXY_Dr25Dphi3Dz400-9.root'
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/totem/offline/Reco/2011/T2/May_PHYS_RUNS_Tracks_NTuple_2KeV_AlignmentConfBB_Geo91mm_GhostSuppression/run_run_5530.006.vmea.root__startEvt0.root"),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')                        
)


process.load("Configuration.TotemCommon.LoggerMax_cfi")

#

# random number generators seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")
process.load("TotemAnalysis.T2ValidRawData.T2RecoValidation2_cfi")


process.T2ValidRaw.OutputFile=cms.untracked.string('TestEffi.root') #gunT1T2muAnB.root

process.T2ValidRaw.MaxEvents=95#14000 #1500#15000


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

process.T2ValidRaw.SelectedHalf=0

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



#To be set true only for simulation file
process.T2ValidRaw.simufile=cms.bool(False)
process.T2ValidRaw.LookToRawEvent=cms.bool(False)
process.T2ValidRaw.UseUncorrupetdEventMap=cms.bool(False)

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
process.T2ValidRaw.MaxPadAllowedInQuarter=40
 
	  
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

process.T2ValidRaw.verbosity=cms.bool(True)
process.T2ValidRaw.produceVfatEffiFile=cms.bool(True)
process.T2ValidRaw.xmlfilenameFull=cms.string('TotemAnalysis/T2ValidRawData/test/T2GeoMapIP5_4quarter_vmea_Full.xml')
process.T2ValidRaw.xmlfilenameUsed_NotDead=cms.string('TotemAnalysis/T2ValidRawData/test/T2GeoMapIP5_4quarter_vmea_VFATsAlive_noTMC-14apr10.xml')
process.T2ValidRaw.DeadSectFileName=cms.string('TotemAnalysis/T2ValidRawData/test/3706_DeadSector.txt')
#process.T2ValidRaw.MaxTrkInProcess=1
process.T2ValidRaw.TrackLabel=cms.string('T2TrackColl3')
process.T2ValidRaw.OnlyClusterAnalysis=cms.bool(False)
process.T2ValidRaw.VFATMonitoring=cms.bool(False)

process.T2ValidRaw.MinTrkInQuarter= 1
process.T2ValidRaw.MaxTrkInQuarter= 1000


process.T2ValidRaw.skipSelectedEvents=cms.bool(False)
process.T2ValidRaw.skipEventFileName=cms.string('TotemAnalysis/T2ValidRawData/test/CorruptedEv.txt')


process.p1 = cms.Path(process.T2ValidRaw)

#process.outpath = cms.EndPath(process.o1)
 
