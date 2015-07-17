import FWCore.ParameterSet.Config as cms

import sys, random, stat, os, os.path

import string


process = cms.Process("valT1T2pi")

########################### COMMON PART ##########################################

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:/tmp/berretti/run_prodT2_Pythia8Inelastic_E3.5TeV_ONLYGEANT_Nov2011_91mmSamp4_RELEASE424-95.root__startEvt0.root')
# )

# logging to txt files
#process.load("Configuration.TotemCommon.LoggerMax_cfi")

# random number generators seeds
#process.load("Configuration.TotemCommon.RandomNumbers_cfi")

#pool_output_dir = '/castor/cern.ch/totem/offline/mcdata/Py8_8TeV_Test_REL445SimAndRec-O/' run_97_.root

nameF='QGSJET01_8TeV-Vtx0-424-T-CorrectCalo-TestAsym_424RecoT1T2GhostSuppressionAlignV1_EffiClu_MultSTD'
nameF='Py8_8TeV-Vtx0-424-T-CorrectCalo-TestAsym_424RecoT1T2GhostSuppressionAlignV2_EffiClu_MultSTD'
nameF='EPOS_LHCRetune_8TeV-Vtx0-Pure445-CorrectCalo-LargeSamp0_445PureRecoT1T2GhostSuppressionAlignV2-Tuned8372Effi'

nameF='EPOS_LHCRetune_8TeV-Vtx0-424-CorrectCalo-LargeSamp0'
nameF='EPOS_LHCRetune_8TeV-Vtx0-Pure445-CorrectCalo-LargeSamp0'



#nameF='Phojet-DPE-Vtx0-424-T-CorrectCalo-LargeSamp0_424RecoT1T2GhostSuppressionAlignV2_EffiClu8372_MultSTD'
pool_output_dir = '/castor/cern.ch/user/b/berretti/'+nameF
pool_output_dir = '/castor/cern.ch/totem/offline/mcdata/'+nameF

finalCastordir = pool_output_dir 

commandls = 'nsls ' + finalCastordir + ' | grep .root | grep -v .py' 
print '@@@@@@@@@@@@@@@ \n'
#echo $finalCastordir
import commands
stat, out = commands.getstatusoutput(commandls)
if not stat:
    print out


#Here I split the list of files run_98_.root_reco.root

import string
lines = string.split(out, '\n')

new_list = []
#range(MinRun,MaxRun)

MinRun="0";
MaxRun="60";

MinR2=int(MinRun);
MaxR2=int(MaxRun);
    
for countJob in range(MinR2,MaxR2):
 OnecastorFilePath = 'rfio:'  + finalCastordir + '/run_'+str(countJob)+'_.root'
 print OnecastorFilePath
 new_list += [OnecastorFilePath]



process.source = cms.Source("PoolSource",
 
    fileNames = cms.untracked.vstring(new_list),                        
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')                        
)

DigiSeed=str(random.randrange(10000, 99998));
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",    
    moduleSeeds = cms.PSet(                                               
      T2Digis = cms.untracked.uint32(44444),
      T2MCl = cms.untracked.uint32(44444),
    ),
    sourceSeed = cms.untracked.uint32(44444),                                         
)




DigiMisalFilename="SimTotem/T2Digitizer/data/2012IntGlobAlignV2.dat"

DigiFNameForHit="SimTotem/T2Digitizer/data/2012IntGlobAlignV2.dat"

# process.load("SimTotem.T2Digitizer.T2DigisOnlySec_cfi")
process.load("SimTotem.T2Digitizer.T2DigisSTDGeoEffi_cfi")#T2DigisOnlyPrim_cfi
#process.load("SimTotem.T2Digitizer.T2DigisOnlyPrim_cfi")
process.T2Digis.saveDigiVFAT=cms.bool(True)
process.T2Digis.Misalignment.simulatemisalign=cms.untracked.bool(True)
#Decomment for alignment multi-scenario validation
process.T2Digis.Misalignment.inputFileNameMisal=DigiMisalFilename




process.load("RecoTotemT1T2.T2MakeCluster.T2MakeCluster_cfi")



process.T2MCl.TakeCleanEventOnly=cms.bool(False) #IMPORTANT
#process.T2MCl.maskvect = cms.vuint32(21,23,25,27,29)
process.T2MCl.SimuClusterEfficiency=cms.bool(True)
process.T2MCl.EffiGeoRootFileName= cms.string("RecoTotemT1T2/T2MakeCluster/data/Geom_effiOutput_All8372_Effi_V2Pl.root")
#8314_EffiGeoSamp20K.root" May2012RunTestH0123.root




process.load("RecoTotemT1T2.T2RecHit.T2RecHit_cfi")
#from RecoTotemT1T2.T2RecHit.T2RecHit_cfi import *
process.T2Hits.Cl1MaxPad = cms.uint32(25) #Tune better
process.T2Hits.Cl1MaxStrip = cms.uint32(25)
process.T2Hits.IncludeClass0Hits = True
#T2Hits.inputFileNameMisal=cms.untracked.string('/afs/cern.ch/exp/totem/scratch/berretti/tmp/testSplitMerge/CMSSW_3_1_1/src/$
#process.T2Hits.inputFileNameMisal=cms.untracked.string('SimTotem/T2Digitizer/data/run_5525-5535_IntAlignHX50000Evt_XovY0.3HIP_ANDGLOB_BBConf.dat')
process.T2Hits.inputFileNameMisal=cms.untracked.string("SimTotem/T2Digitizer/data/2012IntGlobAlignV2.dat")

process.T2Hits.useTXTfile=cms.bool(True)
process.T2Hits.InsertAlignmentbyCFG=cms.bool(True) # True for data
process.T2Hits.verbosity=cms.untracked.bool(False)
process.T2Hits.CorrectWithResolution=cms.bool(True) #False:Old Strategy True:New Strategy

#from RecoTotemT1T2.T2RoadPadFinder.NewLabelT2RoadPadFinder_cfi import *
process.load("RecoTotemT1T2.T2RoadPadFinder.NewLabelT2RoadPadFinder_cfi")#T2RoadPadFinder_cfi

process.T2RoadPadFinder.HitLabel=cms.string("T2Hits")#T2HitsSTD for effi reco |  #T2Hits for std reco
process.T2RoadPadFinder.CluLabel=cms.string("T2MCl")#T2MClSTD for effi reco |   T2MCl  for std reco

process.T2RoadPadFinder.verbosity = 0
process.T2RoadPadFinder.TwoPointsTubesAngularCollinearity=0.09#0.07 default
process.T2RoadPadFinder.MinCluSize_considered_asBlobs = cms.int32(5)
process.T2RoadPadFinder.MinimumNumCl1Hit= 3
process.T2RoadPadFinder.chi2XYProb_Thr= 0.01
process.T2RoadPadFinder.Nmin_padsFinal= 4
process.T2RoadPadFinder.T2RoadCollProdName="NewRoadFinderRELOAD"



process.T2RoadPadFinder.AllowsPadReAssociation=False
process.T2RoadPadFinder.AllowsConcurrentBranches=False
process.T2RoadPadFinder.useStraightPadTowers= cms.bool(True)#False for internal alignment studies
process.T2RoadPadFinder.ResolveOverlapDoubleCount = cms.bool(False) #Default is True, False for shadow alignment and dndeta An
process.T2RoadPadFinder.OverlapDoubleCountDR = cms.double(2.0) #Depend on your alignment Resol
process.T2RoadPadFinder.OverlapDoubleCountDPhi =cms.double(3.5)
process.T2RoadPadFinder.OverlapDoubleCountDTheta =  cms.double(0.01)


###########################################################################
#process.T2RoadPadFinder.VplaneToExclude = cms.vint32(9,19,29,39)
###########################################################################


process.T2RoadPadFinder.QuartedSelected = cms.vint32(0,1,2,3)
process.T2RoadPadFinder.BiggestTubeAngleConsidered =cms.double(0.3)
process.T2RoadPadFinder.NumSigma= cms.double(2.)#Important for ALignment (6) istead of 2
#TolleranceDX
process.T2RoadPadFinder.NumPadCluOccupancyAlert= cms.double(50.)#Default is 50
process.T2RoadPadFinder.InefficiencyMaxJump= cms.int32(3)#3  is default
process.T2RoadPadFinder.Nmin_padsFinal= 4

process.load("RecoTotemT1T2.T2TrackProducer3.T2TrackColl3_cfi")
process.load("TotemT1T2Validation.T2RecoValidation.T2RecoValidation2_cfi")

#from RecoTotemT1T2.T2TrackProducer3.T2TrackColl3_cfi import *
process.T2TrackColl3.StripFitting=cms.bool(False)
process.T2TrackColl3.RoadModuleLabel="T2RoadPadFinder"
process.T2TrackColl3.RoadInstanceLabel="NewRoadFinderRELOAD"
process.T2TrackColl3.verbosity=False
process.T2TrackColl3.RemoveOutliers=True #False for Internal  ALignment studies
process.T2TrackColl3.GhostSuppression=True
process.T2TrackColl3.PickUpDisplacedHit=False
process.T2TrackColl3.PickUpRadius=2.0

outputfilename = '/tmp/berretti/'+nameF+'_'+str(MinRun)+'_'+str(MaxRun)+'.root'
process.o1 = cms.OutputModule("PoolOutputModule",                 
    fileName = cms.untracked.string(outputfilename),
    outputCommands = cms.untracked.vstring('keep *', 'drop *_*_*TrackerHits*_*', 'drop *_*_*Muon*_*', 'drop *_*_*Ecal*_*', 'drop *_*_*Hcal*_*', 'drop *_*_*Calo*_*', 'drop *_*_*Castor*_*', 'drop *_*_*FP420SI_*', 'drop *_*_*ZDCHITS_*', 'drop *_*_*BSCHits_*', 'drop *_*_*ChamberHits_*', 'drop *_*_*FibreHits_*', 'drop *_*_*WedgeHits_*')                           
   #    outputCommands = cms.untracked.vstring('drop TotemRawEvent_*_*_*','keep T2DetIdT2DigiVfatTotemDigiCollection_*_*_*','keep *_*_*T2VfatInformation*_*','keep T1T2Tracks_*_*_*','keep T2Hits_*_*_*','keep T2Roads_*_*_*','keep T2DetIdT2Clustersstdmap_*_*_*','keep T2DetIdT2StripDigiTotemDigiCollection_*_*_*','keep T2DetIdT2PadDigiTotemDigiCollection_*_*_*')    
    )

 
process.load("TotemT1T2Validation.T2RecoValidation.T2RecoValidation2_cfi")
process.T2RecoAn2.OutputFile = nameF+'_'+str(MinRun)+'_'+str(MaxRun)+'.root'
#'QGSJET01_8TeV-Vtx0-424-T-CorrectCalo-TestAsym_424RecoT1T2GhostSuppressionAlignV1_EffiClu_MultSTD.root'
#process.p1 = cms.Path(process.T2RecoAn2)

process.p1 = cms.Path(process.T2Digis*process.T2MCl*process.T2Hits*process.T2RoadPadFinder*process.T2TrackColl3*process.T2RecoAn2)
   

process.outpath = cms.EndPath(process.o1)
  
