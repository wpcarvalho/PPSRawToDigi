# No pile up for the mixing module
# from SimGeneral.MixingModule.mixNoPU_cfi  import *
from Configuration.TotemCommon.mixNoPU_cfi import *

########################### DIGI + RECO T2 ##########################################

from SimTotem.T2Digitizer.T2Digis_TuneG_5525_5535_May2011Effi_Internal_GlobalMisalBBConf_cfi import *
#2011-Tune Digitized data.
#use T2Digis for generic purpose simulation

from RecoTotemT1T2.T2MakeCluster.T2MakeCluster_cfi import *
T2MCl.maskvect = cms.vuint32()  #no masking by default
#Use T2MCl.maskvect = cms.vuint32(21,23,25,27,29) For the 2011 (After June) Run


from RecoTotemT1T2.T2RecHit.T2RecHit_cfi import *

T2Hits.inputFileNameMisal=cms.untracked.string('SimTotem/T2Digitizer/data/run_5525-5535_IntAlignHX50000Evt_XovY0.3HIP_ANDGLOB_BBConf.dat')# put useTXTfile and InsertAlignmentbyCFG to False in order to avoid Hit alignment correction.
T2Hits.useTXTfile=cms.bool(False)   
T2Hits.InsertAlignmentbyCFG=cms.bool(False)  
T2Hits.verbosity=cms.untracked.bool(False)
T2Hits.CorrectWithResolution=cms.bool(True) #False:Old Strategy

from RecoTotemT1T2.T2RoadPadFinder.NewLabelT2RoadPadFinder_cfi import *

T2RoadPadFinder.T2RoadCollProdName="NewRoadFinderRELOAD"
T2RoadPadFinder.useStraightPadTowers= cms.bool(True)#False: Old Strategy
T2RoadPadFinder.ResolveOverlapDoubleCount = cms.bool(False) # False for shadow alignment and dndeta An
T2RoadPadFinder.BiggestTubeAngleConsidered =cms.double(0.3)


from RecoTotemT1T2.T2TrackProducer3.T2TrackColl3_cfi import *
T2TrackColl3.RoadModuleLabel="T2RoadPadFinder"
T2TrackColl3.RoadInstanceLabel="NewRoadFinderRELOAD"
T2TrackColl3.RemoveOutliers=True
T2TrackColl3.GhostSuppression=True


