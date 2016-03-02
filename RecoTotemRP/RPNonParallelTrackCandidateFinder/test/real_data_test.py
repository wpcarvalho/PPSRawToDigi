import FWCore.ParameterSet.Config as cms

process = cms.Process("NonParallelTrackFinderTest")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
#    input = cms.untracked.int32(100)
)


# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")
#process.load("Configuration.TotemCommon.LoggerMax_cfi")


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",    
    moduleSeeds = cms.PSet(                                               
      T2Digis = cms.untracked.uint32(14142)
    ),
    sourceSeed = cms.untracked.uint32(173205)                                               
)


process.load("Configuration/TotemOpticsConfiguration/OpticsConfig_3500GeV_1p5_120urad_thin_cfi")


process.TotemRPIncludeAlignments = cms.ESProducer("TotemRPIncludeAlignments",
    RealFiles = cms.vstring('TotemAlignment/RPData/LHC/2011_05_18/sr+hsx/45_220.xml', 
        'TotemAlignment/RPData/LHC/2011_05_18/sr+hsx/56_220.xml'),
    MisalignedFiles = cms.vstring(),  MeasuredFiles = cms.vstring()
)


process.load('TotemRawData.Readers.RawDataSource_cfi')
process.source.verbosity = 1
process.source.printProgressFrequency = 0

process.source.fileNames = cms.untracked.vstring()
process.source.fileNames.append('/castor/cern.ch/totem/LHCRawData/2011/Physics/run_5601.000.vmea')


# raw to digi conversion
process.load('TotemCondFormats/DAQInformation/DAQMappingSourceXML_cfi')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_220.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_147.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/t1_all_run1.xml')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/t2_4quarters.xml')


process.load('TotemRawData.RawToDigi.Raw2DigiProducer_cfi')

# clusterization
process.load("RecoTotemRP.RPClusterSigmaService.ClusterSigmaServiceConf_cfi")
process.load("RecoTotemRP.RPClusterizer.RPClusterizationConf_cfi")

# reco hit production
process.load("RecoTotemRP.RPRecoHitProducer.RPRecoHitProdConf_cfi")


# geometry
process.load("Configuration/TotemCommon/geometryGlobal_real_cfi")
toberemoved = []
for xmlfile in process.XMLIdealGeometryESSource.geomXMLFiles:
    if xmlfile.endswith("RP_Dist_Beam_Cent.xml"):
       toberemoved.append(xmlfile)
for xmlfile in toberemoved:
    process.XMLIdealGeometryESSource.geomXMLFiles.remove(xmlfile)
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/2011_05_18/RP_Dist_Beam_Cent.xml")
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

# track search/pattern recognition
process.load("RecoTotemRP.RPNonParallelTrackCandidateFinder.RPNonParallelTrackCandidateFinder_cfi")
process.NonParallelTrackFinder.verbosity = 0
process.NonParallelTrackFinder.maxHitsPerPlaneToSearch = 4

# track fitting
process.load("RecoTotemRP.RPTrackCandidateCollectionFitter.RPSingleTrackCandCollFitted_cfi")
process.RPSingleTrackCandCollFit.Verbosity = 0

process.load("RecoTotemRP.RPInelasticReconstruction.Rec_3500GeV_beta_1p5_120urad_220_2Arm_cfi")
process.RP2202ArmReconst.BeamProtTransportSetup = process.BeamProtTransportSetup
process.RP2202ArmReconst.ExpectedRPResolution = 0.020 #mm
process.RP2202ArmReconst.Verbosity = 0

process.load("RecoTotemRP.RPInelasticReconstruction.Rec_3500GeV_beta_1p5_120urad_220_cfi")
process.RP220Reconst.BeamProtTransportSetup = process.BeamProtTransportSetup
process.RP220Reconst.ExpectedRPResolution = 0.020 #mm
process.RP220Reconst.Verbosity = 0
process.RP220Reconst.ElasticScatteringReconstruction = False

process.TriggerBits = cms.EDProducer("RPTriggerBitsProducer",
    verbose = cms.bool(False)
)

process.load('L1TriggerTotem.CoincidenceChip.RPCoincidenceProducer_cfi')
#process.load("TotemRawData.RawToDigi.RPDataCCProducer_cfi")


process.load("RecoTotemRP.RPMulCandidateTrackFinder.RPMulTrackCandFindConf_cfi")
process.load("RecoTotemRP.RPMulTrackCandidateCollectionFitter.RPMulTrackCandCollFitter_cfi")


# T2 ##############
#Fill T2 digi and vfat object
#process.RawToDigi = cms.EDProducer("T2XMLDataDigiProducer",
#  verbosity = cms.untracked.uint32(0),  #was 10
#  discardHighOccupancyVfatverbosity= cms.untracked.bool(False)#IMPORTANT
#)


process.load("SimTotem.T2Digitizer.T2Digis_TuneG_5525_5535_May2011Effi_Internal_GlobalMisalBBConf_cfi")
process.T2Digis.saveDigiVFAT=cms.bool(True) #False DEF


process.load("RecoTotemT1T2.T2MakeCluster.T2MakeCluster_cfi")
process.T2MCl.TakeCleanEventOnly=cms.bool(False) #IMPORTANT
process.load("RecoTotemT1T2.T2RecHit.T2RecHit_cfi")
process.T2Hits.Cl1MaxPad = cms.uint32(25) #Tune better
process.T2Hits.Cl1MaxStrip = cms.uint32(25)
process.T2Hits.IncludeClass0Hits = True
process.T2Hits.inputFileNameMisal=cms.untracked.string('SimTotem/T2Digitizer/data/run_5525-5535_IntAlignHX50000Evt_XovY0.3HIP_ANDGLOB_BBConf.dat')

process.T2Hits.useTXTfile=cms.bool(True)   #True for data
process.T2Hits.InsertAlignmentbyCFG=cms.bool(True) # True for data 
process.T2Hits.verbosity=cms.untracked.bool(False)
process.T2Hits.CorrectWithResolution=cms.bool(True) #False:Old Strategy


process.load("RecoTotemT1T2.T2RoadPadFinder.NewLabelT2RoadPadFinder_cfi")
process.T2RoadPadFinder.HitLabel=cms.string("T2Hits")
process.T2RoadPadFinder.CluLabel=cms.string("T2MCl")
process.T2RoadPadFinder.verbosity = 0
process.T2RoadPadFinder.TwoPointsTubesAngularCollinearity=0.07
process.T2RoadPadFinder.MinCluSize_considered_asBlobs = cms.int32(5)
process.T2RoadPadFinder.MinimumNumCl1Hit= 3
process.T2RoadPadFinder.chi2XYProb_Thr= 0.01
process.T2RoadPadFinder.Nmin_padsFinal= 4
process.T2RoadPadFinder.T2RoadCollProdName="NewRoadFinderRELOAD"
process.T2RoadPadFinder.AllowsPadReAssociation=False
process.T2RoadPadFinder.AllowsConcurrentBranches=False
process.T2RoadPadFinder.useStraightPadTowers= cms.bool(True)#False

#################################################################################################################
process.T2RoadPadFinder.ResolveOverlapDoubleCount = cms.bool(True) #Default is True, False for shadow alignment
#################################################################################################################

process.T2RoadPadFinder.OverlapDoubleCountDR = cms.double(2.0) #Depend on your alignment Resol  
process.T2RoadPadFinder.OverlapDoubleCountDPhi =cms.double(3.5)
process.T2RoadPadFinder.OverlapDoubleCountDTheta =  cms.double(0.01)

process.T2RoadPadFinder.QuartedSelected = cms.vint32(0,1,2,3)
process.T2RoadPadFinder.BiggestTubeAngleConsidered =cms.double(0.3)
process.T2RoadPadFinder.NumSigma= cms.double(2.)
process.T2RoadPadFinder.NumPadCluOccupancyAlert= cms.double(50.)
process.T2RoadPadFinder.InefficiencyMaxJump= cms.int32(3)#2 is default


process.load("RecoTotemT1T2.T2TrackProducer3.T2TrackColl3_cfi")
process.T2TrackColl3.StripFitting=cms.bool(False)
process.T2TrackColl3.RoadModuleLabel="T2RoadPadFinder"
process.T2TrackColl3.RoadInstanceLabel="NewRoadFinderRELOAD"
process.T2TrackColl3.verbosity=False
process.T2TrackColl3.RemoveOutliers=True



process.p = cms.Path(
      process.Raw2DigiProducer 
    * process.TriggerBits 
    * process.RPCC 
    * process.RPClustProd 
    * process.RPRecoHitProd
    * process.NonParallelTrackFinder 
    * process.RPSingleTrackCandCollFit
    * process.RP220Reconst
    * process.RP2202ArmReconst
    * process.RPMulTrackCandFind
    * process.RPMulTrackCandCollFit
    * process.T2MCl
    * process.T2Hits
    * process.T2RoadPadFinder
    * process.T2TrackColl3
)

# store desired results
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("file:output.root"),
    outputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.outpath = cms.EndPath(process.output)

