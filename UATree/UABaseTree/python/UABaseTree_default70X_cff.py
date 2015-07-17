import FWCore.ParameterSet.Config as cms


process = cms.Process("UABaseTree")

process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(-1)
)

# initialize MessageLogger and output report ----------------------------------------

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.cout.placeholder = cms.untracked.bool(False)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
   SkipEvent = cms.untracked.vstring('ProductNotFound')

)

# Data source -----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
		'root://xrootd.ba.infn.it//store/mc/Spring14dr/MinBias_TuneMonash13_13TeV-pythia8/AODSIM/castor_NoPileUp_POSTLS170_V6-v1/00000/00ADED10-63F7-E311-AB96-003048F0EBBE.root'
   )
   )                               


# configure modules via Global Tag --------------------------------------------------
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START70_V7::All'


#Geometry --------------------------------------------------------------------------
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")


#process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

#HLT Stuff
process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

from UATree.UABaseTree.UABaseTree_cfi import *
process.load("UATree.UABaseTree.UABaseTree_cfi")

#FileName of the output of the UATree. Needed here for AutoCrab.
process.uabasetree.outputfilename = cms.untracked.string("UABaseTree_MC_TuneMonash13_13TeV-pythia8.root")

process.path = cms.Sequence()


#  ----------------------------   Filters   ----------------------------
#process.load('UATree.UABaseTree.UABaseTree_filters_cfi')

if not(isMonteCarlo):
   process.path = cms.Sequence(process.path * process.hltPhysicsDeclared)

if useMITFilter:
   process.path = cms.Sequence(process.path * process.mitfilter)
else:
   if not(isMonteCarlo):
     process.path = cms.Sequence(process.path * process.noscraping)


 
#  ----------------------------  Vertex  ----------------------------

if doDAvertex:
  process.load('UATree.UABaseTree.UABaseTree_AnnealingVertex_cfi')
  process.path = cms.Sequence(process.path * process.offlinePrimaryVerticesDA )



#  ----------------------------   Tracking   ----------------------------
if doMBTracking:
   if not(useMITFilter):
      process.load('UATree.UABaseTree.UABaseTree_tracking_cfi')
      process.path = cms.Sequence(process.path * process.redoSiHits )
   process.pixelTertTracks.OrderedHitsFactoryPSet.maxElement = cms.uint32( 3000 )
   process.path = cms.Sequence(process.path  * process.fulltracking)


#  ----------------------------   Jets   ----------------------------
if storeJets:
   process.load("UATree.UABaseTree.UABaseTree_jets_cfi")
#   process.path = cms.Sequence(process.path * process.L1FastJet )



#  ----------------------------   Castor   ----------------------------
#if storeCastor:
   #process.castorDigis.InputLabel = 'source'
   #process.load("RecoLocalCalo.Castor.CastorCellReco_cfi")    #-- redo cell
   #process.load("RecoLocalCalo.Castor.CastorTowerReco_cfi")   #-- redo tower
   #process.load("RecoJets.JetProducers.ak7CastorJets_cfi")    #-- redo jet
   #process.load("RecoJets.JetProducers.ak7CastorJetID_cfi")   #-- redo jetid
   #process.path = cms.Sequence(process.path  * process.castorDigis*process.castorreco*process.CastorCellReco*process.CastorTowerReco*process.ak7BasicJets*process.ak7CastorJetID)



#  ----------------------------  UE Stuff
process.load('UATree.UABaseTree.UEAnalysisTracks_cfi')
process.load('UATree.UABaseTree.UEAnalysisJetsSISCone_cfi')
process.load("UATree.UABaseTree.UEAnalysisJetsAk_cfi")
process.load("QCDAnalysis.UEAnalysis.UEAnalysisParticles_cfi")
process.ueSisCone5TracksJet500.DzTrVtxMax = cms.double(1000)
process.ueSisCone5TracksJet500.DxyTrVtxMax = cms.double(1000)
process.chargeParticles.cut = cms.string('charge != 0 & pt > 0.5 & status = 1')

process.path = cms.Sequence(process.path * process.UEAnalysisTracks * process.ueSisCone5TracksJet500 * process.UEAnalysisJetsAkOnlyReco)
#process.path = cms.Sequence(process.path * process.ueSisCone5TracksJet500)
if isMonteCarlo:
  process.path = cms.Sequence(process.path * process.UEAnalysisParticles * process.ueSisCone5ChgGenJet500 * process.UEAnalysisJetsAkOnlyMC)

#process.load('UATree.UABaseTree.UEJetChecker_cfi')
#process.path = cms.Sequence(process.path * process.uejetchecker)


# --------------------------- Write CMS data -------------------------------

if keepCMSData:
   process.out = cms.OutputModule("PoolOutputModule",
       fileName = cms.untracked.string('/tmp/katsas/cmsdata.root')
   )

#  ----------------------------   Final Path   ----------------------------
process.path = cms.Sequence(process.path * process.uabasetree)

process.p = cms.Path( process.path )

if keepCMSData:
   process.outpath = cms.EndPath(process.out)
