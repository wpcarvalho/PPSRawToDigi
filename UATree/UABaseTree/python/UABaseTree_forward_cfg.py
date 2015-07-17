import FWCore.ParameterSet.Config as cms

#from UATree.UABaseTree.TightPFJetID_Parameters_cfi   import TightPFJetID_Parameters as   TightPFJetID_Parameters_Ref
#from UATree.UABaseTree.LooseCaloJetID_Parameters_cfi import LooseCaloJetID_Parameters as LooseCaloJetID_Parameters_Ref
#from UATree.UABaseTree.TightCaloJetID_Parameters_cfi import TightCaloJetID_Parameters as TightCaloJetID_Parameters_Ref

from UABaseTree_forward_options import *

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
                            fileNames = cms.untracked.vstring()
   )                               

# configure modules via Global Tag --------------------------------------------------
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '::All'

#Geometry --------------------------------------------------------------------------
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

#process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
#process.load('Configuration.StandardSequences.GeometryExtended_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

#HLT Stuff
process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

process.load("UATree.UABaseTree.UABaseTree_forward_cfi")

#FileName of the output of the UATree. Needed here for AutoCrab.
process.uabasetree.outputfilename = cms.untracked.string("UABaseTree_CMS-TOTEM.root")

#  ----------------------------   Filters   ----------------------------
process.load('UATree.UABaseTree.UABaseTree_filters_cfi')

#  ----------------------------   Jets   ----------------------------
if storeJets:
   process.load("UATree.UABaseTree.UABaseTree_jets_cfi")

#  ----------------------------  UE Stuff
process.load('UATree.UABaseTree.UEAnalysisTracks_cfi')
process.load('UATree.UABaseTree.UEAnalysisJetsSISCone_cfi')
process.load("UATree.UABaseTree.UEAnalysisJetsAk_cfi")
process.load("QCDAnalysis.UEAnalysis.UEAnalysisParticles_cfi")
process.ueSisCone5TracksJet500.DzTrVtxMax = cms.double(1000)
process.ueSisCone5TracksJet500.DxyTrVtxMax = cms.double(1000)
process.chargeParticles.cut = cms.string('charge != 0 & pt > 0.5 & status = 1')

# --------------------------- Write CMS data -------------------------------
if keepCMSData:
   process.out = cms.OutputModule("PoolOutputModule",
       fileName = cms.untracked.string('cmsdata.root')
   )

#  ----------------------------   Final Path   ----------------------------
process.selection_sequence = cms.Sequence()
process.reco_sequence = cms.Sequence()
#process.selection_sequence = cms.Sequence(process.selection_sequence * process.hltPhysicsDeclared)

if useMITFilter:
   process.selection_sequence = cms.Sequence(process.selection_sequence * process.mitfilter)
else:
   process.selection_sequence = cms.Sequence(process.selection_sequence * process.noscraping)

#  ----------------------------  Vertex  ----------------------------
if doDAvertex:
  process.load('UATree.UABaseTree.UABaseTree_AnnealingVertex_cfi')
  process.reco_sequence = cms.Sequence( process.reco_sequence * process.offlinePrimaryVerticesDA )

#  ----------------------------   Tracking   ----------------------------
if doMBTracking:
    if not(useMITFilter):
        process.load('UATree.UABaseTree.UABaseTree_tracking_cfi')
        process.reco_sequence= cms.Sequence( process.reco_sequence * process.redoSiHits )

    process.pixelTertTracks.OrderedHitsFactoryPSet.maxElement = cms.uint32( 3000 )
    process.reco_sequence = cms.Sequence( process.reco_sequence * process.fulltracking )

#  ----------------------------   Jets   ----------------------------
#if storeJets:
#   process.reco_sequence = cms.Sequence( process.reco_sequence * process.L1FastJet )

#  ----------------------------   Castor   ----------------------------
#if storeCastor:
   #process.castorDigis.InputLabel = 'source'
   #process.load("RecoLocalCalo.Castor.CastorCellReco_cfi")    #-- redo cell
   #process.load("RecoLocalCalo.Castor.CastorTowerReco_cfi")   #-- redo tower
   #process.load("RecoJets.JetProducers.ak7CastorJets_cfi")    #-- redo jet
   #process.load("RecoJets.JetProducers.ak7CastorJetID_cfi")   #-- redo jetid
   #process.reco_sequence = cms.Sequence(process.reco_sequence  * process.castorDigis*process.castorreco*process.CastorCellReco*process.CastorTowerReco*process.ak7BasicJets*process.ak7CastorJetID)

#process.load('UATree.UABaseTree.UEJetChecker_cfi')
#process.reco_sequence = cms.Sequence(process.reco_sequence * process.uejetchecker)

#  ----------------------------   ZDC   ----------------------------
if storeZDC:
    ##-- Customized conditions
    from CondCore.DBCommon.CondDBSetup_cfi import *
    process.hcalZDCASCIICalib = cms.ESSource("HcalTextCalibrations",
	    input = cms.VPSet(
		cms.PSet(
		    object = cms.string('Gains'),
		    file =
		    cms.FileInPath('UATree/UABaseTree/data/Gains_Run210498_211831.txt')
		    ),
		cms.PSet(
		    object = cms.string('LongRecoParams'),
		    file =
		    cms.FileInPath('UATree/UABaseTree/data/LongRecoParams_Runs210737_211831.txt')
		    ),
		cms.PSet(
		    object = cms.string('MCParams'),
		    file =
		    cms.FileInPath('UATree/UABaseTree/data/DumpMCParams_Run211831.txt')
		    )
		)
	    )
    process.es_prefer_hcalZDCASCIICalib = cms.ESPrefer('HcalTextCalibrations','hcalZDCASCIICalib')
    process.zdcreco.lowGainFrac = 8.6
    process.reco_sequence = cms.Sequence(process.reco_sequence * process.zdcreco)


process.reco_sequence = cms.Sequence(process.reco_sequence *
                                     process.UEAnalysisTracks * 
                                     process.ueSisCone5TracksJet500 * process.UEAnalysisJetsAkOnlyReco)

process.p = cms.Path(process.selection_sequence * process.reco_sequence * process.uabasetree)

if keepCMSData:
   process.outpath = cms.EndPath(process.out)
