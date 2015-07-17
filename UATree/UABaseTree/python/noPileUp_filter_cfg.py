import FWCore.ParameterSet.Config as cms

process = cms.Process("noPileUpFilter")

process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(100)
)


# initialize MessageLogger and output report ----------------------------------------

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.cout.placeholder = cms.untracked.bool(False)
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(False)
   #SkipEvent = cms.untracked.vstring('ProductNotFound')

)

process.source = cms.Source("PoolSource",
   #fileNames = cms.untracked.vstring('dcap:///pnfs/iihe/cms/ph/sc4/store/mc/Winter10/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6/GEN-SIM-RECO/E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/0000/8C1DA622-5D26-E011-9F60-001D0967DC42.root')
   fileNames = cms.untracked.vstring('dcap:///pnfs/iihe/cms/store/user/xjanssen/data//CMSSW_3_9_7/DataCopy_397/__GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6__Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1__GEN-SIM-RECO/DataCopy_397__CMSSW_3_9_7__GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6__Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1__GEN-SIM-RECO_1_1_RJI.root')
)
#process.load("UATree.UABaseTree.noPileUp_filter_cfi")
process.nopileupfilter = cms.EDFilter('noPileUp')


process.out = cms.OutputModule("PoolOutputModule",
       fileName = cms.untracked.string('noPileUpFiltered.root')
)

process.p = cms.Path(process.nopileupfilter)
process.outpath = cms.EndPath(process.out)
