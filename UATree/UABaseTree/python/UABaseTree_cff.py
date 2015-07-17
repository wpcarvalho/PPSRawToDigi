import FWCore.ParameterSet.Config as cms

#from UATree.UABaseTree.TightPFJetID_Parameters_cfi   import TightPFJetID_Parameters as   TightPFJetID_Parameters_Ref
#from UATree.UABaseTree.LooseCaloJetID_Parameters_cfi import LooseCaloJetID_Parameters as LooseCaloJetID_Parameters_Ref
#from UATree.UABaseTree.TightCaloJetID_Parameters_cfi import TightCaloJetID_Parameters as TightCaloJetID_Parameters_Ref


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
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/00A2F6F1-8895-E111-9ADD-0025901D5DB8.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/00E7D27D-9095-E111-B172-0025B324400C.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/02762A56-9395-E111-AD5C-BCAEC532971D.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/02DA5B36-8C95-E111-A7CF-003048F117EA.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/049C146D-8D95-E111-ACC1-E0CB4E553673.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/0663B8AA-8D95-E111-AC1A-003048F117B6.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/08185ED6-8C95-E111-AFE5-0025901D5D90.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/0C62EC3B-8C95-E111-8367-001D09F276CF.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/0E350A56-9395-E111-A87D-BCAEC5329702.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/106405FD-8D95-E111-9663-003048D37524.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/10C65D96-8A95-E111-AC3C-001D09F29533.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/168B7D78-8695-E111-A1C5-0025901D5C80.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/1ACADEFE-8D95-E111-A9DA-BCAEC5364CFB.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/1CD3D3CA-8F95-E111-9851-003048F24A04.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/204C75D7-8C95-E111-80F6-485B3977172C.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/2095B8F5-9195-E111-B31C-5404A63886B2.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/2242F938-8795-E111-A2AF-E0CB4E5536AE.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/26130FD9-9295-E111-BAF0-00237DDC5CB0.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/2A352A5B-8995-E111-8E75-003048F024FA.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/2A518787-8E95-E111-9DD8-00237DDC5CB0.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/2EC9A1F5-8E95-E111-9463-003048D37694.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/32D17CA1-9295-E111-8BBB-001D09F2525D.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/34D0F938-8795-E111-9BB9-BCAEC518FF89.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/38367A56-9395-E111-9433-BCAEC5364C42.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/3A17450E-8895-E111-BCE7-0025901D624E.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/3A1A2DF6-9195-E111-9CEF-BCAEC5329709.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/3A91FE11-8895-E111-8D38-5404A640A643.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/3E754A39-8C95-E111-A39F-001D09F29146.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/401B2988-8E95-E111-9A85-001D09F2424A.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/40A2D038-8C95-E111-857E-001D09F29533.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/40A971F4-9195-E111-B72F-003048F11DE2.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/40E47F0D-9495-E111-883E-001D09F290CE.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/427D1A56-9395-E111-B47B-BCAEC518FF7C.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/42B58C35-8C95-E111-8010-003048F117B4.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/4643783A-8C95-E111-898D-001D09F24EE3.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/7E348F0D-8895-E111-B302-5404A63886A0.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/80109BFF-8895-E111-B417-5404A63886D6.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/80428984-8E95-E111-B147-003048D3C932.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/8424E1F1-8895-E111-AB95-0025901D624A.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/8613107C-8E95-E111-AFC1-0025901D631E.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/868C5FFF-8D95-E111-A2BA-002481E0DEC6.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/8A1273F4-9195-E111-B6D2-003048F024DA.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/8C628C00-8995-E111-9120-003048D37560.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/8C7A4DC1-9495-E111-AF31-002481E0D958.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/8CB756F6-9195-E111-ABC3-002481E0D90C.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/8E02C71D-8F95-E111-9C72-003048F024DE.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/8E346414-8895-E111-B8C1-002481E0D83E.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/8EFD9073-8D95-E111-8019-002481E0D790.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/9AE9B475-8D95-E111-8134-BCAEC532970D.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/9C99B10B-9495-E111-8572-5404A63886AF.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/9E331E6E-8D95-E111-A49B-BCAEC518FF52.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/A085BC1F-8F95-E111-B71A-001D09F27067.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/A0B92AA2-9295-E111-802C-0019B9F72CE5.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/A2256CD5-8C95-E111-B3F2-BCAEC5329716.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/A2286E11-8895-E111-86F7-BCAEC532970D.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/A2350C54-8B95-E111-BE94-001D09F29114.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/A6DA21A2-9295-E111-8F13-001D09F28F25.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/A8310E7B-8D95-E111-A6E3-5404A63886B2.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/AA45CA0D-8895-E111-BD2C-BCAEC5329717.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/AE568A0E-9495-E111-B958-001D09F29114.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/AEAB6D73-8D95-E111-9958-003048D374CA.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/B02E64F7-8895-E111-9E12-5404A63886AE.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/B05E5A37-8795-E111-91AB-BCAEC5364CFB.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/B2647072-8D95-E111-A331-BCAEC532971C.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/B4F068F6-9195-E111-AC0A-0025901D5DB8.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/B6917487-8E95-E111-9744-001D09F2983F.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/B8D7BC01-8995-E111-AEFB-001D09F295A1.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/BAF91739-8C95-E111-BC77-001D09F24DA8.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/BEE10A6D-8D95-E111-AD59-5404A63886AF.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/C09F363C-8795-E111-B919-003048D2C01E.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/C237B9D8-8C95-E111-9011-0025901D62A0.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/C2F123F5-9195-E111-9D5E-0025901D6268.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/C411E2A0-9295-E111-BBB2-BCAEC532971D.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/C83B1DA1-9295-E111-BC09-003048F11112.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/C89C7BFF-8D95-E111-BAC4-003048D2C174.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/CAFAD4A0-9295-E111-91CC-BCAEC5329713.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/CE31A0D7-8C95-E111-A049-BCAEC518FF76.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/D08727F7-8895-E111-8CE2-BCAEC518FF7C.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/D2EFA97B-9095-E111-BC96-00237DDC5CB0.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/D4F68512-8895-E111-B8B4-003048F1BF68.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/D86FE338-8795-E111-AC6F-BCAEC532971E.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/DAFE89FF-8D95-E111-A81A-002481E0D7D8.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/DC05C8D7-8C95-E111-8A85-0025901D625A.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/DC648C7E-9095-E111-905C-003048F024DA.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/DEA9FEF9-8895-E111-A951-003048D2C01E.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/DEB99437-8795-E111-A552-BCAEC518FF40.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/DED8EB39-8C95-E111-9B63-001D09F26C5C.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/E0937CCA-8F95-E111-AE58-003048F024DA.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/E094EDFF-8D95-E111-88FA-BCAEC518FF8A.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/E4062181-9095-E111-A5CA-003048D2BB58.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/E65425CE-8F95-E111-9DEF-003048673374.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/E83BA4D7-8C95-E111-B290-0025901D5DB8.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/EEF876A0-9295-E111-B977-BCAEC53296F4.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/F035843C-8C95-E111-91C6-001D09F24FEC.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/F091AB76-8D95-E111-BC70-BCAEC532971D.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/F25AACCB-8F95-E111-ADD2-003048D2C1C4.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/F433C075-8D95-E111-A73A-5404A640A648.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/F4862BD6-8C95-E111-8EBE-5404A640A63D.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/F4F616F7-8895-E111-85DD-0025901D5DB2.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/FA1AF863-8995-E111-8EE8-BCAEC5329702.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/FA921AF7-8895-E111-BD77-003048F11C5C.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/FADC4784-8E95-E111-A845-002481E0D90C.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/FC0F5100-8E95-E111-8760-001D09F2305C.root',
   '/store/data/Run2012A/LP_MinBias2/RECO/PromptReco-v1/000/193/092/FCD076A1-9295-E111-8EC7-001D09F295FB.root')
   )                               
  #dcap://pnfs/iihe/cms/store/user/xjanssen/data/CMSSW_3_6_2/DataCopy_36x/__MinimumBias__Commissioning10-Jun14thReReco_v1__RECO/DataCopy_36x__CMSSW_3_6_2__MinimumBias__Commissioning10-Jun14thReReco_v1__RECO_1_1_xom.root'),
                            #fileNames = cms.untracked.vstring('dcap:///pnfs/iihe/cms/ph/sc4/store/mc/Summer10/MinBias_7TeV-pythia8/GEN-SIM-RECODEBUG/START36_V10_SP10-v1/0002/1C57F8AA-8E74-DF11-ADCE-0017A4771024.root')
                            #lumisToProcess = cms.untracked.VLuminosityBlockRange('124009:1-124009:68')
                            #)


# configure modules via Global Tag --------------------------------------------------
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_R_36X_V12A::All'
process.GlobalTag.globaltag = 'START50_V16A::All'


#Geometry --------------------------------------------------------------------------
process.load("Configuration.StandardSequences.Geometry_cff")
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
#process.mypath=cms.Path(process.hltPhysicsDeclared+<my_other_things_here>)


from UATree.UABaseTree.UABaseTree_cfi import *
process.load("UATree.UABaseTree.UABaseTree_cfi")

#FileName of the output of the UATree. Needed here for AutoCrab.
process.uabasetree.outputfilename = cms.untracked.string("UABaseTree_MC_Z2star_8TeV_CMS-TOTEM.root")

process.path = cms.Sequence()


#  ----------------------------   Filters   ----------------------------
process.load('UATree.UABaseTree.UABaseTree_filters_cfi')

if not(isMonteCarlo):
   process.path = cms.Sequence(process.path * process.hltPhysicsDeclared)

if useMITFilter:
   process.path = cms.Sequence(process.path * process.mitfilter)
else:
   if not(isMonteCarlo):
     process.path = cms.Sequence(process.path * process.noscraping)


#  ----------------------------   if MC   ----------------------------
if isMonteCarlo:
   process.load('UATree.UABaseTree.UABaseTree_MC_cfi')
   process.myTTRHBuilderWithoutAngle4PixelTriplets.ComputeCoarseLocalPositionFromDisk = True
   process.source.fileNames = cms.untracked.vstring('file:/tmp/katsas/STEP2_RAW2DIGI_L1Reco_RECO_VALIDATION_DQM_421_1_tVN.root')
   
   process.GlobalTag.globaltag = 'START50_V16A::All'
   #process.source.fileNames = cms.untracked.vstring('dcap:///pnfs/iihe/cms/store/user/xjanssen/data//CMSSW_3_9_7/DataCopy_397/__GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6__Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1__GEN-SIM-RECO/DataCopy_397__CMSSW_3_9_7__GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6__Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1__GEN-SIM-RECO_1_1_RJI.root')
   #process.source.fileNames = cms.untracked.vstring('file:/user/selvaggi/step2_RAW2DIGI_L1Reco_RECO_9_1_XqK.root')
   
 
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
