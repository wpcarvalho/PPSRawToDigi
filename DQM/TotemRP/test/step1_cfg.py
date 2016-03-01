import FWCore.ParameterSet.Config as cms

process = cms.Process('RECODQM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# load DQM
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")

# my analyzer
#process.DQMExample_Step1 = cms.EDAnalyzer("DQMExample_Step1",
#    electronCollection       = cms.InputTag("gsfElectrons"),
#    caloJetCollection        = cms.InputTag("ak5CaloJets"),
#    pfMETCollection          = cms.InputTag("pfMet"),
#    conversionsCollection    = cms.InputTag("allConversions"),
#    PVCollection             = cms.InputTag("offlinePrimaryVerticesWithBS"),
#    beamSpotCollection       = cms.InputTag("offlineBeamSpot"),
#
#    TriggerEvent             = cms.InputTag('hltTriggerSummaryAOD','','HLT'),
#    TriggerResults           = cms.InputTag('TriggerResults','','HLT'),
#    #last filter of HLTEle27WP80Sequence
#    TriggerFilter            = cms.InputTag('hltEle27WP80TrackIsoFilter','','HLT'),
#    TriggerPath              = cms.string('HLT_Ele27_WP80_v13'),
#
#    PtThrL1 = cms.untracked.double(30.0),
#    PtThrL2 = cms.untracked.double(10.0),
#    PtThrJet = cms.untracked.double(20.0),
#    PtThrMet = cms.untracked.double(20.0),
#)


process.TotemRPDQMSource = cms.EDAnalyzer("TotemRPDQMSource"
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/j/jkaspar/software/offline/800_pre5/user/reco_test/reco_test.root'
        )
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
                                     fileName = cms.untracked.string("OUT_step1.root"))

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')  #for MC


# Path and EndPath definitions
#process.dqmoffline_step = cms.Path(process.DQMExample_Step1)
process.dqmoffline_step = cms.Path(process.TotemRPDQMSource)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)


# Schedule definition
process.schedule = cms.Schedule(
    process.dqmoffline_step,
    process.DQMoutput_step
)
