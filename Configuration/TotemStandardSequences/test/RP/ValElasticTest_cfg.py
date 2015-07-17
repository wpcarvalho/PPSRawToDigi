import FWCore.ParameterSet.Config as cms

process = cms.Process("valRPelasticBetaXXXEnergyYYYTeV")

# logging to txt files
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# Configure the input module (read the events from a file)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:prodRPelasticBetaXXXEnergyYYYTeV.root'),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    skipBadFiles = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)     # -1 means to take all the events from the source as input
)

process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

#process.load("Configuration.TotemCommon.geometryRP_cfi")
#process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90_150_out/RP_Dist_Beam_Cent.xml')

#process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_90_cfi")

# module RPHitDists
process.load("TotemRPValidation.HitDistributions.HitDistributions_cfi")
process.RPHitDists.outputFile = cms.string('file:valRPelasticBeta90Energy7TeV_RPHitDists.root')
process.RPHitDists.verbosity = 0

# module RPRecTracksVal
#process.load("TotemRPValidation.RPReconstructedTracksValidation.ReconstTracksVal_beta_90_cfi")
#process.RPRecTracksVal.HistogramFileName = 'file:valRPelasticBeta90Energy7TeV_RPRecTracksVal.root'

# module RPAngleVal
#process.load("TotemRPValidation.RPAngleValidation.AngleVal_7000_90_cfi")
#process.RPAngleVal.HistogramFileName = 'file:valRPelasticBeta90Energy7TeV_RPAngleVal.root'

# module ElasticRecoVal

process.load("TotemRPValidation.ElasticReconstruction.ElasticRecoVal_cfi")
process.ElasticRecoVal.verbosity = 0
# process.ElasticRecoVal.validationOutputFile = cms.string('file:valRPelasticBeta90Energy7TeV_ElasticRecoVal_validation.root')
# process.ElasticRecoVal.acceptanceOutputFile = cms.string('file:valRPelasticBeta90Energy7TeV_ElasticRecoVal_acceptance.root')


# module SmearingValidation
process.load("TotemRPValidation.BeamSmearing.SmearingValidation_cfi")
process.SmearingValidation.verbosity = 0
process.SmearingValidation.outputFile = cms.string('file:valRPelasticBeta90Energy7TeV_SmearingValidation.root')


process.load("TotemAnalysis/TotemNtuplizer/TotemNtuplizer_cfi")
process.TotemNtuplizer.outputFileName = 'file:valRPelastic_ntuple.root'
process.TotemNtuplizer.verbosity = 10
process.TotemNtuplizer.ProductLabelSimu = 'RPCCRawBits'
process.TotemNtuplizer.ModulLabelSimu = 'RPDataCCProducer'

'''
#process.load("TotemRPValidation.RPAngleValidation.AngleVal_3500_2p0_cfi")
#process.RPAngleVal.HistogramFileName = 'file:valRPelasticBeta2Energy3.5TeV_RPAngleVal.root'

#process.p = cms.Path(process.RPHitDists*process.RPRecTracksVal*process.RPAngleVal*process.ElasticRecoVal*process.SmearingValidation)
process.load("TotemAnalysis/TotemNtuplizer/TotemNtuplizer_cfi")
process.TotemNtuplizer.outputFileName = 'file:valRPT1T2Default_ntuple.root'
process.TotemNtuplizer.verbosity = 0
process.TotemNtuplizer.ProductLabelSimu = 'RPCCRawBits'
process.TotemNtuplizer.ModulLabelSimu = 'RPDataCCProducer'
'''


#########end of default configuration ################

process.setName_("valRPelasticBeta90Energy7TeV")

process.source.fileNames = cms.untracked.vstring('file:prodRPelasticBeta90Energy7TeV.root')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)     # -1 means to take all the events from the source as input
)

process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90_150_out/RP_Dist_Beam_Cent.xml')

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_90_cfi")

process.RPHitDists.outputFile = cms.string('file:valRPelasticBeta90Energy7TeV_RPHitDists.root')

process.load("TotemRPValidation.RPReconstructedTracksValidation.ReconstTracksVal_beta_90_cfi")
process.RPRecTracksVal.HistogramFileName = 'file:valRPelasticBeta90Energy7TeV_RPRecTracksVal.root'
# process.RPRecTracksVal.Verbosity= 10



process.load("TotemRPValidation.RPAngleValidation.AngleVal_7000_90_cfi")
process.RPAngleVal.SmearedHepMCModuleName = cms.string("generator")
process.RPAngleVal.HistogramFileName = 'file:valRPelasticBeta90Energy7TeV_RPAngleVal.root'
# process.RPAngleVal.Verbosity = 10


process.ElasticRecoVal.validationOutputFile = cms.string('file:valRPelasticBeta90Energy7TeV_ElasticRecoVal_validation.root')
process.ElasticRecoVal.acceptanceOutputFile = cms.string('file:valRPelasticBeta90Energy7TeV_ElasticRecoVal_acceptance.root')

process.SmearingValidation.outputFile = cms.string('file:valRPelasticBeta90Energy7TeV_SmearingValidation.root')

# process.TotemNtuplizer.outputFile = cms.string('file:valRPelasticBeta90Energy7TeV_ntuple.root')


process.RPHitDists.RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit")

process.RPRecTracksVal.RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit")
process.RPRecTracksVal.RPDigClusterSetLabel = cms.InputTag("RPClustProd")
process.RPRecTracksVal.RPDetTriggerSetLabel = cms.InputTag("RPSiDetDigitizer")

process.RPAngleVal.RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit")

process.ElasticRecoVal.RPDetTriggerSetLabel = cms.InputTag("RPSiDetDigitizer")
process.ElasticRecoVal.RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit")
process.ElasticRecoVal.RPRecoElasticEventLabel = cms.InputTag("ElasticReconstruction")

process.TotemNtuplizer.RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit")
process.TotemNtuplizer.RPMulFittedTrackCollectionLabel = cms.InputTag("RPMulTrackCandCollFit")
process.TotemNtuplizer.RPStripDigiSetLabel = cms.InputTag("RPSiDetDigitizer")
process.TotemNtuplizer.RPDigClusterLabel = cms.InputTag("RPClustProd")
process.TotemNtuplizer.RPReconstructedProtonCollectionLabel = cms.InputTag("ElasticReconstruction")
process.TotemNtuplizer.RPReconstructedProtonPairCollectionLabel = cms.InputTag("ElasticReconstruction")
process.TotemNtuplizer.T2PadDigiCollectionLabel = cms.InputTag("T2Digis")
process.TotemNtuplizer.T2StripDigiCollectionLabel = cms.InputTag("T2Digis")
process.TotemNtuplizer.T1DigiWireCollectionLabel = cms.InputTag("T1Digis")
process.TotemNtuplizer.T1DigiVfatCollectionLabel= cms.InputTag("T1Digis")
process.TotemNtuplizer.T1RecHit2DCollectionLabel = cms.InputTag("T1Digis")
process.TotemNtuplizer.RawEventLabel = cms.InputTag("RPSiDetDigitizer")

process.p = cms.Path(process.RPHitDists
                    *process.RPRecTracksVal
                    *process.RPAngleVal
                    *process.ElasticRecoVal
                    *process.SmearingValidation
                    *process.TotemNtuplizer
                     )
