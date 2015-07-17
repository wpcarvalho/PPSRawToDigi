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

# module RPRecTracksVal
#process.load("TotemRPValidation.RPReconstructedTracksValidation.ReconstTracksVal_beta_90_cfi")
#process.RPRecTracksVal.HistogramFileName = 'file:valRPelasticBeta90Energy7TeV_RPRecTracksVal.root'

# module RPAngleVal
#process.load("TotemRPValidation.RPAngleValidation.AngleVal_7000_90_cfi")
#process.RPAngleVal.HistogramFileName = 'file:valRPelasticBeta90Energy7TeV_RPAngleVal.root'

# module ElasticRecoVal
process.load("TotemRPValidation.ElasticReconstruction.ElasticRecoVal_cfi")
process.ElasticRecoVal.verbosity = 0
process.ElasticRecoVal.validationOutputFile = cms.string('file:valRPelasticBeta90Energy7TeV_ElasticRecoVal_validation.root')
process.ElasticRecoVal.acceptanceOutputFile = cms.string('file:valRPelasticBeta90Energy7TeV_ElasticRecoVal_acceptance.root')

# module SmearingValidation
process.load("TotemRPValidation.BeamSmearing.SmearingValidation_cfi")
process.SmearingValidation.verbosity = 0
process.SmearingValidation.outputFile = cms.string('file:valRPelasticBeta90Energy7TeV_SmearingValidation.root')

process.load("TotemAnalysis.TotemNtuplizer.TotemNtuplizer_cfi")
process.TotemNtuplizer.outputFileName = 'file:valRPelastic_ntuple.root'
process.TotemNtuplizer.verbosity = 0
process.TotemNtuplizer.ProductLabelSimu = 'RPCCRawBits'
process.TotemNtuplizer.ModulLabelSimu = 'RPDataCCProducer'

#process.load("TotemRPValidation.RPAngleValidation.AngleVal_3500_2p0_cfi")
#process.RPAngleVal.HistogramFileName = 'file:valRPelasticBeta2Energy3.5TeV_RPAngleVal.root'

#process.p = cms.Path(process.RPHitDists*process.RPRecTracksVal*process.RPAngleVal*process.ElasticRecoVal*process.SmearingValidation)
