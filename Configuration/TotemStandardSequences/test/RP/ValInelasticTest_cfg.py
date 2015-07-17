import FWCore.ParameterSet.Config as cms

process = cms.Process("valRPinelasticBetaXXXEnergyYYYTeV")

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
process.RPHitDists.outputFile = cms.string('file:valRPinelasticBetaXXXEnergyYYYTeV_RPHitDists.root')

# module RPRecTracksVal
#process.load("TotemRPValidation.RPReconstructedTracksValidation.ReconstTracksVal_beta_90_cfi")
#process.RPRecTracksVal.HistogramFileName = 'file:valRPinelasticBetaXXXEnergyYYYTeV_RPRecTracksVal.root'

# module RPAngleVal
#process.load("TotemRPValidation.RPAngleValidation.AngleVal_7000_90_cfi")
#process.RPAngleVal.HistogramFileName = 'file:valRPinelasticBetaXXXEnergyYYYTeV_RPAngleVal.root'

# module RPInelProtRecVal
#process.load("TotemRPValidation.InelasticReconstructionValidation.InelasticReconstVal_90_cfi")
#process.RPInelProtRecVal.HistogramFileName = 'file:valRPinelasticBetaXXXEnergyYYYTeV_RPInelProtRecVal.root'

# module SmearingValidation
process.load("TotemRPValidation.BeamSmearing.SmearingValidation_cfi")
process.SmearingValidation.verbosity = 0
process.SmearingValidation.generatorLabel = cms.string('generator')
process.SmearingValidation.outputFile = cms.string('file:valRPinelasticBeta90Energy7TeV_SmearingValidation.root')

#process.load("TotemRPValidation.RPAngleValidation.AngleVal_3500_2p0_cfi")
#process.RPAngleVal.HistogramFileName = 'file:valRPinelasticBeta2Energy3.5TeV_RPAngleVal.root'



#process.p = cms.Path(process.RPHitDists*process.RPRecTracksVal*process.RPAngleVal*process.RPInelProtRecVal*process.SmearingValidation)

#### end of default configuration ###

process.setName_("valRPinelasticBeta90Energy7TeV")

process.source.fileNames = cms.untracked.vstring('file:prodRPinelasticBeta90Energy7TeV.root')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)     # -1 means to take all the events process.load("The source as input
)

process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90/RP_Dist_Beam_Cent.xml')

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_90_cfi")

process.RPHitDists.outputFile = cms.string('file:valRPinelasticBeta90Energy7TeV_RPHitDists.root')

process.load("TotemRPValidation.RPReconstructedTracksValidation.ReconstTracksVal_beta_90_cfi")
process.RPRecTracksVal.HistogramFileName = cms.string('file:valRPinelasticBeta90Energy7TeV_RPRecTracksVal.root')

process.load("TotemRPValidation.RPAngleValidation.AngleVal_7000_90_cfi")
process.RPAngleVal.HistogramFileName = 'file:valRPinelasticBeta90Energy7TeV_RPAngleVal.root'

process.load("TotemRPValidation.InelasticReconstructionValidation.InelasticReconstVal_90_cfi")
process.RPInelProtRecVal.HistogramFileName = cms.string('file:valRPinelasticBeta90Energy7TeV_RPInelProtRecVal.root')
process.RPInelProtRecVal.Verbosity = 10


process.SmearingValidation.outputFile = cms.string('file:valRPinelasticBeta90Energy7TeV_SmearingValidation.root')



process.RPInelProtRecVal.RPReconstructedProtonCollectionLabel = cms.InputTag("RP220Reconst")



process.p = cms.Path(
                    process.RPHitDists
                    *process.RPRecTracksVal
                    *process.RPAngleVal
                    *process.RPInelProtRecVal
                    *process.SmearingValidation
                     )













