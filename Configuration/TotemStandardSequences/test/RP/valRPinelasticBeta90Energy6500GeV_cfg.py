import FWCore.ParameterSet.Config as cms

from Configuration.TotemStandardSequences.valRPinelasticDefault_cfg import *

process.setName_("valRPinelasticBeta90Energy6500GeV")

process.source.fileNames = cms.untracked.vstring('file:prodRPinelasticBeta90Energy6500GeV.root')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)     # -1 means to take all the events process.load("The source as input
)

process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90/RP_Dist_Beam_Cent.xml')

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_6500GeV_90_cfi")

process.load("TotemRPValidation.RPReconstructedTracksValidation.ReconstTracksVal_beta_90_6500GeV_cfi")
process.load("TotemRPValidation.RPAngleValidation.AngleVal_6500_90_cfi")
process.load("TotemRPValidation.InelasticReconstructionValidation.InelasticReconstVal_90_6500GeV_cfi")

process.RPHitDists.outputFile = cms.string('file:'+process.name_()+'_RPHitDists.root')
process.RPRecTracksVal.HistogramFileName = cms.string('file:'+process.name_()+'_RPRecTracksVal.root')
process.RPAngleVal.HistogramFileName = 'file:'+process.name_()+'_RPAngleVal.root'
process.RPInelProtRecVal.HistogramFileName = cms.string('file:'+process.name_()+'_RPInelProtRecVal.root')
process.SmearingValidation.outputFile = cms.string('file:'+process.name_()+'_SmearingValidation.root')
process.TotemNtuplizer.outputFileName = cms.untracked.string('file:'+process.name_()+'_ntuple.root')

process.p = cms.Path(process.RPHitDists*process.RPRecTracksVal*process.RPAngleVal*process.RPInelProtRecVal*process.SmearingValidation*process.TotemNtuplizer)
