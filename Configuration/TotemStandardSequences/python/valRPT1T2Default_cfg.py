import FWCore.ParameterSet.Config as cms

process = cms.Process("valT1T2RPDefault")

process.load("Configuration.TotemCommon.LoggerMin_cfi")

# Configure the input module (read the events from a file)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:prodT1T2RPDefault.root'),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    skipBadFiles = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)     # -1 means to take all the events from the source as input
)

# random number generators seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# Geometry
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

#process.load("Configuration.TotemCommon.geometryGlobal_cfi")
#process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_1180_Beta_11_220/RP_Dist_Beam_Cent.xml')

#process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_1180GeV_11_cfi")

######################## RP Validation modules 

process.load("TotemRPValidation.HitDistributions.HitDistributions_cfi")
process.RPHitDists.outputFile = cms.string('file:valRPT1T2Default_RPHitDists.root')

process.load("TotemRPValidation.RPReconstructedTracksValidation.ReconstTracksVal_beta_90_6500GeV_cfi")
process.RPRecTracksVal.SmearedHepMCModuleName = cms.string("generator")

process.load("TotemRPValidation.RPAngleValidation.AngleVal_1180_11_cfi")
process.RPAngleVal.SmearedHepMCModuleName = cms.string("generator")

process.load("TotemRPValidation.InelasticReconstructionValidation.InelasticReconstVal_90_6500GeV_cfi")
process.RPInelProtRecVal.SmearedVertexHepMCModuleName = cms.string("generator")

process.load("TotemRPValidation.BeamSmearing.SmearingValidation_cfi")
process.SmearingValidation.verbosity = 0
process.SmearingValidation.outputFile = cms.string('file:valRPT1T2Default_SmearingValidation.root')

######################## T1 Validation modules 

from Configuration.TotemStandardSequences.T1_Validation_cfi import *

process.T1Val = t1valid


######################## T2 Validation modules 

from Configuration.TotemStandardSequences.T2_Validation_cfi import *

process.T2RecoAn = T2RecoAn
process.T2RecoAn.HepMCProductLabel = cms.string("generator")

process.load("TotemAnalysis.TotemNtuplizer.TotemNtuplizer_cfi")
process.TotemNtuplizer.verbosity = 0
process.TotemNtuplizer.ProductLabelSimu = 'RPCCRawBits'
process.TotemNtuplizer.ModulLabelSimu = 'RPDataCCProducer'
process.TotemNtuplizer.T1TrackLabel = 't1tracks2'

#process.p1 = cms.Path(process.T2RecoAn*process.RPHitDists*process.RPRecTracksVal*process.RPAngleVal*process.RPInelProtRecVal*process.SmearingValidation*process.TotemNtuplizer)
