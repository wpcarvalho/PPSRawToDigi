import FWCore.ParameterSet.Config as cms

process = cms.Process("EdgeEfficiency")

process.maxEvents = cms.untracked.PSet(
		input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
		fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/user/f/fnemes/CMSSW_3_1_1_analysis/ver1.5/3720//run_3720.000-1.root"),
		noEventSort = cms.untracked.bool(True),
		skipBadFiles = cms.untracked.bool(True),
		duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_3500GeV_3p5_100urad_cfi")

# raw to digi conversion
process.load('TotemCondFormats.DAQInformation.DAQMappingSourceXML_cfi')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_220.xml')

process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/2010_10_29_temp/RP_Dist_Beam_Cent.xml")
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

process.load("TotemRPValidation.EdgeEfficiency.EdgeEfficiency_cfi")
process.EdgeEfficiency.Verbosity = 0 
process.EdgeEfficiency.selectedTest = 2

process.p = cms.Path(
		process.EdgeEfficiency 
)

