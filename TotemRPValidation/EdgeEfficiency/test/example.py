import FWCore.ParameterSet.Config as cms

process = cms.Process("valRPinelasticBeta3.5Energy3.5TeV220100urad")

########################### COMMON PART ##########################################

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-1.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-2.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-3.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-4.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-5.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-6.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-7.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-8.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-9.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-10.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-11.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-12.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-13.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-14.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-15.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-16.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-17.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-18.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-19.root',
	'file:./input/RPinelasticBeta3.5Energy3.5TeV220100urad-20.root'),
    noEventSort = cms.untracked.bool(True),
    skipBadFiles = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_3500_Beta_3.5_100urad_220/RP_Dist_Beam_Cent.xml')

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_3500GeV_3p5_100urad_cfi")

# logging to txt files 
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# random number generators seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

########################### RP VALIDATION ##########################################

process.load("TotemRPValidation.EdgeEfficiency.EdgeEfficiency_cfi")

# TODO
process.p1 = cms.Path(process.EdgeEfficiency)




