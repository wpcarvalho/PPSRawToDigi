import FWCore.ParameterSet.Config as cms

process = cms.Process("RealGeometryPositionValidation")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.source = cms.Source("EmptySource")


# minimum of logs
process.load("Configuration.TotemCommon.LoggerMax_cfi")

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_3500GeV_3p5_100urad_cfi")


#process.TotemRPIncludeAlignments = cms.ESProducer("TotemRPIncludeAlignments",
#    RealFiles = cms.vstring('Geometry/TotemRPAlignmentData/data/tb_profiles_45_56_2010_07_15.xml'),
#    MisalignedFiles = cms.vstring()
#)

# geometry
process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/2010_07_15/RP_Dist_Beam_Cent.xml")
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

process.load("TotemRPValidation.RPPositionValidation.RPPositionValidation_cfi")
process.RPPosVal.Verbosity = cms.untracked.uint32(10)

process.p = cms.Path(
     process.RPPosVal
)
