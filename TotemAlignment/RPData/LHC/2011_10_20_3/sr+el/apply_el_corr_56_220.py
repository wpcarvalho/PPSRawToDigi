import FWCore.ParameterSet.Config as cms

process = cms.Process("ModifySingularModes")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# empty source
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

input = './../sr+hsx/56_220.xml'

# geometry
#process.load("Configuration.TotemCommon.geometryRP_real_cfi")
process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/2011_10_20_3/RP_Dist_Beam_Cent.xml")
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

# include alignments, if any
process.load("Geometry.TotemRPGeometryBuilder.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.RealFiles = cms.vstring(input)

# worker
process.load("TotemAlignment.RPTrackBased.ModifySingularModes_cfi")
process.ModifySingularModes.inputFile = input
process.ModifySingularModes.outputFile = './56_220.xml'

# near unit
process.ModifySingularModes.z1 = +214628
process.ModifySingularModes.de_rho1 = 2E-3
process.ModifySingularModes.de_x1 = +15.7E-3
process.ModifySingularModes.de_y1 = 22E-3

# far unit
process.ModifySingularModes.z2 = +220000
process.ModifySingularModes.de_rho2 = -0.8E-3
process.ModifySingularModes.de_x2 = +10.8E-3
process.ModifySingularModes.de_y2 = 49E-3

process.p = cms.Path(process.ModifySingularModes)
