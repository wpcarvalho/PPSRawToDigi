import FWCore.ParameterSet.Config as cms

process = cms.Process("ModifySingularModes")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# empty source
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

input = './../sr+hsx/45_220.xml'

# geometry
#process.load("Configuration.TotemCommon.geometryRP_real_cfi")
process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/2011_10_20_1/RP_Dist_Beam_Cent.xml")
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

# include alignments, if any
process.load("Geometry.TotemRPGeometryBuilder.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.RealFiles = cms.vstring(input)

# worker
process.load("TotemAlignment.RPTrackBased.ModifySingularModes_cfi")
process.ModifySingularModes.inputFile = input
process.ModifySingularModes.outputFile = './45_220.xml'

# far unit
process.ModifySingularModes.z2 = -220000
process.ModifySingularModes.de_rho2 = -0.05E-3
process.ModifySingularModes.de_x2 = +7.2E-3
process.ModifySingularModes.de_y2 = +26E-3

# near unit
process.ModifySingularModes.z1 = -214628
process.ModifySingularModes.de_rho1 = -1.69E-3
process.ModifySingularModes.de_x1 = +8.2E-3
process.ModifySingularModes.de_y1 = +97E-3

process.p = cms.Path(process.ModifySingularModes)
