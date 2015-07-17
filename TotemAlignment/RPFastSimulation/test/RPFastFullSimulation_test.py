import FWCore.ParameterSet.Config as cms

process = cms.Process("FastAlignmentProfiles")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# set number of events
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# random seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# declare optics parameters
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_1535_cfi")

# elastic generator
process.load("SimGeneral.HepPDTESSource.pdt_cfi")
process.load("IOMC.Elegent.ElegentSource_cfi")
process.generator.fileName = 'IOMC/Elegent/data/7000GeV_0_20_4E3.root'
process.generator.verbosity = 1
process.generator.t_min = 0
process.generator.model = "PH/ppp3"

# smearing
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")
process.SmearingGenerator.verbosity = 2

# ideal geometry
process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/' + 'RP_Beta_1535_150_out' + '/RP_Dist_Beam_Cent.xml')

# misalignments
process.load("Geometry.TotemRPGeometryBuilder.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.MisalignedFiles = cms.vstring()

# fast simulation
process.load("TotemAlignment.RPFastSimulation.RPFastFullSimulation_cfi")
process.RPFastFullSimulation.verbosity = 1

# one-RP track fitting
process.load("RecoTotemRP.RPTrackCandidateCollectionFitter.RPSingleTrackCandCollFitted_cfi")
process.RPSingleTrackCandCollFit.Verbosity = 0

# alignment profiles
process.load("TotemAlignment.RPProfiles.RPAlignmentProfiles_cfi")
process.RPAlignmentProfiles.verbosity = 0

process.p = cms.Path(
    process.generator *
    process.SmearingGenerator *
    process.RPFastFullSimulation
)
