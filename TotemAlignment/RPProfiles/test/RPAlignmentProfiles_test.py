import FWCore.ParameterSet.Config as cms

process = cms.Process("FastAlignmentProfiles")

# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# set number of events
process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

# random seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# station simulation
process.load("TotemAlignment.RPFastSimulation.RPFastStationSimulation_cfi")

# ideal geometry
process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/' + 'RP_Beta_1535_150_out' + '/RP_Dist_Beam_Cent.xml')

# misalignments
process.load("Geometry.TotemRPGeometryBuilder.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.MisalignedFiles = cms.vstring()

# one-RP track fitting
process.load("RecoTotemRP.RPTrackCandidateCollectionFitter.RPSingleTrackCandCollFitted_cfi")
process.RPSingleTrackCandCollFit.Verbosity = 0

# alignment profiles
process.load("TotemAlignment.RPProfiles.RPAlignmentProfiles_cfi")
process.RPAlignmentProfiles.verbosity = 1

process.p = cms.Path(
    process.RPFastStationSimulation *
    process.RPSingleTrackCandCollFit *
    process.RPAlignmentProfiles
)
